# AED light & extinction standard — revised plan

Status: **DRAFT for review** (2026-09-27). Supersedes `LIGHT_PLAN.md`, which is
kept for history. This revision folds in a holistic review of the first plan:
three claims corrected against the code (D14 latent, `extc_coef` host-dependent,
the "local frame" not normalising orientation), gaps closed (oasim's role, Q1
availability, surface shading feedback, GLM heating bands), and the staging
collapsed into a single schedule (§10).

Owner of the standard: `libaed-api`. Band registry and spectra: `libaed-water`.
Atmosphere and spectral light models: `libaed-light`.

---

## 1. Purpose and scope

The AED API (`aed_api.F90`) couples libaed to four hosts — **GLM, TUFLOW-FV,
SCHISM, ELCOM**. The API defines what each environmental quantity means; each
host interface (`glm_api_aed.F90`, `fv_aed.F90`, `schism_aed.F90`,
`elcom_aed.F90`) converts host data to conform. Light is the least well-defined
part of that contract today.

This plan:

1. defines every light quantity once (surface, in-water, extinction, feedback);
2. makes light **multi-band**, with band-specific bio-shading feedback to hosts
   that carry bands (ELCOM, TUFLOW-FV, GLM heating bands);
3. defines the roles of the three light components — the API light manager,
   `aed_atmosphere` and `aed_oasim` — so there is one light field biology uses;
4. restructures `libaed-light` around an `aed_light_utils` collection;
5. stages the work in git so master only receives verified checkpoints, and
   every behaviour change is registered with a continuity switch.

Out of scope: host hydrodynamic heating schemes themselves (only what AED
returns to them), and calibration changes to module coefficients.

## 2. Status at 2026-09-30

| Item | State |
|---|---|
| Run `_12` segfault (D1) | **Fixed.** Uncommitted Sep-24 edit (`Light` used `idata%rad(top)`, null in TUFLOW-FV) reverted; `Light` uses `I_0` |
| A1 — libaed-api Stage 0 comment + `LIGHT_PLAN.md` | Committed (`83c7b81`) |
| A2 — retire libaed-light `aed_light.F90` | Committed; GLM (both drivers) identical in Docker |
| A3 — `fv_aed` `rad_` → `irr_host_` | Edited, uncommitted, awaiting VM build |
| A4 — `output=` flag on diagnostics | Committed; GLM identical |
| A5 — D14 accumulate (macrophyte, vegetation) | Committed; GLM identical, also with the Q7 binding enabled |
| A7 — drop `slingo_const.h` from vfproj | Committed with A2 |
| A8 — this plan | Committed |
| A9 — D4/D15 check | **TUFLOW-FV confirmed** on Conjimba `_11` (§3.5). ELCOM and the baselines still to do |
| B1 — light tidy | Committed on `standardize_light` (libaed-water, libaed-light). GLM suite (7 cases: CGM, macrophyte ×2 drivers, Rotorua oasim `sky_model` 1 and 2 ×2 drivers) bit-identical to baseline. Needs libaed-water linked again after the plus libs (static link order): fixed in GLM, zeroD, ELCOM, libaed-fv Makefiles; schism-aed CMake still to do |
| Everything else | Not started |

Unrelated local changes to keep out of light commits:
`libaed-benthic/src/aed_bivalve.F90`, `schism-aed/build_local.sh`.

---

## 3. Current-state audit

### 3.1 Light quantities

| # | Quantity | Units | Used for |
|---|----------|-------|----------|
| **Q1** | Incident shortwave *above* the surface (before reflection) | W/m² | oasim (applies its own reflectance); dry cells, canopy/riparian modules; atmosphere scaling |
| **Q2** | Net shortwave *just below* the surface (after albedo, and ice/snow where modelled) | W/m² | top of the in-water profile; today's `par_sf` consumers |
| **Q3** | In-water irradiance per layer, per band | W/m² | photosynthesis (`par`), UV processes (`uva`, `uvb`), `nir` |
| **Q4** | Total net shortwave per layer | W/m² | informational only (no module uses API `rad`) |
| **Q5** | Kd per layer, per band — background + AED contribution | m⁻¹ | in-water light; layer-integrated photosynthesis; host heating via feedback |
| **Q6** | Surface shading fraction from AED (canopy/riparian) | – | host reduces incoming SW (`sw_factor`) |

Today's names conflate them: `I_0`/`par_sf` (called PAR, used as total SW,
Q1 or Q2 depending on host), `rad` (Q4 in GLM, band array in TUFLOW-FV, zeros
in SCHISM, band sum in ELCOM), `extc`/`extc_coef` (total Kd in GLM, bio-only
`bioshade` in TUFLOW-FV, band 1 in ELCOM).

### 3.2 Host inventory

| | GLM | TUFLOW-FV | SCHISM | ELCOM |
|---|---|---|---|---|
| Interface | `GLM/src/glm_api_aed.F90` | `libaed-fv/src/fv_aed.F90` (+ `tuflowfv_external_wq/`) | `schism-aed/src/AED/schism_aed.F90` | `ELCOM/elcom_aed/src/elcom_aed.F90` |
| `I_0` bound to | `MetData%ShortWave` (:268) = **Q1 gross** | `wq%I_0`, "NET SURFACE SHORTWAVE" (`tuflowfv_wq_api.f90:129`) = **Q2** | `srad` (:522) — Q2?, **verify** | met `SOLAR_RAD_INDEX` (`elcom_wq.f90:279`) — **verify** |
| Q1 available? | yes (`MetData%ShortWave`) | not in WQ API 1.0 | **verify** | raw met, **verify** |
| `rad` bound to | `theLake%Light` (:282) = Q4 | not bound | zeros (:199, :679) | Σ bands (:636, :809) |
| Host band irradiance | — | `wq%par(nband,NC3)`, band 1 used | none | `wq%par/nir/uva/uvb` (`elcom_wq.f90:410-413`) |
| Host band Kd | `n_bands`, `light_extc` (heating) | per band (host) | none | `wq%extc(col,lev,4)` |
| `extc` (→ `extc_coef`) | `ExtcCoefSW` (total SW Kd) | `extcoeff_` = `bioshade` (bio-only) | – | `extc(:,:,1)` |
| Bio feedback | `ExtcCoefSW` | `wq%bioshade(NC3)` | none | `wq%bioextc` added to all 4 bands (`elcom_wq.f90:672-676` TODO) |
| `link_ext_par` | nml | nml | nml (F) | **forced F** (:366) |
| `par_fraction` | 0.43 | nml 0.43 | nml 0.43 | **hard-coded 0.45** (:393) |
| Host layer order | bottom-up | top-down | bottom-up | top-down |
| `cpl%Kw` | pointer (:219) | pointer (:317) | pointer (:441) | pointer (:415) — all constant after init |

### 3.3 Current API light flow (`aed_run_column`, `aed_api.F90:1616-1703`)

1. `aed_run_model` maps each column to 1..`col_lev` (`:1386-1397`). **It
   renumbers but keeps host orientation**: a top-down host gets `top=1`,
   `bot=col_lev`.
2. `Light` (unless `link_ext_par` or particles): `BioExtinction` sums each
   module's scalar `light_extinction`; `extc = bio + Kw`; writes
   `idata%extc(bot:top-1)` into the host array; `par(top) = par_fraction·I_0`;
   attenuates downward.
3. `nir/uva/uvb = par/par_fraction · frac` — overwrites host arrays.
4. `calculate_fluxes` (modules read `par_sf`, `par`, `extc_coef`, …).
5. After substeps: `BioExtinction → bioextc`; if not `link_ext_par`,
   `extc(bot:top-1) = bioextc + Kw`; zero `bioextc` if not `bioshade_feedback`.
6. Zone/riparian path: `aed_light_shading` → `sw_factor` if
   `link_solar_shade` (`:1953`).

The `bot:top-1` rule is GLM's `update_light` behaviour embedded in shared code.

### 3.4 Defects

| ID | Severity | Defect | Resolution |
|----|----------|--------|-----------|
| **D1** | Blocker | `Light` used `rad(top)`, null in TUFLOW-FV → segfault | **Fixed** (Stage 0) |
| **D2** | High | GLM binds Q1 (gross) to `I_0`; API path over-lights ~7% | host fix, §7 GLM |
| **D3** | High, latent | SCHISM binds zeros to `rad` | host fix, §7 |
| **D4** | High — **confirmed for 1-layer TUFLOW-FV** (§3.5); top-down 3D and ELCOM to verify | `bot:top-1` is empty for top-down hosts (`nlev:0`) and for 1-layer columns (`1:0`): feedback and `extc` write-back never reach TUFLOW-FV/ELCOM or 2D columns | orientation normalised in `link_column` (§5.4) |
| **D5** | High, latent | naive D4 fix would write `Kw+bio` into TUFLOW-FV `bioshade` → double background | bio-only feedback (§5.8) |
| **D6** | Medium | API overwrites host `nir/uva/uvb` and derives them by scaling PAR | band views (§5.4) |
| **D7** | Medium | ELCOM computes band light but `link_ext_par` is forced off → two light fields | host mode (§5.3) |
| **D8** | Medium | `par_sf` is misnamed (it is total SW) | rename to `sw_sf` (§5.1) |
| **D9** | Medium | `aed_phytoplankton.F90:1586` `par > 0.95*Io` compares PAR with total SW; never true | module fix (§8) |
| **D10** | Medium | oasim applies reflection to *net* SW on hosts without Q1 → reflected twice | `sw_inc_supplied` rule (§5.1) |
| **D11** | Low | `fv_aed` `rad_` naming | **A3** |
| **D12** | Low | ELCOM ignores nml `par_frac` | **A6** |
| **D13** | Low | ELCOM `I_0` may be raw met (gross) | host verify (§7) |
| **D14** | Low, **latent** | macrophyte (`:1551/1558`) and vegetation (`:926`) *assign* to the extinction accumulator. Their bindings are commented out (`aed_macrophyte.F90:187`, disabled in `37797cd` "csiem 1.5 temp updates", 2025-10-22; `aed_vegetation.F90:150` never bound; `aed_habitat.F90:98` likewise, and its body is commented out), so the code never runs today | **A5**: accumulate; bit-identical. Re-enabling is a separate decision (Q7) |
| **D15** | **High — confirmed** (§3.5) | `extc_coef` points at the host's `extc` (`aed_api.F90:1181`) and so means total SW Kd (GLM), bio-only `bioshade` (TUFLOW-FV; never written by the API because of D4) or band-1 Kd (ELCOM). Consumers: phytoplankton, carbon, methane, organic_matter, habitat_water, totals. Also only provided `IF (BSSOCIATED(extc))` | API-owned PAR Kd (§5.5); behaviour change C4 |
| **D16** | Low, **new** | global availability is inferred from `BSSOCIATED()` (`aed_api.F90:939-955`, `:1047-1063`), not declared | declarations (§5.11) |

### 3.5 A9 evidence: Conjimba `_11` (TUFLOW-FV, 2D)

`ass_model/output/ERA_S4_097_WQ_EX_AED_11*.nc`: 32,935 cells, one layer each
(`NumCells2D = NumCells3D`), 625 two-hourly outputs.

| Output | Result | Meaning |
|---|---|---|
| `WQ_DIAG_TOT_EXTC` (`aed_totals` copy of `extc_coef`) | **0 in every cell at every output** | D15: `extc_coef` is TUFLOW-FV's `extcoeff_`, which nothing fills; D4: the API write-back `extc(1:0)` is empty for 1-layer columns |
| `WQ_DIAG_TOT_PAR` / `TOT_LIGHT` | nonzero in 310/625 outputs (daytime), max 590 / 1371 W/m², ratio 0.43 | surface light path is fine (`par = par_fraction·I_0`) |
| `WQ_DIAG_PHY_TCHLA`, `TSS` | nonzero, up to 14 µg/L and 308 mg/L | bio/particle Kd *should* be nonzero, so shading feedback is also lost (D4) |

**Consequence for biology.** Conjimba phytoplankton use `lightModel = 0`
(`aed_phyto_pars.csv`: `I_K = 150`, `KePHY = 0.0005`). In
`photosynthesis_irradiance` (`aed_bio_utils.F90:570`) `extc = 0` gives
`par_b = par_t`, so `z1 = z2` and `fI = 1 + 0/1e-3 = 1`: **phytoplankton are
never light-limited in daylight**, whatever the depth, turbidity or
irradiance (apart from the `par_t < 5e-5` cut-off). Other light models see
surface light throughout the layer, and `lightModel = 7` evaluates `0/0`
(NaN). Because every top-down TUFLOW-FV column also has an empty write-back
range, this probably applies to **every TUFLOW-FV + AED API run** unless
TUFLOW-FV itself fills `extcoeff_` (ask BMT, Q6). Any calibration done on
TUFLOW-FV through the API has been done with light limitation effectively off
by day; C4 is therefore a large change for those models.

### 3.6 libaed-light defects (review of `src/`)

| ID | Defect |
|----|--------|
| **LR1** | `vH_O3` defined twice with different units: `aed_slingo` (used by oasim) returns kg m⁻²; `aed_atmosphere` copy returns matm-cm (missing `*0.4462e-3*48/1000`) |
| **LR2** | ~100-line atmospheric block duplicated between `aed_calculate_column_oasim` (:741-850) and `aed_calculate_dry_oasim` (:1163-1266) |
| **LR3** | `oasim.inc` is `#include`d into both `aed_oasim` and `curtin_light` → every table twice, incl. non-`PARAMETER` `coeff_curtin(5,2,140,9)` |
| **LR4** | `deg2rad` in three files; `Planck`/`lightspeed`/`Avogadro` local to oasim |
| **LR5** | `aed_calculate_surface_oasim` is an empty placeholder for the atmosphere |
| **LR6** | `win/libaed-lighting.vfproj` references nonexistent `src\slingo_const.h` |

---

## 4. Design principles

1. **The API defines meaning; hosts conform.** One physical definition, unit
   and sign per field. Albedo, ice, unit changes, layer order and band mapping
   are host-interface work.
2. **Declared, never inferred.** Every choice (light mode, bands, which
   surface quantities exist, feedback targets) is declared at init and checked
   by `validate()`. No behaviour depends on `ASSOCIATED()` of an optional
   pointer at run time.
3. **One light field.** Biology reads exactly one in-water irradiance field per
   band, whose producer is declared by `light_mode` (§5.3).
4. **AED returns only its own contribution.** Feedback is bio/particle Kd and
   surface shading. Background water Kd belongs to whoever computes the light.
5. **Orientation-agnostic.** The API works surface→bottom internally; hosts
   are mapped at the boundary; no index-range arithmetic.
6. **Calibrations preserved.** Module coefficients (`KePHY`, `KeDOM`, …) stay
   the PAR-band specific Kd; spectral information supplies only band ratios.
7. **Checkpoint identity + change register.** Every checkpoint is
   bit-identical to the previous one *except* for entries in the change
   register (§9). Each entry states its effect, the hosts affected, and the
   switch that reproduces old results where continuity runs need it.
8. **One name per quantity.** Renames go through a deprecation alias with a
   log line, and in-tree code moves to the new name in the same checkpoint.

---

## 5. The standard

### 5.1 Surface quantities

| API field | Module global | Definition | Declared by host |
|---|---|---|---|
| `I_0` | **`sw_sf`** (`par_sf` = deprecated alias) | **Q2** net SW just below the surface, W/m² | required |
| `I_inc` | **`sw_inc`** | **Q1** incident SW above the surface, W/m² | optional; `sw_inc_supplied = .TRUE./.FALSE.` |
| `rad` | `rad` | **Q4** net SW per layer, informational | optional; never used to compute light |

**Q1 when the host has none.** If `sw_inc_supplied = .FALSE.`, `sw_inc`
points at `I_0` and the flag is published as a global
(`sw_inc_is_net`). Consumers that model reflection (oasim, atmosphere
`sky_model=2` scaling) or dry-cell exposure read the flag in `define` and,
when set, **skip their own reflection** and treat the input as net. This
replaces the old "fall back with a warning" rule, which left D10 unfixed on
TUFLOW-FV/ELCOM.

**Rename.** `sw_sf` is the only provided global. `aed_locate_sheet_global`
maps `par_sf` → `sw_sf` and logs a one-time deprecation per module. All
in-tree consumers (~20) move to `sw_sf` in the same checkpoint (B4).
`photosynthesis_irradiance`'s `Io` argument is documented as Q2.

### 5.2 Band registry — owned by libaed-water

Modules read bands in `define` and can only `USE` libaed-water, so the registry
lives there, in a new `libaed-water/src/aed_bands.F90`. It is the **only**
definition; `aed_coupling_t` carries no band fields.

```fortran
TYPE aed_light_band_t
   CHARACTER(len=8) :: name   = 'PAR'   ! canonical: PAR, NIR, UVA, UVB; or host-defined
   AED_REAL :: lambda_lo      = 400.    ! nm
   AED_REAL :: lambda_hi      = 700.    ! nm
   AED_REAL :: sw_fraction    = 0.43    ! share of I_0 (api mode)
   AED_REAL :: Kw             = 0.1     ! background Kd, m-1 (api mode); copied at init
END TYPE

! aed_bands: registry + read-only accessors
INTEGER :: n_bands = 1, par_band = 1
PUBLIC aed_n_bands, aed_band, aed_par_band, aed_band_index(name)
PUBLIC aed_bands_set      ! called only by the API light manager
```

- **Filled** by the API `aed_light` manager inside `aed_set_coupling`, which
  all four hosts call before `aed_configure_models` (`fv_aed.F90:323/325`,
  `glm_api_aed.F90:225/229`, `schism_aed.F90:458/462`, `elcom_aed.F90:434/437`).
  The registry is therefore complete before any module `define`. The legacy
  `glm_aed.F90` driver gets the default single PAR band.
- **Frozen** after `validate()`; `aed_bands_set` after that is a fatal error.
- Exactly one PAR band (by name, or covering 400–700 nm).
- `Kw` is **copied** at init. All four hosts' `cpl%Kw` targets are constant
  after init; if a host later needs time-varying background Kd it becomes a
  per-band pointer — noted, not built.
- No `bio_scale`: modules that do not declare spectral behaviour are grey.

### 5.3 Light modes — who produces the in-water field

| `light_mode` | Producer of `irr(b)` | Host supplies | Intended hosts |
|---|---|---|---|
| `'api'` (default) | API: `irr_b(surf) = sw_fraction_b·I_0`, attenuate with `Kw_b + bioKd_b` | `I_0` | GLM, SCHISM, TUFLOW-FV without band light |
| `'host'` | host | `I_0`, `irr(b)` per layer, **and PAR Kd per layer** (declared) | ELCOM, TUFLOW-FV with band light |
| `'module'` (**reserved**) | a registered light module (e.g. oasim) fills `irr(b)` from its own spectral calculation | `I_0` (+ `I_inc`) | later |

The mode is chosen explicitly in configuration; it is never inferred from what
the host happens to provide. TUFLOW-FV's interface sets it from its own
configuration (a declared switch), not from the presence of `wq%par`.

**oasim's role now:** a diagnostic spectral model. It writes `OAS_`
diagnostics; it does not change `par`, `extc_coef` or feedback. This is
stated in its log at init. `'module'` mode is how oasim becomes the producer
later: the registry, `irr(b)` storage and the manager's `update` hook are
designed so that only a provider interface needs adding. When that happens
oasim would use the modules' banded `light_extinction` (plus its water IOPs)
rather than a separate optics set, so there is one source of Kd.

`light_mode` decides only who produces the in-water field that biology reads.
Who produces the Kd returned to the host is a separate declaration,
`feedback_source` (§5.8(c)), so combinations such as "ELCOM produces light,
oasim produces band Kd" are expressible.

`link_ext_par = .TRUE.` in old nml files maps to `'host'` with a deprecation
line; in that case the host must still declare its PAR Kd (§5.5) or
`validate()` fails with a named message.

### 5.4 In-water storage and orientation

```fortran
TYPE(AED_DPTR),DIMENSION(:),POINTER :: irr      ! irr(b)%p(1:nlev), W/m2, surface = 1
TYPE(AED_DPTR),DIMENSION(:),POINTER :: bioextc  ! bioextc(b)%p(1:nlev), m-1, AED only
```

- **Orientation is normalised in `link_column`**: the manager's views run
  surface = 1 → bottom = `nlev` regardless of host order. For top-down hosts
  this is the host order; for bottom-up hosts the manager keeps a reversed
  index map (or copy-in/out for light arrays only). This fixes D4 at its root
  — §3.3 step 1 shows the existing "local frame" does not.
- Per-band pointers accommodate TUFLOW-FV `par(b,:)` and ELCOM separate
  arrays without copying.
- Module globals are views: `par => irr(par_band)`; `nir/uva/uvb => irr(b)` if
  the band is registered, else derived from PAR (`/par_fraction·frac`) and
  logged as **derived**. The API never writes host-owned irradiance.

### 5.5 `extc_coef` — API-owned PAR Kd (D15)

`extc_coef` is the **total PAR-band Kd** per layer, owned by the API:

- `api` mode: `Kw_PAR + bioKd_PAR`, computed by the manager each substep.
- `host` mode: the host's declared PAR Kd (`bind_band_kd(par_band, kd)`); if
  the host declares none, `validate()` fails. (No "if supplied" inference.)

It is always provided (every mode has a value), and it is never the host's
`extc`/`bioshade` array. The host-dependent meaning today is change **C4**.

### 5.6 Extinction method — one method, array-valued

```fortran
PROCEDURE :: light_extinction => aed_light_extinction      ! name unchanged

SUBROUTINE aed_light_extinction(data, column, layer_idx, extinction)
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction(:)   ! SIZE = aed_n_bands(), default 1
```

- **Accumulate only; never assign.** `aed_common::aed_light_extinction` zeroes
  the array and loops models.
- **Grey modules:** the one-line declaration change; `extinction = extinction
  + Ke*C` broadcasts over bands. The 14 bound implementations: macroalgae,
  macroalgae2 (benthic); chlorophylla, cladophora, macroalgae_clump,
  mixoplankton, knu_cohesive (dev); noncohesive, organic_matter, pathogens,
  pesticides, phyto_abm, phytoplankton, tracer (water). The 3 unbound ones
  (macrophyte, vegetation, habitat) and oasim's commented stub are updated too,
  so re-enabling them later is a one-line change. Out-of-tree modules fail to
  compile until they change the same line.
- **Spectral modules** (Stage 4) query `aed_n_bands()`/`aed_band(b)` in
  `define` and precompute `r(b) = <k(λ)>_b / <k(λ)>_PAR` via `aed_spectra`:
  `extinction(:) = extinction(:) + data%KePHY*C*data%r_phy(:)`.
- Callers: API manager, `glm_aed.F90:1091/1557/1565` (size-1 array).

Declarative alternative (per-variable `light_extinction` attribute, core
computes `Σ Ke·C·r(b)`): recorded as Q8, not now — the method is still needed
for nonlinear cases (canopy caps, cladophora `/dz`).

### 5.7 Spectra — `libaed-water/src/aed_spectra.F90`

The single source of optical data and helpers for modules and for light models:

- `interp`, `calculate_integral_weights` (from `aed_maths`)
- water `a_w`/`b_w`; phyto/sediment/seagrass/CDOM IOPs (from `oasim.inc`)
- `kd_from_iop(a, b_b) = (a + r_s·b_b)/mcosthetas` (oasim :967)
- reference surface spectrum for band weighting (a table only), `band_average`,
  `band_ratio`
- spectral **grid registry** shared by `aed_atmosphere` and `aed_oasim` (§8.2)

Curtin's direct/diffuse split and the atmospheric tables stay in libaed-light
(§8.1). Only the weighting table moves.

### 5.8 Feedback to the host

Two feedbacks, both AED-only, both declared with `bind_feedback`:

**(a) Bio-extinction** `bioextc(b)` per band, every layer (surface included).
Hosts add their own background. Reductions (in the manager, selected by the
host's declaration):

| Host | Accepts | Declared reduction |
|---|---|---|
| ELCOM | per band (after host TODO) | none: `wq%bioextc(col,lev,b)` |
| TUFLOW-FV (WQ API 1.0) | single `bioshade` | `'par'`; per band later via BMT |
| GLM | `ExtcCoefSW`, which drives heating | GLM **registers its heating bands** (`n_bands`, `light_extc`, fractions) and uses `'sw_weighted'`: `Kd_eff = -ln(Σ f_b e^{-Kd_b Δz})/Δz`. With GLM's default single band this reduces to today's PAR value. `glm_api_aed` adds `Kw` and applies GLM's own "not the surface layer" rule |
| SCHISM | none | none |

`bioshade_feedback = .FALSE.` zeroes feedback (a), and only (a); it no longer
touches the internal light field.

**(b) Surface shading** (Q6): `aed_light_shading` → `sw_factor`, today wired
ad hoc in the zone path (`aed_api.F90:1953`, `link_solar_shade`). It becomes a
declared feedback target handled by the manager, with the same
declared-or-absent rule. Semantics: fraction of Q1 that reaches the water; the
host applies it before computing Q2.

**(c) Feedback source** (B7). `feedback_source = 'modules'` (default) returns
the modules' `light_extinction` sum, as in (a). `feedback_source = 'oasim'`
returns band Kd from oasim's spectral calculation instead:

- **Band mapping.** oasim reads the registry in `define` and precomputes
  λ→band overlap weights (`calculate_integral_weights`) from its grid onto the
  host bands (e.g. 33 OASIM bands → ELCOM's UVB/UVA/PAR/NIR).
- **Irradiance-weighted effective Kd**, not a mean of Kd(λ) (the spectrum
  hardens with depth inside a band):
  `Kd_b = -ln(E_b,bot / E_b,top)/Δz`, `E_b = Σ_λ w_λb E(λ)`. This is the same
  routine as GLM's `'sw_weighted'` reduction, `aed_light_utils::effective_kd`.
- **Coverage.** A host band not fully inside oasim's grid (check NIR to
  2800 nm, UVB from 280 nm) falls back to the module sum for that band, logged.
- **Total, declared.** oasim's Kd includes water, so the host declares the
  feedback kind `'total'` and *replaces* its band Kd rather than adding its
  own background. This is an explicit, declared exception to principle 4;
  hosts that can only add (TUFLOW-FV `bioshade`) cannot select
  `'oasim'`, and `validate()` says so. The bio-only alternative
  (`Kd_b(total) − Kd_b(water)`) is Q12.
- **Optics source.** oasim's Kd comes from its IOPs, not from module
  `KePHY`/`KeDOM` calibrations. This is the declared meaning of choosing
  `'oasim'`, logged at init. In `'module'` light mode (Q9) `extc_coef` also
  comes from oasim's PAR band, so biology and host share one optics set.
- **Timing.** oasim computes in `calculate_column` (start of step); feedback is
  returned at end of step, so the host sees Kd from pre-update state — a
  one-substep lag, documented. The manager may instead ask oasim to recompute
  in `return_feedback` through the provider hook shared with `'module'` mode.
- **Transport.** Band Kd reaches the manager through the provider hook (or,
  before it exists, `output=.FALSE.` diagnostics `OAS_Kd_host_<b>` located by
  the manager).

### 5.9 Zones and particles

- Zone aggregation (`fv_aed.F90:1392-1480`; API `zcolm`) handles per-band
  arrays; `z_par/z_nir/…` become views. `z_extc` aggregates the API-owned
  `extc_coef`, not the host array.
- The particle path (`aed_api.F90:1430-1434`) calls the manager's `update`.
- `ATM_` sheet diagnostics are per column, so zones and particles see them via
  the normal diagnostic path; no special case.

### 5.10 Configuration

```
&aed_coupling (or host nml)
light_mode        = 'api'          ! 'api' | 'host'   ('module' reserved)
n_bands           = 1
band_names        = 'PAR'          ! e.g. 'UVB','UVA','PAR','NIR'
band_lambda_lo    = 400.
band_lambda_hi    = 700.
band_sw_fraction  = 0.43           ! api mode
band_Kw           = 0.1            ! api mode
sw_inc_supplied   = .false.        ! normally set by the host interface
bioextc_reduction = 'par'          ! 'par' | 'sw_weighted' | 'none'
feedback_source   = 'modules'      ! 'modules' | 'oasim' (B7; host must accept 'total')
! deprecated, still read: link_ext_par, par_frac, nir_frac, uva_frac, uvb_frac,
!                         base_par_extinction
! continuity switches: see §9
```

Hosts that own their optics fill bands from their own configuration in the
interface; nml values are the fallback.

### 5.11 The `aed_light` manager (`libaed-api/src/aed_light.F90`)

All of §5 is implemented in one API module (pattern of `aed_zones.F90`,
`aed_ptm.F90`). `aed_api.F90` keeps coupling and the column loop and contains
no light arithmetic. The name is free once A2 lands.

| Method | Phase | Role |
|---|---|---|
| `configure(cpl, nml)` | `aed_set_coupling` | read mode, bands, legacy keys, continuity switches |
| `add_band(name, lo, hi, sw_fraction, Kw)` | host interface, init | declare bands |
| `bind_surface(I_0, I_inc, sw_inc_supplied)` | init | Q2 required; Q1 declared present or absent |
| `bind_band(b, irr)` / `bind_band_kd(b, kd)` | init / per column | host mode irradiance and Kd |
| `bind_feedback(kind, b, target, reduce)` | init | `kind` = `'bioextc'` (host adds background), `'total'` (host replaces its Kd; required for `feedback_source='oasim'`) or `'shading'` |
| `validate()` | end of `aed_set_coupling` | freeze registry; check one PAR band, `Σ sw_fraction ≤ 1`, host-mode bindings, feedback targets; **fatal with a named message**; log the full light configuration |
| `provide_globals()` | replaces the light rows of the `BSSOCIATED` table (D16) | `sw_sf` (+`par_sf` alias), `sw_inc`, `sw_inc_is_net`, `par`, `nir/uva/uvb`, `extc_coef`, `rad` if declared |
| `link_column(col, …)` | per column | surface→bottom views (§5.4) |
| `update(column)` | per substep | api: attenuate bands via `aed_light_utils::attenuate_bands`, set `extc_coef`; host: set `extc_coef` from host PAR Kd |
| `return_feedback(column)` | end of step | bio-only `bioextc(:)` every layer, reduce, write; shading |
| `aggregate_zones(…)` | zone mode | per-band zone averaging |

`bind_from_env(env)` reproduces today's behaviour from the existing
`aed_env_t` fields, **including the `bot:top-1` ranges and the host-`extc`
`extc_coef`**, so B2 is a pure extraction. Hosts move to explicit
declarations in B4/B5, and `bind_from_env` is removed at the release.

**Physics split.** Plumbing and policy live in the manager; the Beer-Lambert
kernel is `aed_light_utils::attenuate_bands` (pure, per band, per layer,
surface→bottom) so `glm_aed.F90` shares it. This is the one libaed-light
dependency of libaed-api; the build order (light before api) already allows it.

---

## 6. Light-component roles (summary)

| Component | Library | Owns | Does not |
|---|---|---|---|
| `aed_bands`, `aed_spectra` | water | band registry; optical data; band averaging; spectral grid | compute light |
| `aed_light` manager | api | modes, bindings, validation, globals, in-water field (api mode), `extc_coef`, feedback | spectral physics |
| `aed_light_utils` | light | stateless physics: geometry, atmosphere, surface split, reflectance, `attenuate_bands` | hold state |
| `aed_atmosphere` model | light | incident spectra above the surface (`ATM_`) | in-water light |
| `aed_oasim` model | light | reflection + spectral in-water light as diagnostics; from B7, host-band Kd feedback (`feedback_source='oasim'`) and producer in `'module'` mode | feed biology before B7 |

---

## 7. Host conformance

**GLM (`glm_api_aed.F90`)** — `light_mode='api'`
- `I_0` → SAVEd scalar set each step to net SW below surface/ice (`rad(wlev)`)
  (D2); `I_inc => MetData%ShortWave`, `sw_inc_supplied=.TRUE.`.
- Register GLM heating bands; feedback `'sw_weighted'`; `Kw` addition and the
  surface-layer rule for `ExtcCoefSW` move here.
- Q4: whether legacy `glm_aed.F90` also switches `par_sf` to net.

**TUFLOW-FV (`fv_aed.F90`)** — declared switch between `'api'` and `'host'`
- `irr_host_` naming done (A3).
- Host mode: register bands from `wq%par` first dimension (BMT to confirm
  identity), declare PAR Kd (BMT: is it available?).
- Feedback: `bioshade = reduce(bioextc)`, bio-only (D5); `extc` no longer
  aliased to `extc_coef` (D15).
- `sw_inc_supplied=.FALSE.` unless WQ API exposes gross SW.

**SCHISM (`schism_aed.F90`)** — `light_mode='api'`
- Remove zero `rad` (D3). Verify `srad` is Q2 and whether Q1 is available.
- Remove the interface's own `par/nir/…` calculation (:525-529).

**ELCOM (`elcom_aed.F90` + `elcom_wq.f90`)** — `light_mode='host'`
- Stop forcing `link_ext_par=.false.` (D7); register the four bands with
  wavelength bounds; declare `extc(:,:,1..4)` as band Kd.
- Stop binding the band-sum `rad`. `par_frac` honoured (A6).
- Verify `I_0` (D13); apply albedo in the interface if raw.
- Host side: `bioextc(NC2,NLZ,4)`; complete the `elcom_wq.f90:672-676` TODO.

---

## 8. Module and library changes

### 8.1 Modules

- Extinction signature (B5): §5.6, all implementations incl. unbound ones.
- D14 (A5): macrophyte, vegetation, habitat bodies accumulate. Bit-identical
  (unbound). Re-binding is Q7.
- `par_sf` → `sw_sf` in all consumers (B4).
- `aed_oasim` and dry-cell/canopy users: `sw_inc` + `sw_inc_is_net` (B4, D10).
- D9: `aed_phytoplankton.F90:1586` compares with `par_fraction·Io` or uses a
  depth test (B4).
- Spectral ratios (B6): phytoplankton, mixoplankton, phyto_abm,
  organic_matter (CDOM/POM), noncohesive, knu_cohesive. Canopy modules stay
  grey.

### 8.2 libaed-light

**Target layout** (same concept as `aed_bio_utils`: flat, stateless,
explicit `PUBLIC` list):

```
libaed-water/src/aed_spectra.F90     §5.7
libaed-water/src/aed_bands.F90       §5.2

libaed-light/src/
   aed_light_utils.F90   PURE where possible
      constants   pi, deg2rad, rad2deg, Planck, lightspeed, Avogadro, W->umol   (LR4)
      geometry    zenith_angle
      atmosphere  vH_O3 (one copy, kg m-2; LR1), navy_aerosol_model, slingo,
                  clear_sky_transmittance, sky_spectrum(sky_model, ...)      (LR2)
      surface     direct_diffuse_bands, direct_diffuse_curtin,
                  calculate_lambda_bounds, reflectance
      in-water    attenuate_bands                       (API manager, glm_aed)
                  effective_kd  -ln(E_bot/E_top)/dz over weighted bands
                                                        (GLM 'sw_weighted', oasim B7)
      band sums   band_weights(lo, hi)
   aed_light_tables.F90  solar/atmospheric tables incl. coeff_curtin, compiled once (LR3)
   aed_atmosphere.F90    model: incident surface spectra
   aed_oasim.F90         model: reflection + in-water spectral light
   aed_lighting.F90      factory: 'aed_atmosphere', 'aed_oasim'
```

Removed: `aed_light.F90` (A2), `aed_maths.F90`, `aed_slingo.F90`,
`curtin_light.F90`, `oasim.inc`. `Makefile` and `win/libaed-lighting.vfproj`
follow (and drop the stale `slingo_const.h`, LR6). GLM keeps its C
`zenith_angle` (`glm_util.c:46`); cross-reference comment only.

**`aed_atmosphere` model:**

- Owns everything oasim does before its water-column loop: geometry, ozone,
  water vapour, aerosols, Rayleigh, clouds → direct/diffuse spectra.
- Output: sheet diagnostics `ATM_Ed_dir_<l>`, `ATM_Ed_dif_<l>` (W m⁻² nm⁻¹,
  incident, `rezero=.FALSE.`, **`output=.FALSE.`** via A4), plus `ATM_sw`,
  `ATM_par`, `ATM_uv`, `ATM_O3`, `ATM_zenith` (output).
- Grid: registers the wavelength grid in `aed_spectra`; oasim uses it, so no
  interpolation.
- Hook: `calculate_column` (runs before `calculate_surface` in both drivers —
  `aed_api.F90:2076/2085` vs `:2103`, `glm_aed.F90:1310` vs `:1486` — and
  surface is skipped under ice), plus `calculate_dry` sharing one private
  routine.
- `sky_model=2` scaling input is `sw_inc`, honouring `sw_inc_is_net`.
- Namelist: `&aed_atmosphere` holds `sky_model`, grid, cloud/aerosol/ozone,
  `tzshift`.

**`aed_oasim` linking:**

- `sky_model = 0` = linked: `define` locates the `ATM_` diagnostics; if
  `aed_atmosphere` is absent or listed after oasim the lookup fails with a
  message naming the required order.
- `sky_model = 1, 2` = stand-alone, **reading its existing `&aed_oasim` sky
  keys unchanged** (existing nml files keep working); it calls the same
  `sky_spectrum`, so linked and stand-alone give identical results.
- Reflection applied once, in oasim, skipped when `sw_inc_is_net` (D10).
- Log at init: "diagnostic light model; does not drive biology" (until B7,
  then the declared role: diagnostic, feedback source and/or light producer).
- B7: registry-mapped band Kd for the host (§5.8(c)), using the existing
  `save_Kd` attenuation calculation (`aed_oasim.F90:957-963`), which already
  derives Kd from actual attenuation over the layer.

---

## 9. Change register (behaviour changes)

Everything not listed here is bit-identical at its checkpoint.

| ID | Change | Hosts affected | Checkpoint | Continuity switch |
|---|---|---|---|---|
| **C1** | GLM API `I_0` becomes net SW (D2) | GLM API driver | B4 | `light_compat_surface = .TRUE.` (binds gross, as today) |
| **C2** | oasim reflection skipped on net input (D10) | oasim on TUFLOW-FV, ELCOM, (SCHISM) | B4 | oasim `assume_incident = .TRUE.` |
| **C3** | D9 comparison fixed | all, phytoplankton | B4 | none (bug) |
| **C4** | `extc_coef` = API PAR Kd (D15) | TUFLOW-FV (**large**: switches on daytime light limitation, §3.5), ELCOM (band 1 → total PAR), GLM (SW → PAR; small) | B5 | `light_compat_extc = .TRUE.` (host `extc`, as today) |
| **C5** | Bio-shading feedback reaches top-down and 1-layer columns (D4) | TUFLOW-FV, ELCOM, all 2D | B5 | existing `bioshade_feedback = .FALSE.` |
| **C6** | ELCOM host mode: biology uses ELCOM's band light (D7) | ELCOM | B5 | `light_mode = 'api'` |
| **C7** | Band-specific Kd from spectral modules | multi-band hosts | B6 | grey = don't declare spectra (module `spectral = .FALSE.`) |
| **C8** | GLM heating feedback SW-weighted | GLM with >1 heating band | B6 | `bioextc_reduction = 'par'` |
| **C9** | Macrophyte/vegetation/habitat shading re-enabled (if Q7 = yes) | configs listing them | separate commit | leave unbound / module flag |
| **C10** | A6 `par_frac` honoured | ELCOM with nml ≠ 0.45 | A6 | set `par_frac = 0.45` |
| **C11** | Host band Kd from oasim (total, replaces host Kd) | hosts selecting `feedback_source='oasim'` (ELCOM first) | B7 | `feedback_source = 'modules'` |
| **C12** | `extc_coef` and biology light from oasim (`'module'` mode, Q9) | configs selecting it | B7 | `light_mode` back to `'api'`/`'host'` |
| **C13** | oasim `sky_model=1` with `compute_mean_wind=.false.`: mean wind read from host `mean_wind` (was uninitialised) | configs using it | B1 | none (bug) |
| **C14** | oasim dry cells: `zen` diagnostic uses `tzshift` (was hard-coded 8 h); dead sky block after `RETURN` removed | oasim on wetting/drying hosts, `tzshift` ≠ 8 | B1 | none (bug) |

Continuity switches are removed no earlier than one release after
`light-std-1.0`.

---

## 10. Delivery — the single schedule

Separate repos; default branches `main` (libaed-light, libaed-api, GLM,
ELCOM) and `master` (others). Commits authored as the user, no AI attribution.
Nothing is pushed without asking.

### Phase A — master (small, independent, bit-identical unless in §9)

| # | Repo | Commit | Results |
|---|---|---|---|
| A1 | libaed-api | Stage 0 comment; `LIGHT_PLAN.md` | done (`83c7b81`) |
| A2 | libaed-light | retire `aed_light.F90` (file, factory CASE, Makefile, vfproj) | identical |
| A3 | libaed-fv | D11 `rad_` → `irr_host_` | identical |
| A4 | libaed-water | `output` flag on `aed_variable_t` and both diag definers | identical |
| A5 | libaed-benthic, libaed-riparian | D14: accumulate in the (unbound) macrophyte and vegetation bodies (habitat's body is already fully commented out) | identical |
| A6 | ELCOM | D12 honour `par_frac` | C10 |
| A7 | libaed-light | LR6 drop `slingo_const.h` from vfproj | identical |
| A8 | libaed-api | add `LIGHT_PLAN_REVISED.md` | – |
| A9 | — | **Verify** D4/D15 on Conjimba `_11`/`_12` and ELCOM (feedback ≡ 0; what `extc_coef` holds in TUFLOW-FV). **Capture baselines** for four hosts + oasim (`sky_model` 1, 2; wet and dry) | – |
| A10 | libaed-api | test harness skeleton (§11.2) driving one column through `aed_run_model`; runs against master | – |

### Phase B — branch `standardize_light`

Created from Phase-A master in all affected repos: libaed-water,
libaed-benthic, libaed-dev, libaed-riparian, libaed-demo, libaed-light,
libaed-api, libaed-fv, GLM, schism-aed, ELCOM. `switch_branch.sh <branch>` at
the AED_Tools_Private level checks out a branch in all of them and refuses if
any is dirty. Branch name spelling to confirm (`standardize_light` vs
`standarize_light`).

| CP | Content | Gate | Merge |
|---|---|---|---|
| **B1** light tidy | `aed_spectra` (water); `aed_light_utils`, `aed_light_tables` (light); LR1–LR4; one `sky_spectrum`. Verbatim moves | oasim ≡ A9 baseline | **early** |
| **B2** extract | `aed_bands` (default single band, no readers yet); API `aed_light` manager with `bind_from_env`; `aed_api.F90` delegates; harness targets the manager | all hosts ≡ A9 | **early** |
| **B3** atmosphere | `aed_atmosphere` model; oasim `sky_model=0`; stand-alone oasim nml unchanged | linked ≡ stand-alone; wrong order → named error | **early** |
| **B4** surface | `sw_sf` rename + alias; `sw_inc`/`sw_inc_is_net`; host `I_0`/`I_inc`/`rad` (D2, D3, D13); oasim D10; D9; shading as declared feedback | C1–C3 only | hold |
| **B5** bands & modes | registry filled by hosts; array `light_extinction` in every module repo at once; orientation in `link_column`; `extc_coef` (D15); bio-only feedback (D5); ELCOM host mode; declared globals (D16); `attenuate_bands` in `glm_aed.F90`; remove `bind_from_env` | C4–C6 only; harness: identical `par`/`bioextc` in three orientations | hold |
| **B6** spectral | band ratios in modules; per-band ELCOM feedback; GLM SW-weighted; BMT request | C7–C8 only; ratios vs oasim; CDOM test UV ≫ PAR; GLM energy check | with B4–B5 |

After the release, on a fresh `standardize_light` round (same repo set):

| CP | Content | Gate | Merge |
|---|---|---|---|
| **B7** oasim as provider | provider hook in the manager; `feedback_source='oasim'` with `'total'` feedback kind; `aed_light_utils::effective_kd` (shared with B6 `'sw_weighted'` if not already there); coverage fallback; `'module'` light mode (Q9) | C11–C12 only; oasim band Kd for ELCOM's 4 bands ≡ stand-alone band-integrated calculation; band energy conserved; uncovered band falls back and logs | tag `light-std-1.1` |

B1–B3 are pure refactors or additive; merged early to keep the branch short.
The branch is rebased on master after each early merge. B4–B6 merge together
once all four hosts are regressed and C1–C8 are documented: tag every repo
`light-std-1.0`, release notes from §9.

---

## 11. Verification

1. **A9 first**: D4 and D15 observed in `_11`/ELCOM output before any design
   work depends on them.
2. **Harness** (A10, extended in B2/B5): `libaed-api/test/` builds a synthetic
   column with a tracer that has a known Kd; asserts `par` profile, `extc_coef`
   and `bioextc` for top-down, bottom-up and single-layer columns, and for
   api/host modes. Runs on the Linux VM (`make test`). No test
   infrastructure exists today, so this is a deliverable, not an assumption.
3. **Regression per checkpoint** against A9 baselines: GLM `PoolModel_CGM` +
   one lake; TUFLOW-FV Conjimba short window; SCHISM test; ELCOM test; oasim.
   Differences allowed only for register entries enabled at that checkpoint,
   and each of those must vanish with its continuity switch on.
4. **Init log** (every run): mode, band table, which quantities the host
   supplied/declared, derived bands, `sw_inc_is_net`, feedback targets and
   reductions, continuity switches in force, deprecated keys and globals used.

---

## 12. Open decisions

| # | Question | Recommendation |
|---|---|---|
| Q1 | Band identity: arbitrary host-defined intervals vs fixed four | arbitrary, canonical names recognised |
| Q2 | Who owns `sw_fraction`/`Kw_b` | host, falling back to AED nml |
| Q3 | `I_0`/`sw_sf` = Q2 net | yes (3 of 4 hosts already) |
| Q4 | Legacy `glm_aed.F90` also switches to net | yes, so the two GLM drivers stay identical (C1 applies to both) |
| Q5 | C5 default: feedback on for TUFLOW-FV/ELCOM once fixed | on, with `bioshade_feedback` as the switch |
| Q6 | BMT: identity of `wq%par` bands; is PAR Kd and gross SW available; does `bioshade` affect heating | ask |
| Q7 | Re-enable macrophyte/vegetation/habitat extinction (disabled in csiem 1.5 "temp updates") | decide with the csiem team; separate commit, C9 |
| Q8 | Declarative per-variable extinction | not now; revisit after B6 |
| Q9 | `'module'` light mode (oasim as producer) | implement in B7 with the provider hook |
| Q10 | Continuity switches: individual (as §9) vs one `light_compat` level | individual; clearer in release notes |
| Q11 | Branch name spelling | `standardize_light` |
| Q12 | oasim feedback: total Kd replacing the host's (declared `'total'`) vs bio-only `Kd_b(total) − Kd_b(water)` added to host background | total — the subtraction is inexact for effective Kd and mixes oasim water optics with host `Kw` |

---

## Appendix — key references

| Item | Location |
|---|---|
| API env/data types | `libaed-api/src/aed_api.F90:138-197, 234-290` |
| Coupling type (`Kw` pointer :92) | `aed_api.F90:71-99`, module `Kw` :317, bound :680 |
| Column mapping (orientation kept) | `aed_api.F90:1386-1397` |
| Particle light call | `aed_api.F90:1430-1434` |
| `aed_run_column` light flow | `aed_api.F90:1616-1703` |
| `Light` / `BioExtinction` | `aed_api.F90:2263-2330` |
| Global provision by `BSSOCIATED` (D16) | `aed_api.F90:939-955, 1047-1063` |
| `extc_coef` → host `extc` (D15) | `aed_api.F90:1181` (column), `:1298` (zones) |
| Surface shading feedback | `aed_api.F90:1953` |
| Module extinction interface | `libaed-water/src/aed_core.F90:90, 116, 1246`; `aed_common.F90:534` |
| `aed_light_shading` | `aed_common.F90:597`; `aed_core.F90:1228` |
| Unbound extinction (D14) | `aed_macrophyte.F90:187, 1531-1559`; `aed_vegetation.F90:150, 910-938`; `aed_habitat.F90:98, 1115` |
| `extc_coef` consumers | `aed_phytoplankton.F90:792`, `aed_carbon.F90:407`, `aed_methane.F90:432`, `aed_organic_matter.F90:620`, `aed_habitat_water.F90:241`, `aed_totals.F90:239` |
| Coupling → configure order | `fv_aed.F90:323/325`; `glm_api_aed.F90:225/229`; `schism_aed.F90:458/462`; `elcom_aed.F90:434/437` |
| GLM API bindings | `GLM/src/glm_api_aed.F90:219, 268, 282, 358, 372` |
| GLM legacy driver light | `GLM/src/glm_aed.F90:580, 721, 1091, 1557-1565` |
| TUFLOW-FV | `tuflowfv_wq_api.f90:129-146`; `fv_aed.F90:317, 665, 674, 731-781, 1392-1480` |
| SCHISM | `schism_aed.F90:199, 441, 522-529, 679` |
| ELCOM | `elcom_aed.F90:366, 393, 415, 605, 636-648`; `elcom_wq.f90:279, 410-418, 672-676` |
| oasim | `aed_oasim.F90:125, 223, 438-452, 484-491, 500, 630, 741-850, 885, 967, 1163-1266, 1292` |
| `par > 0.95*Io` (D9) | `aed_phytoplankton.F90:1586` |
