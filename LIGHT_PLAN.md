# AED API light & extinction standard — plan

Status: **DRAFT for review** (2026-09-26). Nothing in this plan is implemented
yet except where marked. Owner of the standard: `libaed-api`.

## 1. Purpose

The AED API (`aed_api.F90`) couples libaed to four hosts — **GLM, TUFLOW-FV,
SCHISM, ELCOM**. Its intent is a *standard*: the API defines what each
environmental quantity means, and each host's interface module
(`glm_api_aed.F90`, `fv_aed.F90`, `schism_aed.F90`, `elcom_aed.F90`) is
responsible for massaging host data to conform. The API must never infer
meaning from which pointers happen to be associated, and host-specific quirks
must live in the host interface, not in the API.

Light is currently the least well-defined part of that contract. This plan:

1. defines each light quantity once, unambiguously (surface irradiance,
   in-water irradiance, extinction);
2. makes the API **multi-band aware**, so hosts that already carry several
   wavebands (ELCOM, TUFLOW-FV, GLM's heating bands) can receive
   **band-specific bio-shading feedback** from AED;
3. lists the defects found in the current implementation and a staged,
   regression-tested path to fix them.

## 2. Trigger: Conjimba run `_12` segfault

`ERA_S4_097_WQ_EX_AED_12` (TUFLOW-FV 2025.2, clean build) died on the first AED
call:

```
START do_aed_models :
forrtl: severe (174): SIGSEGV ... libtuflowfv_exter  aed_apiaed_run_mo
```

Cause: an uncommitted edit (2026-09-24) in `aed_api.F90` `Light` changed the
surface PAR source from `idata%I_0` to `idata%rad(top)`. That was a correct
diagnosis for GLM (see D2) but applied in the shared API; `fv_aed` never
associates `rad` (`fv_aed.F90:754` is commented out), so `rad(top)` is a null
dereference. Run `_11` (same inputs, earlier library) completed normally.

The edit fixed a host-conformance problem by changing the API's semantics for
every host — the pattern this plan is designed to prevent.

## 3. Current-state audit

### 3.1 The distinct light quantities

| # | Quantity | Units | Needed by |
|---|----------|-------|-----------|
| **Q1** | Incident shortwave *above* the water surface (before reflection) | W/m² | `aed_oasim` (applies its own surface reflectance, `aed_oasim.F90:484-491, 885`); terrestrial/dry-cell modules (riparian, vegetation) — to confirm |
| **Q2** | Net shortwave *just below* the surface (after albedo, and ice/snow if the host models them) | W/m² | starting point of the in-water light profile; `par_sf` consumers |
| **Q3** | In-water irradiance per layer, per band (PAR, NIR, UVA, UVB, …) | W/m² | photosynthesis (`par`), photo-degradation / pathogen inactivation (`uva`, `uvb`), `nir` |
| **Q4** | Total net shortwave per layer | W/m² | nobody (API provides global `rad`; no module locates it) |
| **Q5** | Extinction coefficient per layer, per band: background (water) + AED bio/particle contribution | m⁻¹ | in-water light calculation; host heating (via feedback) |

The current names conflate these:

- `I_0` / module global **`par_sf`** — the name says "PAR", but every one of the
  ~20 consumers treats it as **total surface shortwave** (Q1 or Q2 depending on
  host), e.g. `aed_phytoplankton.F90:836` `! surface short wave radiation`,
  `aed_oasim.F90:723` `! Incoming shortwave flux`.
- `rad` — GLM: Q4. TUFLOW-FV: the name `rad_` is used in `fv_aed` for the
  host's **band irradiance array** (Q3). SCHISM: a zero-filled dummy.
  ELCOM: sum of the four band arrays.
- `extc` — sometimes the host's total Kd (GLM `ExtcCoefSW`), sometimes a
  bio-only feedback array (TUFLOW-FV `bioshade`), sometimes band 1 of the
  host's Kd array (ELCOM `extc(:,:,1)`).

### 3.2 Host inventory

| | GLM | TUFLOW-FV | SCHISM | ELCOM |
|---|---|---|---|---|
| Interface module | `GLM/src/glm_api_aed.F90` | `libaed-fv/src/fv_aed.F90` (+ `tuflowfv_external_wq/`) | `schism-aed/src/AED/schism_aed.F90` | `ELCOM/elcom_aed/src/elcom_aed.F90` |
| `I_0` bound to | `MetData%ShortWave` (:268) = **Q1** (gross) | `wq%I_0`, declared "NET SURFACE SHORTWAVE IRRADIANCE" (`tuflowfv_wq_api.f90:129`) = **Q2** | `srad` (:522) — believed Q2, **verify** | `met_data%values(:,SOLAR_RAD_INDEX)` (`elcom_wq.f90:279`), declared "NET" — **raw met input, verify** |
| `rad` bound to | `theLake%Light` (:282) = **Q4** ✅ | not bound (:754 commented) | zeros (:199 "dummied", :679) | Σ nir+par+uva+uvb (:636, :809) |
| Host band irradiance | — (GLM computes `Light` for heating) | `wq%par(nband,NC3)`; `fv_aed` uses band 1 only (:674 `rad_(1,:)`) | none | `wq%par/nir/uva/uvb(col,lev)` (`elcom_wq.f90:410-413`) |
| Host band Kd | `n_bands`, `light_extc` (glm.nml) for heating | per band (host) | none | `wq%extc(col,lev,4)` (`elcom_wq.f90:414-418`); AED reads band 1 only (:605) |
| Bio-shading returned | `ExtcCoefSW` (single; also drives heating) | `wq%bioshade(NC3)` (single) | none | `wq%bioextc(col,lev)` (single) — **added identically to all 4 bands**, `elcom_wq.f90:672-676` "TO DO – MAKE THESE BIOEXTC SPECIFIC TO THE WAVELENGTHS" |
| `link_ext_par` | nml | nml (true → `lpar => rad_(1,:)`) | nml (default F) | **forced F** (:366) although ELCOM computes its own band light |
| `par_fraction` | cpl default 0.43 | nml `par_frac` 0.43 | nml 0.43 | **hard-coded 0.45** (:393), nml ignored |
| Layer order | bottom-up (`bot=1, top=wlev`) | **top-down** (surface cell first) | bottom-up (`top=n_layers`) | top-down (`surf_map`) |

### 3.3 Current API light flow (`aed_run_column`, `aed_api.F90:1616-1703`)

1. `Light(...)` (unless `link_ext_par` or particles): `BioExtinction` sums
   every module's scalar `light_extinction`; `extc = bio + Kw`;
   **writes `idata%extc(bot:top-1)`** (host array); surface
   `par(top) = par_fraction * I_0` (currently `rad(top)` — the `_12` bug);
   attenuate down.
2. `nir/uva/uvb = par / par_fraction * {nir,uva,uvb}_fraction` —
   **overwrites** whatever the host put in those arrays, and assumes all bands
   attenuate like PAR.
3. `calculate_fluxes` …
4. After the substeps: `BioExtinction → idata%bioextc`; if not
   `link_ext_par`, **`idata%extc(bot:top-1) = bioextc + Kw`**; if not
   `bioshade_feedback`, zero `bioextc`.

The `bot:top-1` range and "skip the surface layer" rule were written to match
GLM's `update_light` (comments at :2289-2293, :1693-1696). They are GLM
behaviour embedded in the shared API.

### 3.4 Defects and hazards found

| ID | Severity | Defect |
|----|----------|--------|
| **D1** | **Blocker** | `Light` uses `idata%rad(top)`; null in TUFLOW-FV → segfault (run `_12`). |
| **D2** | High | GLM binds gross shortwave (Q1) to `I_0`, so the API path over-lights by ~7% (albedo never removed). Root cause of the Sep-24 edit. The *host* is non-conformant; the API was right. |
| **D3** | High (latent) | SCHISM binds a zero array to `rad`. Under the Sep-24 edit SCHISM would run silently with **PAR = 0 everywhere**. |
| **D4** | **High — verify** | `extc(bot:top-1)` is an **empty range for top-down hosts** (TUFLOW-FV, ELCOM: `top=1`, `bot=nlev` → `nlev:0`) **and for every single-layer (2D) column** (`1:0`). Bio-shading feedback therefore probably **never reaches TUFLOW-FV or ELCOM**, including 2D Conjimba. Same range in `Light` (:2293). |
| **D5** | High (latent) | If D4 were fixed naively, TUFLOW-FV would get `Kw + bio` written into `bioshade` (`fv_aed.F90:665` `extc => extcoeff_`), double-counting background extinction the host already applies. |
| **D6** | Medium | API overwrites host `nir/uva/uvb` arrays (ELCOM binds them to `wq%nir…`) and derives them by scaling PAR, i.e. assumes band-independent attenuation. |
| **D7** | Medium | ELCOM computes its own band light, but `link_ext_par` is forced off, so the API computes a second, independent PAR field from `I_0` — two inconsistent light fields. |
| **D8** | Medium | `par_sf` is misnamed; it is total shortwave. |
| **D9** | Medium | `aed_phytoplankton.F90:1586` `IF (par > 0.95*Io)` compares PAR with total shortwave — can never be true (≈0.45 ratio). Symptom of D8. |
| **D10** | Medium | `aed_oasim` applies its own surface reflectance, but on TUFLOW-FV/ELCOM (and GLM after a naive D2 fix) it receives *net* shortwave → reflection applied twice. |
| **D11** | Low | `fv_aed` names the host band-irradiance array `rad_`; the commented binding `rad => rad_(:,col)` would bind a wrong slice *and* the wrong quantity (PAR, not total SW). |
| **D12** | Low | ELCOM ignores nml `par_frac` (hard-coded 0.45); others default 0.43. |
| **D13** | Low | ELCOM `I_0` may be gross (raw met `SOLAR_RAD_INDEX`) despite the "NET" comment. |
| **D14** | Medium | The extinction accumulator is wiped by some modules. `aed_common.F90:534` zeroes `extinction` once, then each model adds to it, but `aed_macrophyte.F90:1551/1558` and `aed_vegetation.F90:926` **assign** (`extinction = zero_` / `= MIN(...)`), discarding every earlier module's contribution. The result depends on the order of `&aed_models`. (`aed_macroalgae.F90:2169` `extinction = extinction` is a harmless no-op.) |

## 4. Design principles

1. **The API defines meaning; hosts conform.** Every field has one physical
   definition, units and sign convention. Conversions (albedo, ice, unit
   changes, layer reordering, band mapping) happen in the host interface.
2. **Explicit over inferred.** Behaviour is selected by declared configuration
   (light mode, band registry), never by `ASSOCIATED()` checks on optional
   pointers.
3. **API owns AED's contribution, host owns the medium.** The API computes and
   returns only AED's (bio/particle) extinction. Background water extinction
   belongs to whichever side computes the light.
4. **Orientation-agnostic API.** No `bot:top-1` style ranges; the API works in
   its local surface-to-bottom frame and returns full columns. Host-specific
   rules (e.g. GLM not updating the surface layer's `ExtcCoefSW`) move into
   the host interface.
5. **Calibrations are preserved.** Existing module extinction coefficients
   (`KePHY`, `KeDOM`, …) remain the **PAR-band** specific extinction. Spectral
   information supplies only *band ratios* relative to PAR.
6. **Bit-identical default.** With default configuration (one PAR band,
   `light_mode='api'`), results are identical to the pre-Sep-24 library for
   TUFLOW-FV/SCHISM/ELCOM, and identical to the non-API GLM driver for GLM.

## 5. The standard

### 5.1 Surface quantities

| API env field | Module global | Definition | Required |
|---|---|---|---|
| `I_0` (scalar) | `par_sf` (legacy) **and** `sw_sf` (new alias, preferred) | **Q2**: net shortwave just below the water surface, W/m², after albedo and after ice/snow attenuation | **Yes**, both light modes |
| `I_inc` (scalar, new) | `sw_inc` | **Q1**: incident shortwave above the surface, before reflection, W/m² | Optional. If absent the API sets `I_inc => I_0` and logs a one-time warning |
| `rad` (per layer) | `rad` | **Q4**: net total shortwave per layer, W/m² | Optional, informational. **Never** used to compute light |

`par_sf` keeps its name for backward compatibility but is documented as total
net shortwave; new code uses `sw_sf`. `aed_oasim` (and any module that models
surface reflection or terrestrial exposure) switches to `sw_inc`.

### 5.2 Band registry

Declared once at init by the host interface, via the coupling structure:

```fortran
TYPE aed_light_band_t
   CHARACTER(len=8) :: name      = ''      ! 'PAR','NIR','UVA','UVB' or host-defined
   AED_REAL :: lambda_lo = 400.            ! nm
   AED_REAL :: lambda_hi = 700.            ! nm
   AED_REAL :: sw_fraction = 0.43          ! share of I_0 in this band   (light_mode='api')
   AED_REAL :: Kw = 0.1                    ! background Kd, m-1           (light_mode='api')
END TYPE aed_light_band_t

! in aed_coupling_t
INTEGER :: light_mode = LIGHT_MODE_API      ! LIGHT_MODE_API | LIGHT_MODE_HOST
INTEGER :: n_bands    = 1
TYPE(aed_light_band_t),DIMENSION(:),ALLOCATABLE :: bands
```

**Where the registry lives.** Modules must be able to read the band table in
`define`, and modules depend only on libaed-water. So the band type and the
registry (`n_bands`, `bands(:)`, `par_band`, plus read-only accessors)
live in **libaed-water** (`aed_core`, or a small `aed_bands.F90` next to it),
not in the API or libaed-light. The API `aed_light` manager (§5.10) **fills**
the registry. libaed-water only **holds** it.

**When it is filled.** All four hosts call `aed_set_coupling` before
`aed_configure_models` (TUFLOW-FV `fv_aed.F90:323/325`, GLM
`glm_api_aed.F90:225/229`, SCHISM `schism_aed.F90:458/462`, ELCOM
`elcom_aed.F90:434/437`). The registry is therefore complete before any
module's `define` runs. The legacy `glm_aed.F90` driver (define at :284) never
sets bands and gets the libaed-water default: one PAR band.

No `bio_scale`: a module that does not declare spectral behaviour is grey, so
its Kd is the same in every band. Band-dependent effects come only from modules
that declare them (§5.5).

Rules:

- Bands are arbitrary host-defined wavelength intervals. The canonical names
  `PAR` (400–700), `NIR` (700–2800), `UVA` (315–400), `UVB` (280–315) are
  recognised so the legacy module globals (`par`, `nir`, `uva`, `uvb`) can be
  bound to them.
- Exactly one band must be identified as PAR (by name, or the band covering
  400–700 nm). Its index is `par_band`.
- Default when the host declares nothing: the single PAR band with
  `sw_fraction = par_fraction`, `Kw = base_par_extinction` — i.e. today's
  behaviour.
- In `api` mode `Σ sw_fraction ≤ 1` is checked; in `host` mode `sw_fraction`
  and `Kw` are ignored (the host owns them).
- Legacy nml `nir_frac/par_frac/uva_frac/uvb_frac` map onto the canonical bands
  when present.

### 5.3 Light modes

| | `light_mode = 'api'` | `light_mode = 'host'` |
|---|---|---|
| Who computes in-water light | API | host |
| Host supplies | `I_0` (Q2) | `I_0` (Q2) **and** irradiance per band per layer |
| API computes | `irr_b(top) = sw_fraction_b · I_0`; attenuate each band with `Kw_b + bioKd_b` | nothing (reads host irradiance) |
| API returns | `bioextc_b` per band, bio-only | `bioextc_b` per band, bio-only |
| Intended hosts | GLM, SCHISM | ELCOM, TUFLOW-FV |

`link_ext_par = .TRUE.` in existing nml files is accepted as a synonym for
`light_mode = 'host'` (with a deprecation note in the log).

### 5.4 In-water irradiance storage

Per column the API holds one 1-D layer array per band:

```fortran
TYPE(AED_DPTR),DIMENSION(:),POINTER :: irr      ! irr(b)%p(1:nlev), W/m2
TYPE(AED_DPTR),DIMENSION(:),POINTER :: bioextc  ! bioextc(b)%p(1:nlev), m-1, AED contribution only
```

An array of per-band pointers accommodates both host layouts without copying:
TUFLOW-FV `par(b, top:bot)` (strided per band) and ELCOM separate `par/nir/…`
arrays. In `api` mode the API allocates the storage itself.

Module-facing names are **views** onto bands: `par => irr(par_band)`,
`nir/uva/uvb => irr(<band>)` if that band exists. If it does not exist
(e.g. single-band config), the API keeps today's derived behaviour
(`irr_PAR / par_fraction * frac`) and flags it in the log as derived. The API
never writes into host-owned irradiance in `host` mode.

`extc_coef` (module global) is the API's internal **PAR-band total** Kd
(`Kw_PAR + bioKd_PAR` in `api` mode; host PAR Kd if supplied in `host` mode),
used by `photosynthesis_irradiance` etc.

### 5.5 Extinction: module interface (one method, array-valued)

The existing method is kept. **No new method is added.** Its return argument
becomes an array with one element per registered band:

```fortran
! aed_core.F90 -- type-bound, unchanged name
PROCEDURE :: light_extinction => aed_light_extinction

SUBROUTINE aed_light_extinction(data, column, layer_idx, extinction)
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction(:)   ! was scalar; SIZE = n_bands (default 1)
```

- **Semantics:** add this module's Kd contribution (m-1) for each band to
  `extinction(b)`. **Accumulate only; never assign** (fixes D14).
- **Grey modules need no code change beyond the declaration.** Existing
  bodies of the form `extinction = extinction + (Ke*conc)` broadcast the
  scalar over every band under Fortran array semantics. Each of the 17
  implementations changes one line (`AED_REAL,INTENT(inout) :: extinction`
  → `extinction(:)`). The compiler catches any implementation that was missed.
- **Default = today.** With one band, `extinction(1)` is the PAR Kd, which is
  bit-identical to the scalar path.
- **Spectral modules** read the registry in `define` and precompute band
  ratios `r(b) = <k(λ)>_b / <k(λ)>_PAR` for each constituent. The calibrated
  PAR coefficient is unchanged, and runtime cost is one multiply per band:
  `extinction(:) = extinction(:) + data%KePHY*conc*data%r_phy(:)`.
- **A module does not need to know the bands in advance.** It only needs to
  find them before it precomputes. The registry is complete before `define`
  (§5.2), so a module that wants spectral behaviour asks for `n_bands` and
  the wavelength bounds of each band, then builds `r(:)` over whatever bands
  the host declared. A module that does not ask stays grey.
- `aed_common::aed_light_extinction(column, layer_idx, extinction(:))`
  zeroes the array and loops the models, as today. Callers (`aed_api`
  `Light`/`BioExtinction`, `glm_aed.F90:1091/1557/1565`) pass an array of
  size `n_bands`. The legacy GLM driver passes a size-1 array.

Spectral shapes for the modules that opt in (Stage 4):

| Module | Constituent | Spectral shape (source) |
|---|---|---|
| `aed_phytoplankton` (+ mixoplankton, phyto_abm) | Chl-a | per-group chlorophyll-specific absorption spectrum (oasim IOP set) |
| `aed_organic_matter` | CDOM/DOC | `a(λ) = a_ref · exp(-S(λ-λ_ref))`, strong in UV |
| `aed_organic_matter` | POM/detritus | exponential with a shallower slope (oasim POM IOP) |
| `aed_noncohesive`, `knu_cohesive` | mineral particles | roughly flat / weak λ-dependence |
| macrophyte, macroalgae, vegetation, cladophora | canopy | grey (no change) |

- **Where the spectra live.** Band ratios are computed inside modules
  (libaed-water, benthic, …), so the spectra and the band-averaging helper
  must be in **libaed-water** (`aed_spectra.F90`), not libaed-light. This is
  because of build order. `aed_oasim` then `USE`s it from libaed-water, so
  there is one source of optics. The reference surface spectrum for band
  weighting is a small table (from `curtin_light`) that moves along with it.
- Kd is approximated from absorption (plus backscatter where available). This
  ratio-based approximation is adequate for band feedback. `aed_oasim` remains
  the high-fidelity spectral option.

**Alternative considered (declarative).** `aed_variable_t%light_extinction`
already exists (`aed_core.F90:116`) and is never set. Modules could declare a
specific Kd (plus an optional spectral-shape id) when they define a variable,
and the core would compute `Σ Ke·C·r(b)` itself. That is purer for the linear
majority. However, the method must remain for the nonlinear cases (macrophyte
canopy cap, cladophora `/dz`, vegetation), which would leave two paths. Record
this as a possible later simplification (§10 Q7). The array-valued method is the
Stage 3 choice.

### 5.6 Feedback to the host

The API returns `bioextc(b)` — **AED contribution only, per band, every
layer** (surface included). Hosts add their own background. Host interfaces
reduce to what the host accepts:

| Host | Accepts | Reduction (in host interface) |
|---|---|---|
| ELCOM | per band (after host TODO is completed) | direct: `wq%bioextc(col,lev,b)`; `elcom_wq.f90:672-676` adds band-specific values |
| TUFLOW-FV (WQ API 1.0) | single `bioshade` | PAR band (default); later per band via WQ API 2.0 request to BMT |
| GLM | single `ExtcCoefSW` (drives heating too) | `ExtcCoefSW = Kw + bioextc(PAR)` excluding the surface layer — GLM's rule, now in `glm_api_aed`. Option: SW-weighted effective Kd `-ln(Σ f_b e^{-Kd_b Δz})/Δz` to match GLM `light_extc` bands |
| SCHISM | none | not returned (option for future heat coupling) |

A helper `aed_bioextc_reduce(mode, ...)` in the API provides `'par'` and
`'sw_weighted'` reductions so hosts don't reimplement them.

`bioshade_feedback = .FALSE.` continues to mean "return zeros to the host";
it no longer affects the API's internal light field.

### 5.7 Layer orientation

The API already maps each column to a local frame in `aed_run_model`
(:1386-1397). All light code operates on the local frame with explicit
surface→bottom iteration; no index-range arithmetic that assumes bottom-up.
Regression tests cover top-down, bottom-up and single-layer columns (D4).

### 5.8 Zones and particles

- Zone aggregation (`fv_aed.F90:1392-1480`, API `zcolm` for `par_sf`, `nir`,
  `uva`, `uvb`, `bioextc`) is extended to per-band arrays; the legacy
  `z_par/z_nir/…` become views as in §5.4.
- The particle path (`aed_api.F90:1430-1434`) calls the same banded light
  routine.

### 5.9 Configuration summary (`&aed_coupling` or host nml)

```
light_mode        = 'api' | 'host'        ! default 'api'
n_bands           = 1                     ! default
band_names        = 'PAR'                 ! e.g. 'UVB','UVA','PAR','NIR'
band_lambda_lo    = 400.
band_lambda_hi    = 700.
band_sw_fraction  = 0.43                  ! api mode only
band_Kw           = 0.1                   ! api mode only
bioextc_reduction = 'par' | 'sw_weighted' ! single-band hosts
! legacy, still accepted: link_ext_par, par_frac, nir_frac, uva_frac, uvb_frac, base_par_extinction
```

Hosts that own their optics (ELCOM, TUFLOW-FV) fill the registry from their
own configuration in the interface module; the nml values are a fallback.

### 5.10 The `aed_light` API module (manager)

Everything in §5.1-5.9 is implemented in **one new API module,
`libaed-api/src/aed_light.F90` (`MODULE aed_light`)**. It follows the pattern of
`aed_zones.F90` and `aed_ptm.F90`. `aed_api.F90` keeps the coupling and the
column loop, and delegates all light work to it. `aed_api.F90` then no longer
contains any light arithmetic, band logic or orientation assumptions.

The module name is freed by retiring `libaed-light/src/aed_light.F90` (§7.1).

**Object.** A single `aed_light_t` instance (module-private, SAVEd) holds the
mode, the band table (which it pushes into the libaed-water registry), surface
bindings, per-column band pointers, feedback targets and diagnostics state.

| Method | Phase | Role |
|---|---|---|
| `configure(cpl, nml)` | init, from `aed_set_coupling` | read `light_mode` and legacy `link_ext_par` / `*_frac` / `base_par_extinction`; defaults |
| `add_band(name, lo, hi, sw_fraction, Kw)` | init, host interface | declare bands; canonical names map to `par/nir/uva/uvb` |
| `bind_surface(I_0, I_inc)` | init | Q2 is required and Q1 optional (falls back to `I_0` with a one-time warning) |
| `bind_band(b, irr)` | init/column (host mode) | point band `b` at host irradiance, whatever its layout (TUFLOW-FV strided `par(b,:)`, ELCOM separate arrays) |
| `bind_feedback(b, target, reduce)` | init | where `bioextc` goes: per band, `'par'` or `'sw_weighted'` reduction, or none |
| `validate()` | end of `aed_set_coupling`, before `define` | one PAR band exists; `Σ sw_fraction ≤ 1` (api mode); host mode has band irradiance bound; **fatal with a named message** rather than a segfault (D1-class errors become init errors); logs the band table, mode and derived-vs-host status |
| `provide_globals()` | after `define` | point module globals `par_sf`/`sw_sf`/`sw_inc`/`par`/`nir`/`uva`/`uvb`/`extc_coef` at the right storage |
| `link_column(col, top, bot)` | per column | per-band views in the API's local surface→bottom frame (§5.7) |
| `update(column)` | per substep | api mode: attenuate each band with `Kw_b + bioKd_b`; host mode: no-op (reads the host). Derives missing canonical bands only when declared as derived |
| `return_feedback(column)` | end of step | compute bio-only `bioextc(:)` for every layer through `aed_light_extinction(…, ext(:))`, apply the declared reduction, write to host targets; zeros if `bioshade_feedback=.FALSE.` |
| `aggregate_zones(…)` | zone mode | per-band zone averaging (replaces `z_par/z_nir/…` special cases) |

**Fallback rules (all explicit, all logged once at init):**

1. `I_inc` not bound → `I_inc => I_0`, with a warning (`aed_oasim` then applies
   reflection to net light; it is flagged).
2. Canonical band not declared but a module uses `nir/uva/uvb` → derived from
   PAR via `frac/par_fraction`, logged "derived".
3. `link_ext_par=.TRUE.` → `light_mode='host'`, with a deprecation note.
4. Host declares no bands → the single PAR band from the nml fractions (today's
   behaviour).
5. Nothing is ever inferred from `ASSOCIATED()` at run time. Every choice is
   made in `validate()`.

**Compatibility path.** `bind_from_env(env)` builds the bindings from the
current `aed_env_t` fields (`I_0`, `par`, `extc`, `bioextc`, `nir/uva/uvb`).
Existing host interfaces therefore keep working unchanged in Stage 1. Hosts
move to the explicit `add_band`/`bind_*` calls in Stages 2-3.

**Physics split.** The manager owns *plumbing and policy*. The Beer-Lambert
kernel (per band, per layer, surface→bottom) lives in libaed-light as a pure
procedure (`attenuate_bands` in `aed_light_utils`, §7.1) so that `glm_aed.F90`
can call the same code. The manager calls it. Spectra live in libaed-water (§5.5), and full
spectral models (`aed_oasim`) remain libaed-light *modules*.

**Is it a good idea? Yes, because:**

- It gives a single place for the standard, so each host conforms to one
  interface rather than to scattered `aed_api` internals.
- Validation at init converts the whole D1/D3 class of silent or segfault
  failures into clear messages.
- It removes GLM-specific rules from shared code (§3.3).
- It is testable in isolation: the unit harness in §9 drives `aed_light`
  directly.

The cost is one more module and one extraction refactor (Stage 1), which is
kept bit-identical.

## 6. Host conformance changes

**GLM (`glm_api_aed.F90`)** — `light_mode='api'`
- `I_0` → a SAVEd scalar set each step in `api_do_glm` to `rad(wlev)`
  (net, includes ice/snow). Fixes D2 at the source.
- `I_inc => MetData%ShortWave`.
- `ExtcCoefSW` update (excluding surface layer) moves here from the API.
- Decision needed (§10 Q4): also make `glm_aed.F90` bind `par_sf` to net SW so
  the two GLM drivers stay identical.

**TUFLOW-FV (`fv_aed.F90`)** — `light_mode='host'` when TUFLOW-FV provides band light, else `'api'`
- Rename `rad_` → `irr_host_` (D11); delete the commented `rad` binding.
- Register bands from `wq%par` first dimension (identity/wavelengths from
  TUFLOW-FV config — **confirm with BMT**).
- `extc => extcoeff_` replaced: write `bioshade = reduce(bioextc)`, bio-only (D5).
- `I_0` unchanged (already Q2).

**SCHISM (`schism_aed.F90`)** — `light_mode='api'`
- Remove the zero `rad` binding (D3).
- Verify `srad` is net (Q2); if gross, apply SCHISM's albedo in the interface.
- Remove the interface's own `par/nir/…` calculation (:525-529) — the API owns it.

**ELCOM (`elcom_aed.F90` + `elcom_src_v3/src/elcom_wq.f90`)** — `light_mode='host'`
- Stop forcing `link_ext_par=.false.` (D7); register the four ELCOM bands with
  their wavelength bounds.
- Stop binding the band-sum `rad`.
- Honour nml `par_frac` (D12).
- Verify `I_0` is net; if raw met, apply ELCOM albedo in the interface (D13).
- Host side: `bioextc(NC2,NLZ)` → `bioextc(NC2,NLZ,4)`; complete the
  `elcom_wq.f90:672-676` TODO with band-specific addition.

## 7. Module changes

- **Extinction signature (Stage 3):** `extinction` → `extinction(:)` in
  `aed_core` (default), `aed_common` (accumulator) and all 17 implementations:
  macroalgae, macroalgae2, macroalgae_clump, macrophyte, chlorophylla,
  cladophora, habitat, mixoplankton, knu_cohesive, vegetation, noncohesive,
  organic_matter, pathogens, pesticides, phytoplankton, phyto_abm, tracer.
  Out-of-tree modules fail to compile until they change the same one line.
- **D14:** macrophyte and vegetation change from assigning to accumulating
  (`extinction(:) = extinction(:) + MIN(kemac, max_extc)`). This is a
  behaviour change wherever those modules are listed after light-absorbing
  modules, so it is noted in release notes.
- All `par_sf` consumers: documentation only (meaning = net SW); new code uses `sw_sf`.
- `aed_oasim`: `sw_inc` instead of `par_sf` (D10).
- `aed_phytoplankton.F90:1586`: compare against `par_fraction*Io` or use a
  depth test (D9).
- Spectral band ratios inside `light_extinction` (Stage 4): phytoplankton,
  organic_matter, noncohesive/cohesive, mixoplankton, phyto_abm.
- Review riparian/vegetation/macroalgae dry-cell use of `par_sf` vs `sw_inc`.

### 7.1 libaed-light restructure

**Current state (review):**

| File | Lines | Content | Fate |
|---|---|---|---|
| `aed_light.F90` | 297 | Unused stub model (diagnostics with hard-coded 0.43/0.05; `BioExtinction`/`update_light` under `#if 0`) | **Retire.** Frees the name for the API module |
| `aed_atmosphere.F90` | 338 | Stub, a "copy of aed_dummy" writing sine-wave test variables, plus `vH_O3` | **Rewrite** as the real atmosphere model (§7.2) |
| `aed_maths.F90` | 251 | `interp` (generic), `calculate_integral_weights`, `pi`/`deg2rad` | → `aed_spectra` (libaed-water) |
| `aed_slingo.F90` | 393 | Slingo cloud, `vH_O3`, `zenith_angle`, Navy aerosol; all `PUBLIC` by default | → `aed_light_utils` |
| `curtin_light.F90` | 414 | Direct/diffuse spectral split, `calculate_lambda_bounds` | → `aed_light_utils` |
| `aed_oasim.F90` | 1375 | The in-water spectral model | Keep; atmosphere part moves out (§7.2) |
| `oasim.inc` | 2326 | Spectral tables | Split: water/IOP → `aed_spectra`; solar/atmospheric → `aed_light_tables` |
| `aed_lighting.F90` | 114 | Factory | Keep: `'aed_atmosphere'`, `'aed_oasim'` |

**Defects found in the review:**

- **L1:** `vH_O3` is defined twice with different units. The `aed_slingo` copy
  (used by oasim) returns kg m⁻². The `aed_atmosphere` copy returns matm-cm,
  because it lacks the `*0.4462e-3*48/1000` conversion.
- **L2:** The ~100-line atmospheric transmittance block is copied between
  `aed_calculate_column_oasim` (:741-850) and `aed_calculate_dry_oasim`
  (:1163-1266). The copies differ only in the diagnostics they write.
- **L3:** `oasim.inc` is `#include`d into both `aed_oasim` and `curtin_light`,
  so every table, including the non-`PARAMETER` `coeff_curtin(5,2,140,9)`,
  exists twice.
- **L4:** `deg2rad` is defined in three files, and `Planck`/`lightspeed`/`Avogadro`
  locally in oasim.
- **L5:** `aed_calculate_surface_oasim` is empty ("ASPECTS OF ATMOSPHERE MODEL
  … CAN BE CONSIDERED TO COME HERE"), which is the intent realised in §7.2.

**Target layout** (same concept as `aed_bio_utils`: a flat, stateless
collection of shared calculations with an explicit `PUBLIC` list):

```
libaed-water/src/aed_spectra.F90     (optics that biology modules need; §5.5)
   interp, calculate_integral_weights            <- aed_maths
   water a_w/b_w; phyto/sediment/seagrass IOPs   <- oasim.inc
   kd_from_iop(a, b_b) = (a + r_s*b_b)/mcosthetas  <- oasim :967
   band_average, band_ratio                       -> r(b) of §5.5
   spectral grid registry (lambda, bounds)        -> shared ATM/OAS grid, §7.2

libaed-light/src/
   aed_light_utils.F90   PURE where possible
      constants   pi, deg2rad, rad2deg, Planck, lightspeed, Avogadro, W->umol
      geometry    zenith_angle
      atmosphere  vH_O3 (one copy, kg m-2), navy_aerosol_model, slingo,
                  clear_sky_transmittance, sky_spectrum(sky_model, ...)  <- L2
      surface     direct_diffuse_bands, direct_diffuse_curtin,
                  calculate_lambda_bounds, reflectance
      in-water    attenuate_bands   <- the Beer-Lambert kernel (API aed_light
                                       and glm_aed.F90)
      band sums   band_weights(lo, hi)  for PAR / SW / UV
   aed_light_tables.F90  solar/atmospheric tables, compiled once (L3)
   aed_atmosphere.F90    model: incoming surface spectra (§7.2)
   aed_oasim.F90         model: in-water spectral light
   aed_lighting.F90      factory
```

Removed: `aed_light.F90`, `aed_maths.F90`, `aed_slingo.F90`,
`curtin_light.F90`, `oasim.inc`. Build lists (`Makefile`,
`win/libaed-lighting.vfproj`) change to match. There are no references outside
libaed-light apart from `aed_common`'s `USE aed_lighting`, which is unchanged.
Retiring `aed_light` must land **before or with** the API module to avoid a
duplicate `aed_light.mod`.

GLM keeps its C `zenith_angle` (`glm_util.c:46`, different signature). Add a
cross-reference comment only.

### 7.2 `aed_atmosphere` as the sky model

`aed_atmosphere` becomes a real model that owns the "above the water"
part of oasim: everything before the water-column loop. That covers solar
geometry, ozone, water vapour, aerosols, Rayleigh, clouds and the resulting
direct/diffuse spectra at the surface. `aed_oasim` is left with surface
reflectance and in-water propagation.

- **nml:** `sky_model`, `lambda_method`/`nlambda`/`lambda_min`/`lambda_max`,
  cloud/aerosol/ozone options and `tzshift` move from `&aed_oasim` to
  `&aed_atmosphere`.
- **Output:** sheet diagnostics `ATM_Ed_dir_<l>` and `ATM_Ed_dif_<l>` (W m⁻² nm⁻¹,
  above-surface, `rezero=.FALSE.`), plus scalar `ATM_sw`, `ATM_par`, `ATM_uv`,
  `ATM_O3` and `ATM_zenith`.
- **Grid:** atmosphere registers the wavelength grid in the `aed_spectra` grid
  registry. oasim uses the same grid, so no interpolation happens at the interface.
- **Hook: `calculate_column`, not `calculate_surface`.** Both drivers call
  `aed_calculate_column` *before* `aed_calculate_surface` (API
  `aed_api.F90:2076/2085` vs `:2103`; GLM `glm_aed.F90:1310` vs `:1486`), and
  the API skips surface calls under ice. `calculate_column` runs every step in
  `&aed_models` order. `aed_atmosphere` computes once, at the column surface,
  and also implements `calculate_dry` (shared private routine) for dry cells.
- **Order is self-enforcing:** when linked, oasim's `define` locates the `ATM_`
  diagnostics. If `aed_atmosphere` is not listed before `aed_oasim`, that
  lookup fails with a message naming the required order.
- **Linking in oasim (explicit, backward compatible):** `sky_model = 0` means
  "linked from `aed_atmosphere`". The existing values (1 = OASIM clear/cloudy
  sky, 2 = Curtin split of measured SW) keep working stand-alone, but oasim
  then calls the same `aed_light_utils::sky_spectrum` rather than its own
  copy (L2), so the results are identical either way.
- **Other consumers** can link to `ATM_` outputs in the same way, e.g. riparian
  or terrestrial modules wanting UV or PAR above the canopy, or `sw_inc` (Q1) for
  hosts that don't supply it.
- **Surface convention:** `ATM_` spectra are *incident* (Q1-like).
  Reflection is applied once, in oasim (resolves D10 for the linked path).
  With `sky_model = 2` (spectra derived from measured SW) the scaling
  input must be `sw_inc`, not the net `par_sf`.
- **Open item:** `aed_variable_t` has `rezero` but no "don't write to
  output" flag, so 2×`nlambda` diagnostics would all be written by
  hosts. Proposed: an optional `output=.FALSE.` attribute on
  `aed_define_sheet_diag_variable`, a small `aed_core` addition that hosts
  honour when building their output lists (§10 Q9).

## 8. Staged delivery

| Stage | Scope | Acceptance |
|---|---|---|
| **0 — unblock** (done) | Revert `aed_api.F90` `Light` to `I_0` (drop the Sep-24 edit). | Conjimba `_12` runs; results match `_11` over the overlap. |
| **1 — extract** | Retire libaed-light `aed_light` (§7.1). Create API `aed_light` (§5.10) and move today's `Light`, `BioExtinction`, `nir/uva/uvb` derivation and `extc` write into it **unchanged**, including the current `bot:top-1` ranges. `bind_from_env` only; no host changes. | Bit-identical on all four hosts (GLM `PoolModel_CGM`, Conjimba short window, SCHISM and ELCOM tests). |
| **2 — surface semantics** | §5.1; `validate()` with fallbacks; host bindings for `I_0`/`I_inc`/`rad` (GLM D2, SCHISM D3, ELCOM D13 check); oasim → `sw_inc` (D10); D9; D12. | GLM API driver bit-identical to `glm_aed` driver on `PoolModel_CGM` (the Sep-24 check). TUFLOW-FV/SCHISM bit-identical to Stage 1. |
| **3 — bands & modes** | Band registry in libaed-water; array-valued `light_extinction` (§5.5, grey only) + D14; `light_mode`; orientation fix (D4) and bio-only feedback (D5); ELCOM host mode (D6, D7); TUFLOW-FV rename (D11); `aed_light_utils::attenuate_bands` shared with `glm_aed.F90`; explicit `add_band`/`bind_*` in hosts. | Default single-band config bit-identical to Stage 2 (except D14 where it applies). New tests: top-down / bottom-up / 1-layer column feedback non-empty; `bioshade` nonzero in a TUFLOW-FV test. ELCOM 4-band grey reproduces current ELCOM Kd. |
| **4 — spectral physics** | `aed_spectra` in libaed-water; band ratios in modules; ELCOM per-band host feedback; SW-weighted reduction for GLM; BMT request for per-band `bioshade`. | Band ratios checked against oasim for a reference water; CDOM-rich test shows UV Kd ≫ PAR Kd; energy sanity on GLM heating. |
| **L1 — light repo tidy** (parallel; with Stage 1) | §7.1 layout: `aed_light_utils`, `aed_light_tables`, `aed_maths`/interp → `aed_spectra`; fix L1-L4; oasim's column and dry paths call one `sky_spectrum`. Code is moved **verbatim**. | oasim test case bit-identical (`sky_model` 1 and 2, wet and dry cells). |
| **L2 — atmosphere split** | §7.2: `aed_atmosphere` model, `ATM_` outputs, spectral grid registry, oasim `sky_model = 0`; optional non-output diagnostics (Q9). | oasim with `sky_model = 0` + `aed_atmosphere` bit-identical to stand-alone oasim with the same sky settings; wrong model order gives the named error. |

Note on D4: once feedback actually reaches TUFLOW-FV/ELCOM, results **will
change** for existing sims (bio-shading switches on). Treat as a deliberate,
documented behaviour change, with a flag to reproduce the old (no-feedback)
behaviour for continuity runs.

### 8.1 Git staging

**Principle:** small, independent, non-destructive items go to master first,
each as its own commit. The redesign then happens on one branch,
**`standardize_light`**, and master only receives it at verified checkpoints.

**State at 2026-09-27:** Stage 0 put `aed_api.F90` back to the committed code
(the Sep-24 edit was never committed). Master is therefore clean. Unrelated
local changes must stay out of light commits:
`libaed-benthic/src/aed_bivalve.F90` and `schism-aed/build_local.sh`. Default
branches differ: `main` for libaed-light, libaed-api, GLM and ELCOM; `master`
for the rest.

**Phase A — master (prep; each commit bit-identical unless marked):**

| # | Repo | Commit | Results |
|---|---|---|---|
| A1 | libaed-api | Stage 0 comment in `Light`; add `LIGHT_PLAN.md` | identical |
| A2 | libaed-light | Retire unused `aed_light.F90` (§7.1: file, factory CASE, Makefile, vfproj) | identical; no nml uses it |
| A3 | libaed-fv | D11: rename `rad_` → `irr_host_`; delete the commented `rad` binding | identical |
| A4 | libaed-water | Q9: optional `output=` on `aed_define_sheet_diag_variable` (default `.TRUE.`) and the matching `aed_variable_t` field. Hosts honour it later | identical |
| A5 | libaed-benthic, libaed-riparian | D14: macrophyte/vegetation accumulate rather than assign. **Changes results** where these are listed after absorbing modules. Separate commits, with that stated in the message (if Q8 = unconditional) | **changes** |
| A6 | ELCOM | D12: honour nml `par_frac` | changes only if nml ≠ 0.45 |
| A7 | — (no commit) | Verify D4 on Conjimba `_11`/`_12` and an ELCOM case. **Capture regression baselines** for all four hosts plus an oasim case (`sky_model` 1 and 2, wet and dry). Every later "bit-identical" check compares against these | — |

Anything that depends on a design decision (§10 Q3/Q4: D2, D9, D10, `sw_sf`)
stays off master.

**Phase B — branch `standardize_light`**, created from the A-level master in
every affected repo: libaed-water, libaed-benthic, libaed-dev, libaed-riparian,
libaed-light, libaed-api, libaed-fv, GLM, schism-aed, ELCOM. All repos must be
on the same branch to build. A small `switch_branch.sh <branch>` at the
AED_Tools_Private level checks out the branch in all of them and refuses to run
if any repo has uncommitted changes.

Commit order (one logical step per commit; each builds on all hosts):

| Checkpoint | Commits | Gate | Merge to master? |
|---|---|---|---|
| **B1** light tidy (L1) | `aed_spectra` (interp/weights/water+IOP tables) in water; `aed_light_utils` + `aed_light_tables` in light; oasim uses one `sky_spectrum` | oasim bit-identical to A7 | **Yes**: pure refactor |
| **B2** extract (Stage 1) | API `aed_light` module with `bind_from_env`; `aed_api.F90` delegates to it | all hosts bit-identical to A7 | **Yes**: pure refactor |
| **B3** atmosphere split (L2) | `aed_atmosphere` model; oasim `sky_model = 0` | linked = stand-alone oasim | Yes (additive) |
| **B4** surface semantics (Stage 2) | per host: GLM D2, SCHISM D3, ELCOM D13; `sw_sf`/`sw_inc`; oasim D10; D9 | Stage 2 acceptance | Hold |
| **B5** bands & modes (Stage 3) | registry; array `light_extinction` (all module repos in one coordinated set); D4/D5 orientation and feedback; ELCOM host mode; `attenuate_bands` | Stage 3 acceptance | Hold |
| **B6** spectral (Stage 4) | band ratios in modules; per-band feedback | Stage 4 acceptance | Merge B4-B6 together once all four hosts are regressed and the behaviour changes are documented |

B1-B3 don't change results, so merging them early keeps the branch short-lived.
B4-B6 change results and wait on §10 decisions, so they stay on the branch
until they pass as a set. The branch is rebased on master after each merge.
The merge of B4-B6 is the release point: tag every repo (e.g.
`light-std-1.0`), and write release notes listing D2, D4, D9, D10 and D14
effects and the continuity flags.

## 9. Verification plan

1. **D4 confirmation first** (before Stage 3): inspect `bioshade` / extinction
   outputs of Conjimba `_11` and an ELCOM test — expect bio contribution ≡ 0.
2. Unit harness in `libaed-api` driving one synthetic column through
   `aed_light` (and `aed_run_model`) in three orientations; assert `par` profile and
   `bioextc` identical across orientations.
3. Regression per host at each stage (GLM: `PoolModel_CGM` + a lake case;
   TUFLOW-FV: Conjimba short window; SCHISM: test case; ELCOM: test case).
4. Log at init: light mode, band table, which fields the host supplies,
   derived-vs-host status of `nir/uva/uvb`, `I_inc` fallback warning.

## 10. Open decisions

1. **Band identity** — arbitrary host-defined wavelength intervals with
   canonical names as defaults (recommended), or fixed four bands?
2. **Who owns `sw_fraction` / `Kw_b`** — host (recommended, falls back to
   AED nml) or AED nml only?
3. **`par_sf` meaning** — net below-surface SW, Q2 (recommended; 3 of 4 hosts
   already supply it) vs incident Q1.
4. **GLM `glm_aed.F90`** — change `par_sf` to net as well (drivers stay
   identical; GLM results shift slightly) or leave and document?
5. **D4 behaviour change** — enable bio-shading feedback for TUFLOW-FV/ELCOM by
   default once fixed, or opt-in?
6. **TUFLOW-FV band identity** — what is the first dimension of `wq%par`, and
   does TUFLOW-FV apply `bioshade` to PAR only or to its heating? (ask BMT)
7. **Declarative extinction**: later, move linear modules to per-variable
   `light_extinction` attributes (core computes `Σ Ke·C·r(b)`), keeping the
   method only for nonlinear cases? (Recommended: not now. Revisit after
   Stage 4.)
8. **D14 continuity**: fix macrophyte/vegetation accumulation unconditionally
   (recommended), or behind a compatibility flag?
9. **Non-output diagnostics**: add `output=.FALSE.` to
   `aed_define_sheet_diag_variable` for inter-module arrays like the `ATM_`
   spectra (recommended), or accept 2×`nlambda` extra output variables?
10. **oasim link switch**: `sky_model = 0` meaning "linked from
   `aed_atmosphere`" (recommended; backward compatible), or a separate
   `link_atmosphere` flag?

## Appendix — key references

| Item | Location |
|---|---|
| API env/data types (`I_0`, `rad`, `par`, `extc`, `bioextc`) | `libaed-api/src/aed_api.F90:138-197, 234-290` |
| Coupling type (`link_ext_par`, fractions, `Kw`) | `aed_api.F90:71-99` |
| `aed_run_column` light flow | `aed_api.F90:1616-1703` |
| `Light` (Sep-24 edit at :2308) | `aed_api.F90:2263-2319` |
| Column orientation mapping | `aed_api.F90:1386-1397` |
| Module extinction interface | `libaed-water/src/aed_core.F90:90, 116, 1246`; `aed_common.F90:534` |
| Accumulator reset (D14) | `aed_macrophyte.F90:1551, 1558`; `aed_vegetation.F90:926` |
| Coupling → configure order | `fv_aed.F90:323/325`; `glm_api_aed.F90:225/229`; `schism_aed.F90:458/462`; `elcom_aed.F90:434/437` |
| libaed-light `aed_light` (to retire) | `libaed-light/src/aed_light.F90`; `aed_lighting.F90:41, 75`; `Makefile:41`; `win/libaed-lighting.vfproj:35` |
| GLM API bindings | `GLM/src/glm_api_aed.F90:268, 282, 358, 372` |
| GLM legacy driver light | `GLM/src/glm_aed.F90:580, 721, 1560` |
| TUFLOW-FV WQ object | `libaed-fv/tuflowfv_external_wq/tuflowfv_wq_api.f90:129-146` |
| TUFLOW-FV bindings | `libaed-fv/src/fv_aed.F90:665, 674, 731-781` |
| SCHISM bindings | `schism-aed/src/AED/schism_aed.F90:199, 522-529, 679, 753, 763` |
| ELCOM bindings | `ELCOM/elcom_aed/src/elcom_aed.F90:366, 393, 605, 611, 636-648` |
| ELCOM host light/extinction | `ELCOM/elcom_src_v3/src/elcom_wq.f90:279, 325-330, 410-418, 672-676` |
| oasim spectral optics | `libaed-light/src/aed_oasim.F90:188-280, 416-452, 484-491, 885` |
| `par > 0.95*Io` | `libaed-water/src/aed_phytoplankton.F90:1586` |
