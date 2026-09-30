      https://aquatic.science.uwa.edu.au/

AED suite is a collection of community-driven libraries of modules and
algorithms for simulation of "aquatic ecodynamics" - water quality, aquatic
biogeochemsitry, biotic habitat and aquatic ecosystem dynamics.

libaed-api is a new simpler interface to the library suite.

Each module aims to be based on state-of-the-art science sourced from a wide
variety of scientific literature, making the library one of the most
advanced available to aquatic ecosystem modellers.

Column and threading contract (libaed-api, 2026-09-11)
-------------------------------------------------------
* aed_run_model processes columns serially in slot order on the calling thread.
  Concurrent calls, or a host OpenMP parallel region around the column loop, are
  UNSUPPORTED: the module work arrays (flux_pel, flux_pel_ptm, flux_atm, flux_ben,
  flux_zon, the internal light arrays, aed_substep_frac) are shared per process, the
  libaed-water module data and the aed_ptm particle stores are process-global, and the
  orphaned `!$OMP DO` directives in the benthic/zone loops have no enclosing PARALLEL
  region in this library. This is a documented limitation, not a verified property.
* Per-column orientation is free: bot_idx <= top_idx or bot_idx > top_idx, mixed
  within one call. The traversal direction is reset on every helper call (before
  2026-09-11 it was an implicitly SAVEd local: one reversed column silenced the pelagic
  kinetics of every column processed after it).
* Before aed_set_coupling the host MUST associate conf%rain_factor, conf%sw_factor,
  conf%friction and conf%Kw (Light() adds Kw unconditionally). Every env array must have
  its own storage of extent >= n_layers; arrays are sectioned, not copied.
* Verified by cases_heat_wall/tools/test_api_ptm_flux.F90 (two columns, 81 checks):
  different rates per column and layer at split factors 1,2,4,8; linear accumulation;
  inactive column untouched (state and diagnostics); swapped slots give identical
  per-column results; direction does not persist; mixed orientation in one call.
