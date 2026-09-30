/******************************************************************************
 * aed_api_guard.c - detection-only guard for the AED API entry points (plan 5.2, S5 2026-09-15)
 *
 * One compare-and-swap word holds the id of the entry point that is active (0 = none), so the flag and its owner
 * change together. Every public entry point of aed_api.F90, aed_ptm.F90 and aed_zones.F90 becomes a thin wrapper:
 *     CALL aed_api_guard_enter(ID) ; CALL <name>_impl(...) ; CALL aed_api_guard_leave(ID)
 * Internal calls between entry points (aed_run_model -> Particles, aed_calculate_particles, aed_split_particles;
 * aed_set_model_data -> set_ptm_aed_var_num) call the _impl routines, so a model step never trips its own guard.
 *
 * Overlap - a second entry while one is active, from another thread or from a callback (aed_ptm_flux_fn) - stops
 * the run: "FATAL: concurrent call to <routine> while <owner> is active". Documented limits: this DETECTS misuse; it
 * does not make interleaved independent model instances in one process safe, and it cannot protect the public
 * module data (ptm_istat, ptm_env, aedZones, ...) or a pointer after it is returned.
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

/* ids follow STEP5_api_inventory_20260915.md; 0 is "no owner" */
static const char *aed_api_names[] = {
    "(none)",
    "aed_set_coupling", "aed_set_model_env", "aed_configure_models", "aed_set_model_data",
    "aed_check_model_setup", "aed_set_mobility_fn", "aed_run_model", "aed_clean_model", "aed_var_index",
    "aed_ptm_init", "Particles", "aed_calculate_particles", "aed_split_particles", "aed_ptm_set_cell_map",
    "set_ptm_aed_var_num", "aed_init_zones", "api_set_zone_funcs"
};
#define AED_API_NIDS ((int)(sizeof aed_api_names / sizeof aed_api_names[0]))

static volatile long aed_api_owner = 0;

static const char *id_name(long id)
{
    return (id > 0 && id < AED_API_NIDS) ? aed_api_names[id] : "(unknown)";
}

/* Fortran: SUBROUTINE aed_api_guard_enter(id) BIND(C, name="aed_api_guard_enter"); INTEGER(C_INT), VALUE :: id */
void aed_api_guard_enter(int id)
{
    long prev = __sync_val_compare_and_swap(&aed_api_owner, 0L, (long)id);
    if ( prev != 0 ) {
        fprintf(stderr, "FATAL: concurrent call to %s while %s is active\n", id_name(id), id_name(prev));
        fflush(stderr);
        exit(1);
    }
}

/* releases only the owner's claim; a mismatched release is itself misuse */
void aed_api_guard_leave(int id)
{
    long prev = __sync_val_compare_and_swap(&aed_api_owner, (long)id, 0L);
    if ( prev != (long)id ) {
        fprintf(stderr, "FATAL: %s released the API guard held by %s\n", id_name(id), id_name(prev));
        fflush(stderr);
        exit(1);
    }
}

/* for tests: the current owner id (0 when free) */
int aed_api_guard_owner(void)
{
    return (int)aed_api_owner;
}
