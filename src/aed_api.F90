!###############################################################################
!#                                                                             #
!# aed_api.F90                                                                 #
!#                                                                             #
!# A generic interface between models and libaed-xxx                           #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2024 - 2025 - The University of Western Australia                 #
!#                                                                             #
!#  This file is part of libaed (Library for AquaticEco Dynamics)              #
!#                                                                             #
!#  AED is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  AED is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!###############################################################################

!# CAB: need to check order of indexes in flux,cc etc
!# CAB: libaed core is completely index agnostic so we can do what we want here
!# CAB: but better to follow (vars, layers, columns)
!#
!# CAB: Then need to investigate array of aed columns so we dont have to keep
!# CAB: redefining them.  Fluxes should be OK because dta only persists in a
!# CAB: column run anyway

#include "aed_api.h"

#ifndef __GFORTRAN__
#  ifndef isnan
#    define isnan(x) ieee_is_nan(x)
#  endif
#endif

!# Ultimately these should go - the code should be sufficiently generic
!# to handle all cases, but while merging is happening we probably need them
#define GLM_VARIANT   1
#define TFV_VARIANT   2
#define SCH_VARIANT   3
#define ELC_VARIANT   4

!#define CUR_VARIANT   GLM_VARIANT
#define CUR_VARIANT   SCH_VARIANT


!-------------------------------------------------------------------------------
MODULE aed_api
!
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private
!
!
!-------------------------------------------------------------------------------
!
   PUBLIC aed_configure_models,  &
          aed_coupling_t,        &
          aed_set_coupling,      &
          aed_env_t,             &
          aed_set_model_env,     &
          aed_data_t,            &
          aed_set_model_data,    &
          aed_check_model_setup, &
          aed_mobility_fn_t,     &
          aed_set_mobility_fn,   &
          aed_run_model,         &
          aed_var_index,         &
          aed_clean_model

   !#===========================================================#!
   TYPE AED_DPTR
      AED_REAL,DIMENSION(:),POINTER :: p         => null()
   END TYPE AED_DPTR
   !#===========================================================#!

   !#===========================================================#!
   !* A structure to pass configuration values to AED           *!
   !*-----------------------------------------------------------*!
   TYPE aed_coupling_t
      INTEGER  :: MaxLayers
      INTEGER  :: MaxColumns

      LOGICAL  :: mobility_off
      LOGICAL  :: bioshade_feedback
      LOGICAL  :: repair_state
      LOGICAL  :: link_rain_loss
      LOGICAL  :: link_solar_shade
      LOGICAL  :: link_bottom_drag
      LOGICAL  :: ice

      INTEGER  :: split_factor = 1 !# not sure we need this anymore
      INTEGER  :: benthic_mode = 1

      AED_REAL :: rain_factor  = 1.
      AED_REAL :: sw_factor    = 1.
      AED_REAL :: friction     = 1.

      AED_REAL :: Kw
      AED_REAL :: Ksed

      AED_REAL :: nir_fraction = 0.52   ! 0.51
      AED_REAL :: par_fraction = 0.43   ! 0.45
      AED_REAL :: uva_fraction = 0.048  ! 0.035
      AED_REAL :: uvb_fraction = 0.002  ! 0.005
   END TYPE aed_coupling_t
   !#===========================================================#!

   !#===========================================================#!
   !* A structure to pass data array pointers to AED            *!
   !*-----------------------------------------------------------*!
   TYPE aed_data_t
      AED_REAL,DIMENSION(:,:),POINTER :: cc         => null()  !# (n_vars,n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_hz      => null()  !# (n_ben_vars)
      AED_REAL,DIMENSION(:,:),POINTER :: cc_diag    => null()  !# (n_diag_vars,n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_diag_hz => null()  !# (n_diag_ben_vars)
   END TYPE aed_data_t
   !#===========================================================#!

   !#===========================================================#!
   !* A structure to pass environment array pointers to AED     *!
   !*-----------------------------------------------------------*!
   TYPE aed_env_t
      INTEGER :: n_layers
      LOGICAL, POINTER :: active                    => null()

      AED_REAL,POINTER :: yearday                   => null()
      AED_REAL,POINTER :: timestep                  => null()

      AED_REAL,POINTER :: longitude                 => null()
      AED_REAL,POINTER :: latitude                  => null()

      AED_REAL,DIMENSION(:),POINTER :: height       => null() !# layer height (previously "h")
      AED_REAL,DIMENSION(:),POINTER :: area         => null() !# layer area
      AED_REAL,DIMENSION(:),POINTER :: dz           => null() !# layer thickness
      AED_REAL,DIMENSION(:),POINTER :: depth        => null() !# layer_depth (previously "z")

      AED_REAL,DIMENSION(:),POINTER :: temp         => null() !# temperature
      AED_REAL,DIMENSION(:),POINTER :: salt         => null() !# salinity
      AED_REAL,DIMENSION(:),POINTER :: rho          => null() !# density
      AED_REAL,DIMENSION(:),POINTER :: rad          => null()
      AED_REAL,DIMENSION(:),POINTER :: tss          => null() !# total suspended solids
      AED_REAL,DIMENSION(:),POINTER :: ss1          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss2          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss3          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss4          => null()
      AED_REAL,DIMENSION(:),POINTER :: cvel         => null() !# cell velocity
      AED_REAL,DIMENSION(:),POINTER :: ustar_bed    => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_uorb      => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_t         => null()
      AED_REAL,DIMENSION(:),POINTER :: pres         => null()
      AED_REAL,DIMENSION(:),POINTER :: extc         => null() !# extinction coefficient

      !# sedzones are an odd mix - for GLM a zone will be different at
      !# different levels - while for others a column would be all int one zone
      AED_REAL,DIMENSION(:),POINTER :: sed_zones    => null()
      INTEGER, POINTER :: mat_id       => null()
      AED_REAL,POINTER :: sed_zone     => null()

      AED_REAL,POINTER :: wind         => null()
      AED_REAL,POINTER :: air_temp     => null()
      AED_REAL,POINTER :: air_pres     => null()
      AED_REAL,POINTER :: rain         => null()
      AED_REAL,POINTER :: evap         => null()
      AED_REAL,POINTER :: humidity     => null()
      AED_REAL,POINTER :: longwave     => null()
      AED_REAL,POINTER :: layer_stress => null()
      AED_REAL,POINTER :: col_depth    => null() !# col_depth (sheet - total depth of column)

      AED_REAL,POINTER :: I_0          => null() !# par_sf
      !# if missing the following may be calculated from I_0 (par_sf) above
      AED_REAL,DIMENSION(:),POINTER :: par => null()
      AED_REAL,DIMENSION(:),POINTER :: nir => null()
      AED_REAL,DIMENSION(:),POINTER :: uva => null()
      AED_REAL,DIMENSION(:),POINTER :: uvb => null()

      !# feedback data
      AED_REAL,DIMENSION(:),POINTER :: biodrag    => null()
      AED_REAL,DIMENSION(:),POINTER :: bioextc    => null()
      AED_REAL,DIMENSION(:),POINTER :: solarshade => null()
      AED_REAL,DIMENSION(:),POINTER :: windshade  => null()
      AED_REAL,POINTER :: bathy    => null()
      AED_REAL,POINTER :: rainloss => null()
   END TYPE aed_env_t
   !#===========================================================#!

   !#===========================================================#!
   !* A structure defining a water column.                      *!
   !*-----------------------------------------------------------*!
   TYPE api_col_data_t
      INTEGER :: n_layers           = 0 !# number of layers in this column
                                        !# in cases like GLM this may vary each timestep
      AED_REAL :: col_num

      LOGICAL,POINTER :: active     => null()

      AED_REAL,POINTER :: longitude => null()
      AED_REAL,POINTER :: latitude  => null()

      !# Main arrays storing/pointing to the state and diagnostic variables
      AED_REAL,DIMENSION(:,:),POINTER :: cc         => null()  !# (n_vars, n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_hz      => null()  !# (n_ben_vars)
      AED_REAL,DIMENSION(:,:),POINTER :: cc_diag    => null()  !# (n_diag_vars, n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_diag_hz => null()  !# (n_diag_ben_vars)

      AED_REAL,DIMENSION(:),POINTER :: lheights     => null()  !# (n_layers)
      AED_REAL,DIMENSION(:),POINTER :: dz           => null()
      AED_REAL,DIMENSION(:),POINTER :: area         => null()
      AED_REAL,DIMENSION(:),POINTER :: depth        => null()

      AED_REAL,DIMENSION(:),POINTER :: temp         => null()
      AED_REAL,DIMENSION(:),POINTER :: salt         => null()
      AED_REAL,DIMENSION(:),POINTER :: rho          => null()
      AED_REAL,DIMENSION(:),POINTER :: rad          => null()
      AED_REAL,DIMENSION(:),POINTER :: extc         => null()
      AED_REAL,DIMENSION(:),POINTER :: cvel         => null()  !# cell velocity
      AED_REAL,DIMENSION(:),POINTER :: pres         => null()

      AED_REAL,DIMENSION(:),POINTER :: tss          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss1          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss2          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss3          => null()
      AED_REAL,DIMENSION(:),POINTER :: ss4          => null()

      AED_REAL,DIMENSION(:),POINTER :: ustar_bed    => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_uorb      => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_t         => null()

      AED_REAL,DIMENSION(:),POINTER :: sed_zones    => null()

      !# The following are sheet vars
      INTEGER, POINTER :: mat_id       => null()
      AED_REAL,POINTER :: sed_zone     => null()

      AED_REAL,POINTER :: col_depth    => null()
      AED_REAL,POINTER :: wind         => null()
      AED_REAL,POINTER :: air_temp     => null()
      AED_REAL,POINTER :: air_pres     => null()
      AED_REAL,POINTER :: rain         => null()
      AED_REAL,POINTER :: evap         => null()
      AED_REAL,POINTER :: humidity     => null()
      AED_REAL,POINTER :: longwave     => null()
      AED_REAL,POINTER :: layer_stress => null()

      AED_REAL,POINTER :: I_0          => null()
      !# if missing the following may be calculated from I_0 (par_sf) above
      AED_REAL,DIMENSION(:),POINTER :: par => null()
      AED_REAL,DIMENSION(:),POINTER :: nir => null()
      AED_REAL,DIMENSION(:),POINTER :: uva => null()
      AED_REAL,DIMENSION(:),POINTER :: uvb => null()

      !# These are feedback arrays
      AED_REAL,DIMENSION(:),POINTER :: biodrag    => null()
      AED_REAL,DIMENSION(:),POINTER :: bioextc    => null()
      AED_REAL,DIMENSION(:),POINTER :: solarshade => null()
      AED_REAL,DIMENSION(:),POINTER :: windshade  => null()
      AED_REAL,POINTER :: bathy    => null()
      AED_REAL,POINTER :: rainloss => null()
   END TYPE api_col_data_t
   !#===========================================================#!

   !#===========================================================#!
   !* A structure defining a water column.                      *!
   !*-----------------------------------------------------------*!
   TYPE api_env_def_t
      CHARACTER(:),POINTER :: name
      CHARACTER(:),POINTER :: longname
      CHARACTER(:),POINTER :: units
      LOGICAL :: sheet
      AED_REAL,POINTER :: data      => null()
      AED_REAL,DIMENSION(:),POINTER :: datac => null()
   END TYPE api_env_def_t
   !#===========================================================#!

#define BSSOCIATED(x) ( ASSOCIATED(data(1)%x) )

!
!-------------------------------------------------------------------------------
!MODULE DATA

   CHARACTER(len=80) :: cfg_fname = "none"

   AED_REAL :: Kw, Ksed

   AED_REAL :: rain_factor
   AED_REAL :: sw_factor
   AED_REAL :: friction

   AED_REAL :: nir_fraction =  0.52   ! 0.51
   AED_REAL :: par_fraction =  0.43   ! 0.45
   AED_REAL :: uva_fraction =  0.048  ! 0.035
   AED_REAL :: uvb_fraction =  0.002  ! 0.005

   INTEGER :: benthic_mode = 1

   INTEGER :: split_factor = 1
   LOGICAL :: mobility_off = .FALSE.
   LOGICAL :: bioshade_feedback = .TRUE.
   LOGICAL :: repair_state = .TRUE.
   LOGICAL :: do_zone_averaging = .FALSE.
   LOGICAL :: link_solar_shade = .TRUE.
   LOGICAL :: link_rain_loss = .FALSE.
   LOGICAL :: link_bottom_drag = .FALSE.
!  LOGICAL :: link_surface_drag = .FALSE.
!  LOGICAL :: link_water_density = .FALSE.
   LOGICAL :: link_water_clarity = .FALSE.
   LOGICAL :: link_ext_par = .FALSE.
   LOGICAL :: do_2d_atm_flux = .TRUE.

   LOGICAL :: bottom_one = .TRUE.

!  AED_REAL,DIMENSION(:),POINTER :: sed_zones => null()

   !-------------------------------------------------------------
   !# External variables
   !-------------------------------------------------------------
   TYPE(api_col_data_t),DIMENSION(:),ALLOCATABLE,TARGET :: data

   !# Arrays for environmental variables (used if they are not supplied externally)
!  AED_REAL,DIMENSION(:),POINTER :: nir => null()
!  AED_REAL,DIMENSION(:),POINTER :: par => null()
!  AED_REAL,DIMENSION(:),POINTER :: uva => null()
!  AED_REAL,DIMENSION(:),POINTER :: uvb => null()
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: lpar

!  INTEGER, DIMENSION(:,:),POINTER :: mat_id => null()
!  INTEGER, DIMENSION(:,:),POINTER :: col_num => null()
!  LOGICAL, POINTER :: active => null()

!  !# Maps to nearest cell with water (for riparian exchange)
! may be tuflow specific
   AED_REAL,DIMENSION(:),POINTER :: nearest_active => null()
   AED_REAL,DIMENSION(:),POINTER :: nearest_depth => null()
!  INTEGER, DIMENSION(:),POINTER :: route_table => null()

   !# Maps of surface, bottom and wet/dry (active) cells
!  INTEGER,DIMENSION(:),POINTER :: surf_map => null()
!  INTEGER,DIMENSION(:),POINTER :: benth_map => null()
!  LOGICAL,DIMENSION(:),POINTER :: active => null()

   !# Arrays for work, vertical movement (ws), and cross-boundary fluxes
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: ws   !# (n_vars, n_layers)
   AED_REAL,DIMENSION(:),ALLOCATABLE :: min_, max_

   !# To support light
   AED_REAL,POINTER :: yearday => null()
   AED_REAL,POINTER :: timestep => null()
!# These 2 are column based
!  AED_REAL,POINTER :: longitude => null()
!  AED_REAL,POINTER :: latitude => null()

   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: biodrag
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: bioextc
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: solarshade
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: windshade
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: bathy
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: rainloss

   !# Particle groups
!  INTEGER :: num_groups
!  TYPE(partgroup),DIMENSION(:),POINTER :: particle_groups => null()
!  TYPE(partgroup_cell),DIMENSION(:),ALLOCATABLE :: all_particles => null()

   !# Misc variables/options
   LOGICAL :: request_nearest = .FALSE.
   LOGICAL :: have_nearest = .FALSE.
!  INTEGER :: ThisStep = 0
!  INTEGER :: n_cellids = 0

   !#------------------------------------------------------------
   !#   internal variables
   !#------------------------------------------------------------

   !# Integers storing number of variables being simulated
   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet

   CHARACTER(len=48),ALLOCATABLE :: names(:)
   CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)

   INTEGER :: MaxLayers = 0
   INTEGER :: zone_var = 0

   AED_REAL :: dt_eff
   LOGICAL :: reinited = .FALSE.

#ifdef f2003
   USE, intrinsic :: iso_fortran_env, ONLY : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
#else
#  define stdin  5
#  define stdout 6
#  define stderr 0
#endif
   INTEGER :: log = stderr


  !-----------------------------------------------------------------------------
  INTERFACE

    SUBROUTINE aed_mobility_fn_t(N,dt,h,A,ww,min_C,cc)
       INTEGER,INTENT(in)     :: N       !# number of vertical layers
       AED_REAL,INTENT(in)    :: dt      !# time step (s)
       AED_REAL,INTENT(in)    :: h(*)    !# layer thickness (m)
       AED_REAL,INTENT(in)    :: A(*)    !# layer areas (m2)
       AED_REAL,INTENT(in)    :: ww(*)   !# vertical speed (m/s)
       AED_REAL,INTENT(in)    :: min_C   !# minimum allowed cell concentration
       AED_REAL,INTENT(inout) :: cc(*)   !# cell concentration
    END SUBROUTINE aed_mobility_fn_t

  END INTERFACE
  !-----------------------------------------------------------------------------

  PROCEDURE(aed_mobility_fn_t),POINTER :: doMobility => null()

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_show_vars
!-------------------------------------------------------------------------------
! Print names and attributs of all variables
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: i, j
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   print "(5X,'Configured AED variables to simulate:')"

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%sheet .OR. tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            names(j) = TRIM(tvar%name)
            min_(j) = tvar%minimum
            max_(j) = tvar%maximum
            !print *,"     S(",j,") AED pelagic(3D) variable: ", TRIM(names(j))
            print "(7X,'S(',I4,') water column variable     : ',A)",j , TRIM(names(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%sheet .AND. .NOT. (tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            bennames(j) = TRIM(tvar%name)
            min_(n_vars+j) = tvar%minimum
            max_(n_vars+j) = tvar%maximum
            !print *,"     B(",j,") AED benthic(2D) variable: ", TRIM(bennames(j))
            print "(7X,'B(',I4,') bottom variable           + ',A)",j , TRIM(bennames(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            IF ( .NOT.  tvar%sheet ) THEN
               j = j + 1
               print "(7X,'D(',I4,') water column diagnostic   > ',A)",j , TRIM(tvar%name)
               !print *,"     D(",j,") AED diagnostic 3Dvariable: ", TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            IF (tvar%sheet ) THEN
               j = j + 1
               !print *,"     D(",j,") AED diagnostic 2Dvariable: ", TRIM(tvar%name)
               print "(7X,'D(',I4,') bottom/surface diagnostic ~ ',A)",j , TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE aed_show_vars
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed_configure_models(fname, NumWQ_Vars, NumWQ_Ben, NumWQ_Diag, NumWQ_DiagS)
!-------------------------------------------------------------------------------
! Initialize the AED-API driver by reading settings from "fname".
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: fname
   INTEGER,INTENT(out) :: NumWQ_Vars, NumWQ_Ben
   INTEGER,INTENT(out) :: NumWQ_Diag, NumWQ_DiagS
!
!LOCALS
   INTEGER :: status, i, namlst

#if DEBUG
   TYPE(aed_variable_t),POINTER :: tvar
#endif

   CHARACTER(len=64) :: models(64)
   NAMELIST /aed_models/ models
!
!-------------------------------------------------------------------------------
!BEGIN
#ifdef __INTEL_COMPILER
#  ifdef __INTEL_LLVM_COMPILER
   print *,'    aed_api built using intel fortran (ifx) version ', __INTEL_LLVM_COMPILER
#  else
   print *,'    aed_api built using intel fortran version ', __INTEL_COMPILER
#  endif
#else
# ifdef __GNUC__
   print *,'    aed_api built using gfortran version ', __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
# else
#  ifdef __clang__
    print*,"    aed_api built using flang version ", __clang_major__, '.', __clang_minor__, '.', __clang_patchlevel__
#  else
    print*,"    aed_api built using unknow fortran"
#  endif
# endif
#endif

   print *,'    libaed enabled.... aed_configure_models processing: ', TRIM(fname)
   namlst = find_free_lun()

   print*,'     ---------- AED API config : start ----------'

!  IF ( aed_init_core('.') /= 0 ) STOP "     ERROR: Initialisation of aed_core failed"
   CALL aed_print_version

   !# Create model tree
   print *,"     Processing aed_models config from ",TRIM(fname)
   OPEN(namlst,file=fname,action='read',status='old',iostat=status)
   IF ( status /= 0 ) CALL STOPIT("Cannot open file " // TRIM(fname))

   cfg_fname = fname

   models = ''
   READ(namlst, nml=aed_models, iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed_models"

   DO i=1,size(models)
      IF ( models(i) == '' ) EXIT
      IF ( benthic_mode .GT. 1 ) models(i) = TRIM(models(i)) // ':za' ! make all models zone averaged
      CALL aed_define_model(models(i), namlst)
   ENDDO

   !# should be finished with this file
   CLOSE(namlst)
   print *,"      ... nml file parsing completed."

   n_aed_vars = aed_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

#if DEBUG
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         print *,"AED var ", i, tvar%sheet, tvar%diag, tvar%extern, TRIM(tvar%name)
      ELSE
         print *,"AED var ", i, " is empty"
      ENDIF
   ENDDO
#endif

   print "(/,5X,'AED : n_aed_vars  = ',I3,' ; MaxLayers         = ',I4)",n_aed_vars,MaxLayers
   print "(  5X,'AED : n_vars      = ',I3,' ; n_vars_ben        = ',I4)",n_vars,n_vars_ben
   print "(  5X,'AED : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I4,/)",n_vars_diag,n_vars_diag_sheet

   !# names = grab the names from info
   ALLOCATE(names(n_vars),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (names)'
   ALLOCATE(bennames(n_vars_ben),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (bennames)'

   NumWQ_Vars  = n_vars
   NumWQ_Ben   = n_vars_ben
   NumWQ_Diag  = n_vars_diag
   NumWQ_DiagS = n_vars_diag_sheet

   print*,'     ----------  AED API config : end  ----------'

   aed_configure_models = n_aed_vars
END FUNCTION aed_configure_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_mobility_fn(mobl)
!-------------------------------------------------------------------------------
! The pre-kinetics stage of an AED timestep run will attempt to run
! settling/rising code if mobility_off is false;
! This routine saves the pointer to the mobility code
!-------------------------------------------------------------------------------
!ARGUMENTS
   PROCEDURE(aed_mobility_fn_t),POINTER :: mobl
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
!
    doMobility => mobl
END SUBROUTINE aed_set_mobility_fn
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_coupling(conf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_coupling_t), INTENT(in) :: conf
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   MaxLayers = conf%MaxLayers
   mobility_off = conf%mobility_off
   bioshade_feedback = conf%bioshade_feedback
   repair_state = conf%repair_state
   link_rain_loss = conf%link_rain_loss
   link_solar_shade = conf%link_solar_shade
   link_bottom_drag = conf%link_bottom_drag

   split_factor = conf%split_factor
   if (split_factor == 0) split_factor = 1
   benthic_mode = conf%benthic_mode

   rain_factor = conf%rain_factor
   sw_factor = conf%sw_factor
   friction = conf%friction

   Kw = conf%Kw
   Ksed = conf%Ksed
END SUBROUTINE aed_set_coupling
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_data(dat, ncols, nlevs)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: ncols, nlevs
   TYPE(aed_data_t),INTENT(in) :: dat(ncols)
!
!LOCALS
   INTEGER :: av, v, sv, status, col
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (.NOT. ALLOCATED(data) ) ALLOCATE(data(ncols))
   DO col=1,ncols
      data(col)%cc         => dat(col)%cc
!# Eventually this will be properly separated
!     data(col)%cc_hz      => dat(col)%cc_hz
      data(col)%cc_hz      => dat(col)%cc(n_vars:n_vars+n_vars_ben,1)
      data(col)%cc_diag    => dat(col)%cc_diag
      data(col)%cc_diag_hz => dat(col)%cc_diag_hz
   ENDDO
   IF (nlevs > MaxLayers) MaxLayers = nlevs

   ALLOCATE(min_((n_vars + n_vars_ben))) ; ALLOCATE(max_((n_vars + n_vars_ben)))

   CALL aed_show_vars

   !----------------------------------------------------------------------------

   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "     ERROR getting variable info"
      IF ( .NOT. ( tvar%extern .OR. tvar%diag) ) THEN  !# neither global nor diagnostic variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            DO col=1,ncols
               data(col)%cc(n_vars+sv,:) = tvar%initial
            ENDDO
         ELSE
            v = v + 1
            DO col=1,ncols
               data(col)%cc(v,:) = tvar%initial
            ENDDO
         ENDIF
      ENDIF
   ENDDO

   !# Allocate array with vertical movement rates (m/s, positive for upwards),
   !# and set these to the values provided by the model.
   !# allocated for all vars even though only state vars entries will be used
   ALLOCATE(ws(n_aed_vars, MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = zero_

!  !# Trigger an error if we don't have all we need.
!  CALL aed_check_model_setup
END SUBROUTINE aed_set_model_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_env(env, ncols)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: ncols
   TYPE(aed_env_t),INTENT(in) :: env(ncols)
!
!LOCALS
   INTEGER :: tv, col
   LOGICAL :: need_biodg = .FALSE., need_bioex = .FALSE.
   LOGICAL :: need_sshad = .FALSE., need_wshad = .FALSE.
   LOGICAL :: need_rianl = .FALSE., need_bathy = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (.NOT. ALLOCATED(data) ) ALLOCATE(data(ncols))

   IF (.NOT. link_ext_par) ALLOCATE(lpar(MaxLayers,ncols))

!  dt_eff = timestep/FLOAT(split_factor)
   dt_eff = env(1)%timestep/FLOAT(split_factor)

   !# feedback
   need_biodg = (.NOT.ASSOCIATED(env(1)%biodrag))
   need_bioex = (.NOT.ASSOCIATED(env(1)%bioextc))
   need_sshad = (.NOT.ASSOCIATED(env(1)%solarshade))
   need_wshad = (.NOT.ASSOCIATED(env(1)%windshade))
   need_rianl = (.NOT.ASSOCIATED(env(1)%rainloss))
   need_bathy = (.NOT.ASSOCIATED(env(1)%bathy))

   IF (need_biodg) THEN ; ALLOCATE(biodrag(MaxLayers,ncols))    ; biodrag = zero_    ; ENDIF
   IF (need_bioex) THEN ; ALLOCATE(bioextc(MaxLayers,ncols))    ; bioextc = zero_    ; ENDIF
   IF (need_sshad) THEN ; ALLOCATE(solarshade(MaxLayers,ncols)) ; solarshade = zero_ ; ENDIF
   IF (need_wshad) THEN ; ALLOCATE(windshade(MaxLayers,ncols))  ; windshade = zero_  ; ENDIF
   IF (need_rianl) THEN ; ALLOCATE(rainloss(ncols))             ; rainloss = zero_   ; ENDIF
   IF (need_bathy) THEN ; ALLOCATE(bathy(ncols))                ; bathy = zero_      ; ENDIF

   DO col=1,ncols
      yearday  => env(col)%yearday
      timestep => env(col)%timestep

      data(col)%n_layers     =  env(col)%n_layers
      data(col)%col_num      =  col

      data(col)%longitude    => env(col)%longitude
      data(col)%latitude     => env(col)%latitude

      data(col)%temp         => env(col)%temp
      data(col)%salt         => env(col)%salt
      data(col)%rho          => env(col)%rho
      data(col)%dz           => env(col)%dz
      data(col)%lheights     => env(col)%height
      data(col)%area         => env(col)%area
      data(col)%depth        => env(col)%depth
      data(col)%extc         => env(col)%extc
      data(col)%tss          => env(col)%tss
      data(col)%ss1          => env(col)%ss1
      data(col)%ss2          => env(col)%ss2
      data(col)%ss3          => env(col)%ss3
      data(col)%ss4          => env(col)%ss4
      data(col)%cvel         => env(col)%cvel
      data(col)%rad          => env(col)%rad
      data(col)%ustar_bed    => env(col)%ustar_bed
      data(col)%wv_uorb      => env(col)%wv_uorb
      data(col)%wv_t         => env(col)%wv_t
      data(col)%pres         => env(col)%pres

      data(col)%I_0          => env(col)%I_0
      data(col)%wind         => env(col)%wind
      data(col)%air_temp     => env(col)%air_temp
      data(col)%air_pres     => env(col)%air_pres
      data(col)%rain         => env(col)%rain
      data(col)%evap         => env(col)%evap
      data(col)%humidity     => env(col)%humidity
      data(col)%longwave     => env(col)%longwave
      data(col)%col_depth    => env(col)%col_depth
      data(col)%layer_stress => env(col)%layer_stress
      data(col)%sed_zones    => env(col)%sed_zones
      data(col)%sed_zone     => env(col)%sed_zone

      IF (need_biodg) THEN ; data(col)%biodrag    => biodrag(:,col)
      ELSE ; data(col)%biodrag    => env(col)%biodrag    ; ENDIF
      IF (need_bioex) THEN ; data(col)%bioextc    => bioextc(:,col)
      ELSE ; data(col)%bioextc    => env(col)%bioextc    ; ENDIF
      IF (need_sshad) THEN ; data(col)%solarshade => solarshade(:,col)
      ELSE ; data(col)%solarshade => env(col)%solarshade ; ENDIF
      IF (need_wshad) THEN ; data(col)%windshade  => windshade(:,col)
      ELSE ; data(col)%windshade  => env(col)%windshade  ; ENDIF
      IF (need_rianl) THEN ; data(col)%rainloss   => rainloss(col)
      ELSE ; data(col)%rainloss   => env(col)%rainloss   ; ENDIF
      IF (need_bathy) THEN ; data(col)%bathy      => bathy(col)
      ELSE ; data(col)%bathy      => env(col)%bathy      ; ENDIF

!     data(col)%par          => env(col)%par
      IF (link_ext_par) THEN
        data(col)%par        => env(col)%par
      ELSE
        data(col)%par        => lpar(:,col)
      ENDIF
      data(col)%nir          => env(col)%nir
      data(col)%uva          => env(col)%uva
      data(col)%uvb          => env(col)%uvb

      data(col)%mat_id       => env(col)%mat_id
      data(col)%active       => env(col)%active
   ENDDO

   IF (ASSOCIATED(yearday))        tv=aed_provide_sheet_global('yearday',       'yearday',           'day'           )
   IF (ASSOCIATED(timestep))       tv=aed_provide_sheet_global('timestep',      'timestep',          'seconds'       )

   IF (BSSOCIATED(longitude))      tv=aed_provide_sheet_global('longitude',     'longitude',         'radians'       )
   IF (BSSOCIATED(latitude))       tv=aed_provide_sheet_global('latitude',      'latitude',          'radians'       )

   IF (BSSOCIATED(temp))           tv=aed_provide_global('temperature','temperature',           'celsius')
   IF (BSSOCIATED(salt))           tv=aed_provide_global('salinity',   'salinity',              'g/kg'   )
   IF (BSSOCIATED(rho))            tv=aed_provide_global('density',    'density',               'kg/m3'  )
   IF (BSSOCIATED(dz))             tv=aed_provide_global('layer_ht',   'layer heights',         'm'      )
   IF (BSSOCIATED(area))           tv=aed_provide_global('layer_area', 'layer area',            'm2'     )
   IF (BSSOCIATED(depth))          tv=aed_provide_global('depth',      'depth',                 'm'      )
   IF (BSSOCIATED(extc))           tv=aed_provide_global('extc_coef',  'extinction coefficient','/m'     )
   IF (BSSOCIATED(tss))            tv=aed_provide_global('tss',        'tss',                   'g/m3'   )
   IF (BSSOCIATED(ss1))            tv=aed_provide_global('ss1',        'ss1',                   'g/m3'   )
   IF (BSSOCIATED(ss2))            tv=aed_provide_global('ss2',        'ss2',                   'g/m3'   )
   IF (BSSOCIATED(ss3))            tv=aed_provide_global('ss3',        'ss3',                   'g/m3'   )
   IF (BSSOCIATED(ss4))            tv=aed_provide_global('ss4',        'ss4',                   'g/m3'   )
   IF (BSSOCIATED(cvel))           tv=aed_provide_global('cell_vel',   'cell velocity',         'm/s'    )
   IF (BSSOCIATED(pres))           tv=aed_provide_global('pressure',   'pressure',              ''       )

   IF (BSSOCIATED(sed_zones))      tv=aed_provide_global('sed_zones',  'sediment zones',        '-'      )

                                   tv=aed_provide_global('nir',                 'nir',               'W/m2'   )
                                   tv=aed_provide_global('par',                 'par',               'W/m2'   )
                                   tv=aed_provide_global('uva',                 'uva',               'W/m2'   )
                                   tv=aed_provide_global('uvb',                 'uvb',               'W/m2'   )
   IF (BSSOCIATED(I_0))            tv=aed_provide_sheet_global('par_sf',        'par_sf',            'W/m2'   )

   IF (BSSOCIATED(col_depth))      tv=aed_provide_sheet_global('col_depth',     'column water depth','m above bottom')
   IF (BSSOCIATED(wind))           tv=aed_provide_sheet_global('wind_speed',    'wind speed',        'm/s'           )
   IF (BSSOCIATED(air_temp))       tv=aed_provide_sheet_global('air_temp',      'air temperature',   'celsius'       )
   IF (BSSOCIATED(air_pres))       tv=aed_provide_sheet_global('air_pres',      'air pressure',      'Pa'            )
   IF (BSSOCIATED(rain))           tv=aed_provide_sheet_global('rain',          'rainfall',          'm/s'           )
   IF (BSSOCIATED(evap))           tv=aed_provide_sheet_global('evap',          'evaporation',       'm/s'           )
   IF (BSSOCIATED(humidity))       tv=aed_provide_sheet_global('humidity',      'relative humidity', '-'             )
   IF (BSSOCIATED(longwave))       tv=aed_provide_sheet_global('longwave',      'longwave',          'W/m2'          )
   IF (BSSOCIATED(layer_stress))   tv=aed_provide_sheet_global('taub',          'layer stress',      'N/m2'          )

!  IF (ASSOCIATED(col_num))        tv=aed_provide_sheet_global('col_num',       'column number',     '-'             )
                                   tv=aed_provide_sheet_global('col_num',       'column number',     '-'             )
   IF (ASSOCIATED(nearest_active)) tv=aed_provide_sheet_global('nearest_active','nearest active',    '-'             )
   IF (ASSOCIATED(nearest_depth))  tv=aed_provide_sheet_global('nearest_depth', 'nearest depth',     'm'             )

   IF (BSSOCIATED(mat_id))         tv=aed_provide_sheet_global('material',      'material',          '-'             )

   IF (BSSOCIATED(sed_zone))       tv=aed_provide_sheet_global('sed_zone',      'current sediment zone', '-'         )

                                   tv=aed_provide_sheet_global('biodrag',       'biodrag',           ''              )
                                   tv=aed_provide_sheet_global('bioextc',       'bioextc',           ''              )
                                   tv=aed_provide_sheet_global('solarshade',    'solarshade',        ''              )
                                   tv=aed_provide_sheet_global('windshade',     'windshade',         ''              )
                                   tv=aed_provide_sheet_global('bathy',         'bathy',             'm above datum' )
                                   tv=aed_provide_sheet_global('rainloss',      'rain loss',         'm/s'           )
END SUBROUTINE aed_set_model_env
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_check_model_setup
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: av
   INTEGER :: v, d, sv, sd, ev, err_count
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed_vars
      IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = BSSOCIATED(temp)
            CASE ( 'salinity' )    ; tvar%found = BSSOCIATED(salt)
            CASE ( 'density' )     ; tvar%found = BSSOCIATED(rho)
            CASE ( 'layer_ht' )    ; tvar%found = BSSOCIATED(dz)
            CASE ( 'extc_coef' )   ; tvar%found = BSSOCIATED(extc)
            CASE ( 'tss' )         ; tvar%found = BSSOCIATED(tss)
            CASE ( 'cell_vel' )    ; tvar%found = BSSOCIATED(cvel)
            CASE ( 'par' )         ; tvar%found = BSSOCIATED(par)
            CASE ( 'nir' )         ; tvar%found = BSSOCIATED(nir)
            CASE ( 'uva' )         ; tvar%found = BSSOCIATED(uva)
            CASE ( 'uvb' )         ; tvar%found = BSSOCIATED(uvb)
            CASE ( 'pressure' )    ; tvar%found = BSSOCIATED(pres)
            CASE ( 'depth' )       ; tvar%found = BSSOCIATED(depth)
            CASE ( 'sed_zone' )    ; tvar%found = BSSOCIATED(sed_zone)
            CASE ( 'sed_zones' )   ; tvar%found = BSSOCIATED(sed_zones)
            CASE ( 'wind_speed' )  ; tvar%found = BSSOCIATED(wind)
            CASE ( 'par_sf' )      ; tvar%found = BSSOCIATED(I_0)
            CASE ( 'taub' )        ; tvar%found = BSSOCIATED(layer_stress)
            CASE ( 'col_depth' )   ; tvar%found = BSSOCIATED(col_depth)
            CASE ( 'layer_area' )  ; tvar%found = BSSOCIATED(area)
            CASE ( 'rain' )        ; tvar%found = BSSOCIATED(rain)
            CASE ( 'evap' )        ; tvar%found = BSSOCIATED(evap)
            CASE ( 'air_temp' )    ; tvar%found = BSSOCIATED(air_temp)
            CASE ( 'air_pres' )    ; tvar%found = BSSOCIATED(air_pres)
            CASE ( 'humidity' )    ; tvar%found = BSSOCIATED(humidity)
            CASE ( 'longwave' )    ; tvar%found = BSSOCIATED(longwave)
            CASE ( 'longitude' )   ; tvar%found = BSSOCIATED(longitude)
            CASE ( 'latitude' )    ; tvar%found = BSSOCIATED(latitude)
            CASE ( 'yearday' )     ; tvar%found = ASSOCIATED(yearday)
            CASE ( 'timestep' )    ; tvar%found = ASSOCIATED(timestep)

            CASE ( 'rainloss' )    ; tvar%found = BSSOCIATED(rainloss)
            CASE ( 'material' )    ; tvar%found = BSSOCIATED(mat_id)
            CASE ( 'bathy' )       ; tvar%found = BSSOCIATED(bathy)
            CASE ( 'ss1' )         ; tvar%found = BSSOCIATED(ss1)
            CASE ( 'ss2' )         ; tvar%found = BSSOCIATED(ss2)
            CASE ( 'ss3' )         ; tvar%found = BSSOCIATED(ss3)
            CASE ( 'ss4' )         ; tvar%found = BSSOCIATED(ss4)
!           CASE ( 'col_num' )     ; tvar%found = BSSOCIATED(col_num)
            CASE ( 'col_num' )     ; tvar%found = .TRUE.
            CASE ( 'nearest_active' ) ; tvar%found = have_nearest ; request_nearest = have_nearest
            CASE ( 'nearest_depth' )  ; tvar%found = have_nearest ; request_nearest = have_nearest

            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
         ELSE
            d = d + 1
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
         ELSE
            v = v + 1
         ENDIF
      ENDIF
      IF ( .NOT. tvar%found ) THEN
         print *, "ERROR: Undefined variable ", TRIM(tvar%name)
         err_count = err_count + 1
      ENDIF
   ENDDO

   IF ( n_vars < v ) print *,"More vars than expected",v,n_vars
   IF ( n_vars_ben < sv ) print *,"More sheet vars than expected"
   IF ( n_vars_diag < d ) print *,"More diag vars than expected"
   IF ( n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   IF ( err_count > 0 ) CALL STOPIT("*** Errors in configuration")
END SUBROUTINE aed_check_model_setup
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, col, top, flux_pel, flux_atm, flux_ben, flux_rip)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: col, top
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_rip(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed_vars
      IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE ( TRIM(tvar%name) )
            CASE ( 'temperature' ) ; column(av)%cell => data(col)%temp
            CASE ( 'salinity' )    ; column(av)%cell => data(col)%salt
            CASE ( 'density' )     ; column(av)%cell => data(col)%rho
            CASE ( 'layer_ht' )    ; column(av)%cell => data(col)%dz
            CASE ( 'layer_area' )  ; column(av)%cell => data(col)%area
            CASE ( 'depth' )       ; column(av)%cell => data(col)%depth
            CASE ( 'pressure' )    ; column(av)%cell => data(col)%pres
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => data(col)%col_depth

            CASE ( 'extc_coef' )   ; column(av)%cell => data(col)%extc
            CASE ( 'tss' )         ; column(av)%cell => data(col)%tss
            CASE ( 'ss1' )         ; column(av)%cell => data(col)%tss   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss2' )         ; column(av)%cell => data(col)%tss   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss3' )         ; column(av)%cell => data(col)%tss   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss4' )         ; column(av)%cell => data(col)%tss   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'cell_vel' )    ; column(av)%cell => data(col)%cvel

            CASE ( 'par' )         ; column(av)%cell => data(col)%par
            CASE ( 'nir' )         ; column(av)%cell => data(col)%nir
            CASE ( 'uva' )         ; column(av)%cell => data(col)%uva
            CASE ( 'uvb' )         ; column(av)%cell => data(col)%uvb

            CASE ( 'par_sf' )      ; column(av)%cell_sheet => data(col)%I_0
            CASE ( 'longwave' )    ; column(av)%cell_sheet => data(col)%longwave

            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => data(col)%wind
            CASE ( 'rain' )        ; column(av)%cell_sheet => data(col)%rain
            CASE ( 'rainloss' )    ; column(av)%cell_sheet => data(col)%rainloss
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => data(col)%air_temp
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => data(col)%air_pres
            CASE ( 'humidity' )    ; column(av)%cell_sheet => data(col)%humidity
            CASE ( 'evap' )        ; column(av)%cell_sheet => data(col)%evap

            CASE ( 'taub' )        ; column(av)%cell_sheet => data(col)%layer_stress ! CAB? col_taub

            CASE ( 'sed_zones' )   ; column(av)%cell => data(col)%sed_zones
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => data(col)%sed_zone

            CASE ( 'bathy' )       ; column(av)%cell_sheet => data(col)%bathy
            CASE ( 'col_num' )     ; column(av)%cell_sheet => data(col)%col_num

            CASE ( 'nearest_active' ) ; column(av)%cell_sheet => nearest_active(col)
            CASE ( 'nearest_depth' )  ; column(av)%cell_sheet => nearest_depth(col)

     ! CAB: Not handled yet
     !      CASE ( 'material' )    ; IF ( do_zone_averaging ) THEN
     !                                  column(av)%cell_sheet => zone(zm(col))
     !                               ELSE
     !                                  column(av)%cell_sheet => mat(col)
     !                               ENDIF
!           CASE ( 'material' )    ; column(av)%cell_sheet => data(col)%mat_id

            CASE ( 'longitude' )   ; column(av)%cell_sheet => data(col)%longitude
            CASE ( 'latitude' )    ; column(av)%cell_sheet => data(col)%latitude
            CASE ( 'yearday' )     ; column(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column(av)%cell_sheet => timestep

            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd =    sd + 1
            column(av)%cell_sheet => data(col)%cc_diag_hz(sd)
         ELSE
            d = d + 1
            column(av)%cell => data(col)%cc_diag(d,:)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            column(av)%cell_sheet => data(col)%cc_hz(sv)
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
            column(av)%flux_rip => flux_rip(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => data(col)%cc(v,:)
            column(av)%flux_pel => flux_pel(v,:)
            column(av)%flux_ben => flux_ben(v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_rip => flux_rip(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_states(column, col, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: col, wlev
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv
   INTEGER :: i,v,lev
#if DEBUG
   INTEGER :: last_naned
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
#if DEBUG
   last_naned = -1
#endif
   DO lev=1, wlev
      CALL aed_equilibrate(column, lev)    !MH this should be in the main do_glm routine ????!!!
      v = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               IF ( repair_state ) THEN
#if DEBUG
                  IF ( isnan(data(col)%cc(v,lev)) ) last_naned = i
#endif
                  IF ( .NOT. isnan(min_(v)) ) THEN
                     IF ( data(col)%cc(v, lev) < min_(v) ) data(col)%cc(v, lev) = min_(v)
                  ENDIF
                  IF ( .NOT. isnan(max_(v)) ) THEN
                     IF ( data(col)%cc(v, lev) > max_(v) ) data(col)%cc(v, lev) = max_(v)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO

#if DEBUG
   IF ( last_naned > -1 ) THEN
      print*
      IF ( aed_get_var(last_naned, tv) ) THEN
         print*,"NaNs detected in CC in var ", TRIM(tv%name)
      ELSE
         print*,"NaNs detected in CC unidentified var"
      ENDIF
!     STOP
   ENDIF
#endif
END SUBROUTINE check_states
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_run_model(nCols, wlev, doSurface)
!-------------------------------------------------------------------------------
!                        wlev is the number of levels used;
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nCols
   INTEGER,INTENT(in) :: wlev
   LOGICAL,INTENT(in) :: doSurface
!
!LOCALS
   INTEGER :: col, top, bot

   TYPE (aed_column_t),ALLOCATABLE,SAVE :: column(:)
   TYPE (aed_column_t),ALLOCATABLE,SAVE :: column_sed(:)
   AED_REAL,ALLOCATABLE,TARGET,SAVE :: flux_ben(:)
   AED_REAL,ALLOCATABLE,TARGET,SAVE :: flux_atm(:)
   AED_REAL,ALLOCATABLE,TARGET,SAVE :: flux_rip(:)
   AED_REAL,ALLOCATABLE,TARGET,SAVE :: flux_pel(:,:)
   AED_REAL,ALLOCATABLE,TARGET,SAVE :: flux_zon(:,:)
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: flux_pel_pre
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: flux_pel_z
!
!-------------------------------------------------------------------------------
!BEGIN
!# This worked for GLM because its zones were detirmined by the layers - this is not the case for
!# most other systems so it needs a rethink [same calc used further down too)
   IF (.NOT.ALLOCATED(flux_pel))   ALLOCATE(flux_pel(n_vars+n_vars_ben,MAX(wlev,aed_n_zones)))

   IF (.NOT.ALLOCATED(flux_zon))   ALLOCATE(flux_zon(n_vars+n_vars_ben, aed_n_zones))
   IF (.NOT.ALLOCATED(flux_ben))   ALLOCATE(flux_ben(n_vars+n_vars_ben))
   IF (.NOT.ALLOCATED(flux_atm))   ALLOCATE(flux_atm(n_vars+n_vars_ben))
   IF (.NOT.ALLOCATED(flux_rip))   ALLOCATE(flux_rip(n_vars+n_vars_ben))
   IF (.NOT.ALLOCATED(column))     ALLOCATE(column(n_aed_vars))
   IF (.NOT.ALLOCATED(column_sed)) ALLOCATE(column_sed(n_aed_vars))

   IF (.NOT.ALLOCATED(flux_pel_pre)) ALLOCATE(flux_pel_pre(n_vars+n_vars_ben, MAX(wlev, aed_n_zones)))
   IF (.NOT.ALLOCATED(flux_pel_z))   ALLOCATE(flux_pel_z(n_vars+n_vars_ben, MAX(wlev, aed_n_zones)))

   IF ( bottom_one ) THEN
      top = wlev ; bot = 1
   ELSE
      top = 1 ; bot = wlev
   ENDIF

   DO col=1, nCols
      data(col)%cc_diag = 0.
      data(col)%cc_diag_hz = 0.
   ENDDO

   IF ( .NOT. reinited ) THEN
      DO col=1, nCols
      !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
         CALL define_column(column, col, top, flux_pel, flux_atm, flux_ben, flux_rip)
         CALL re_initialize(column)
      ENDDO
      reinited = .TRUE.
   ENDIF

   DO col=1, nCols
   !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
      CALL define_column(column, col, top, flux_pel, flux_atm, flux_ben, flux_rip)
      CALL pre_kinetics(col)
   ENDDO

   DO col=1, nCols
   !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
      CALL define_column(column, col, top, flux_pel, flux_atm, flux_ben, flux_rip)
      CALL aed_run_column(wlev, doSurface)
   ENDDO

!-------------------------------------------------------------------------------
CONTAINS


   !############################################################################
   SUBROUTINE aed_run_column(wlev, doSurface)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in) :: wlev
      LOGICAL,INTENT(in) :: doSurface
   !
   !LOCALS
      INTEGER  :: v, lev, zon, split
   !
   !----------------------------------------------------------------------------
   !BEGIN
      DO split=1,split_factor
         IF (benthic_mode .GT. 1) THEN
            CALL p_calc_zone_areas(aedZones, aed_n_zones, data(col)%area, data(col)%lheights, wlev)
            CALL p_copy_to_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, wlev)
         ENDIF

         !# Time-integrate one biological time step
         CALL calculate_fluxes(column, wlev)

         !# Update the water column layers
         DO v = 1, n_vars
            DO lev = 1, wlev
               data(col)%cc(v, lev) = data(col)%cc(v, lev) + dt_eff*flux_pel(v, lev)
            ENDDO
         ENDDO

         !# Now update benthic variables, depending on whether zones are simulated
         IF ( benthic_mode .GT. 1 ) THEN
            ! Loop through benthic state variables to update their mass
            DO v = n_vars+1, n_vars+n_vars_ben
               ! Loop through each sediment zone
               DO zon = 1, aed_n_zones
                  ! Update the main cc_sed data array with the
                  aedZones(zon)%z_cc(v, bot) = aedZones(zon)%z_cc(v, bot) + dt_eff*flux_zon(v, zon)
               ENDDO
            ENDDO

            !# Distribute cc-sed benthic properties back into main cc array
            CALL p_copy_from_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, wlev)
         ELSE
            DO v = n_vars+1, n_vars+n_vars_ben
               data(col)%cc(v, bot) = data(col)%cc(v, bot) + dt_eff*flux_ben(v)
            ENDDO
         ENDIF

         CALL check_states(column, col, wlev)
      ENDDO
   END SUBROUTINE aed_run_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE pre_kinetics(col)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in) :: col
   !
   !LOCALS
      TYPE(aed_variable_t),POINTER :: tv
      AED_REAL :: min_C
      INTEGER  :: i, lev, v
   !
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( .NOT. mobility_off ) THEN
         ! compute vertical settling/mobility
         v = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tv) ) THEN
               IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
                  v = v + 1
                  ! only for state_vars that are not sheet
                  IF ( .NOT. ieee_is_nan(tv%mobility) ) THEN
                     ! default to ws that was set during initialisation
                     ws(i,:) = tv%mobility
                  ELSE
                     ! zero nan values
                     ws(i,:) = zero_
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         DO lev = 1, wlev
            ! update ws for modules that use the mobility method
            !# direction doesn't seem to matter ? CAB
            CALL aed_mobility(column, lev, ws(:,lev))
         ENDDO

         !# Calculate source/sink terms due to the settling or rising of
         !# state variables in the water column (note that settling into benthos
         !# is done in aed_do_benthos)
         IF ( ASSOCIATED(doMobility) ) THEN
            v = 0
            DO i=1,n_aed_vars
               IF ( aed_get_var(i, tv) ) THEN
                  IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
                     v = v + 1
                     !# only for state_vars that are not sheet, and also non-zero ws
                     IF ( .NOT. isnan(tv%mobility) .AND. SUM(ABS(ws(i,:)))>zero_ ) THEN
                        min_C = tv%minimum
                        CALL doMobility(wlev, timestep, data(col)%dz, data(col)%area, ws(i,:), min_C, data(col)%cc(v, :))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      CALL check_states(column, col, wlev)

#if CUR_VARIANT == GLM_VARIANT
!CAB ??
      IF (benthic_mode .GT. 1) THEN
         CALL p_calc_zone_areas(aedZones, aed_n_zones, data(col)%area, data(col)%lheights, wlev)
         CALL p_copy_to_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, wlev)
      ENDIF
#endif

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer

      !# populate local light/extc arrays one column at a time
      IF (.NOT. link_ext_par) THEN
         IF ( bottom_one ) THEN
            CALL update_light(column, col, top)
         ELSE
            CALL update_light(column, col, bot-top+1)
         ENDIF
      ENDIF

      ! non PAR bandwidth fractions (set assuming single light extinction)
      data(col)%nir = (data(col)%par/par_fraction) * nir_fraction
      data(col)%uva = (data(col)%par/par_fraction) * uva_fraction
      data(col)%uvb = (data(col)%par/par_fraction) * uvb_fraction
   END SUBROUTINE pre_kinetics
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE define_zone_column(column, zon, top, bot)
   !----------------------------------------------------------------------------
   ! Set up the current column pointers
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: column(n_aed_vars)
      INTEGER, INTENT(in)  :: zon, top, bot
   !
   !LOCALS
      INTEGER :: av
      INTEGER :: v, d, sv, sd, ev
      TYPE(aed_variable_t),POINTER :: tvar
   !
   !----------------------------------------------------------------------------
   !BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE ( TRIM(tvar%name) )
            CASE ( 'temperature' ) ; column(av)%cell => aedZones(zon)%z_env(:)%z_temp
            CASE ( 'salinity' )    ; column(av)%cell => aedZones(zon)%z_env(:)%z_salt
            CASE ( 'density' )     ; column(av)%cell => aedZones(zon)%z_env(:)%z_rho
            CASE ( 'layer_ht' )    ; column(av)%cell => aedZones(zon)%z_env(:)%z_dz
            CASE ( 'extc_coef' )   ; column(av)%cell => aedZones(zon)%z_env(:)%z_extc
            CASE ( 'tss' )         ; column(av)%cell => aedZones(zon)%z_env(:)%z_tss
            CASE ( 'cell_vel' )    ; column(av)%cell => aedZones(zon)%z_env(:)%z_vel
            CASE ( 'par' )         ; column(av)%cell => aedZones(zon)%z_env(:)%z_par
            CASE ( 'nir' )         ; column(av)%cell => aedZones(zon)%z_env(:)%z_nir
            CASE ( 'uva' )         ; column(av)%cell => aedZones(zon)%z_env(:)%z_uva
            CASE ( 'uvb' )         ; column(av)%cell => aedZones(zon)%z_env(:)%z_uvb
            CASE ( 'pressure' )    ; column(av)%cell => aedZones(zon)%z_env(:)%z_pres
            CASE ( 'depth' )       ; column(av)%cell => aedZones(zon)%z_env(:)%z_depth
            CASE ( 'layer_area' )  ; column(av)%cell => aedZones(zon)%z_env(:)%z_area
            CASE ( 'sed_zones' )   ; column(av)%cell => aedZones(zon)%z_env(:)%z_sed_zones; zone_var = av
        !   CASE ( 'sed_zone' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(bot)%z_sed_zones; zone_var = av !CAB ??? (bot)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_sed_zone; zone_var = av
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_wind
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_I_0
            CASE ( 'taub' )        ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_layer_stress !CAB ??? (bot)
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_col_depth
            CASE ( 'rain' )        ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_rain
            CASE ( 'evap' )        ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_evap
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_air_temp
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_air_pres
            CASE ( 'humidity' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_humidity
            CASE ( 'longitude' )   ; column(av)%cell_sheet => aedZones(zon)%longitude
            CASE ( 'latitude' )    ; column(av)%cell_sheet => aedZones(zon)%latitude
            CASE ( 'yearday' )     ; column(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column(av)%cell_sheet => timestep
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
         ELSE
            d = d + 1
            column(av)%cell => aedZones(zon)%z_cc_diag(d, :)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
!CAB ??
#if 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(n_vars+sv, bot)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(n_vars+sv, top)
            ENDIF
#else
            column(av)%cell_sheet => aedZones(zon)%z_cc_hz(sv)
#endif
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => aedZones(zon)%z_cc(v, :)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(v, :)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
   END SUBROUTINE define_zone_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE re_initialize(column)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: column(:)
   !
   !LOCALS
      INTEGER :: lev,zon,av,sv,sd
      TYPE(aed_variable_t),POINTER :: tvar
   !
   !----------------------------------------------------------------------------
   !BEGIN
      DO lev=1, wlev
         CALL aed_initialize(column, lev)
      ENDDO

#if CUR_VARIANT == GLM_VARIANT
      !# (1) BENTHIC INITIALISATION
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers

         DO zon=1,aed_n_zones
            CALL define_zone_column(column_sed, zon, aed_n_zones, 1)

            aedZones(zon)%z_env(1)%z_sed_zones = zon  !MH TMP !CAB ???
!           print *,'aedZones(zon)%z_sed_zones',aedZones(zon)%z_sed_zones !MH TMP
            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_env(1)%z_sed_zones !CAB ???

            sv = 0 ; sd = 0

            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%diag ) THEN  !# Diagnostic variable
                  IF ( tvar%sheet ) THEN
                     sd = sd + 1
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
                  ENDIF
               ELSEIF ( .NOT. tvar%extern ) THEN !# State variable
                  IF ( tvar%sheet ) THEN
                     sv = sv + 1
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc_hz(sv)
                  ENDIF
               ENDIF
            ENDDO

            CALL aed_initialize_benthic(column_sed, zon)
         ENDDO
      ENDIF
#else
!# assumes only GLM does the odd zoning stuff
      IF ( .NOT. do_zone_averaging ) THEN
         CALL aed_initialize_benthic(column, 1)
  ! CAB: Not yet
  !   ELSE
  !      CALL aed_initialize_zone_benthic(nCols, active, n_aed_vars, data(col)%cc_diag, benth_map)
      ENDIF
#endif
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE glm_benthics(column, wlev)
   !----------------------------------------------------------------------------
   ! Calculate the benthic fluxes for GLM. This has been seperated from the
   ! general calc because its more complicated.
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: column(:)
      INTEGER, INTENT(in) :: wlev
   !
   !LOCALS
      INTEGER :: lev,zon,v,v_start,v_end,av,sv,sd
      AED_REAL :: scale
      AED_REAL :: localrainl, localshade, localdrag
      LOGICAL :: splitZone
      TYPE(aed_variable_t),POINTER :: tvar
   !
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers, and output flux
         !# needs disaggregating from each zone back to the overlying cells/layers

!$OMP DO
         DO zon=1,aed_n_zones
            CALL define_zone_column(column_sed, zon, aed_n_zones, 1)

            !# Reinitialise flux_ben to be repopulated for this zone
            flux_ben = zero_
            flux_pel_pre = zero_

            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_env(1)%z_sed_zones ! CAB???

            sv = 0 ; sd = 0

            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%diag ) THEN  !# Diagnostic variable
                  IF ( tvar%sheet ) THEN
                     sd = sd + 1
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
                  ENDIF
               ELSEIF ( .NOT. tvar%extern ) THEN !# State variable
                  IF ( tvar%sheet ) THEN
                     sv = sv + 1
                     IF ( tvar%bot ) THEN
                        column_sed(av)%cell_sheet => aedZones(zon)%z_cc(1, n_vars+sv)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            !print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)

            IF ( benthic_mode .EQ. 3 ) THEN
               !# Zone is able to operated on by riparian and dry methods
               CALL aed_calculate_riparian(column_sed, zon, aedZones(zon)%z_env(1)%z_pc_wet) ! CAB???
               IF (aedZones(zon)%z_env(1)%z_pc_wet .EQ. 0. ) CALL aed_calculate_dry(column_sed, zon) !CAB ???

               !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
               CALL aed_rain_loss(column, 1, localrainl);
               IF (link_rain_loss) rain_factor = localrainl

               !# update feedback arrays to shade the water (ie reduce incoming light, Io)
               CALL aed_light_shading(column, 1, localshade)
               IF (link_solar_shade) sw_factor = localshade

               !# now the bgc updates are complete, update links to host model
               CALL aed_bio_drag(column, 1, localdrag)
               IF (link_bottom_drag) friction = localdrag
            ENDIF

            !# Calculate temporal derivatives due to benthic processes.
            !# They are stored in flux_ben (benthic vars) and flux_pel (water vars)
            flux_pel_pre = flux_pel

!           print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
            CALL aed_calculate_benthic_zone(column_sed, 1, zon)  !CAB???

            !# Record benthic fluxes in the zone array
            flux_zon(:, zon) = flux_ben(:)

            !# Now we have to find out the water column flux that occured and
            !# disaggregate it to relevant layers
            flux_pel_z(:, zon) = flux_pel(:, zon)-flux_pel_pre(:, zon)
         ENDDO
!$OMP END DO

         !# Disaggregation of zone induced fluxes to overlying layers
         v_start = 1 ; v_end = n_vars
         zon = aed_n_zones
         DO lev=wlev,1,-1
            IF ( zon .GT. 1 ) THEN
               IF (lev .GT. 1) THEN
                ! splitZone = lheights(lev-1) < aedZones(zon-1)%z_env(lev-1)%z_height !CAB ???
                  splitZone = data(col)%lheights(lev-1) < aedZones(zon-1)%z_env(1)%z_height !CAB ???
               ELSE
                  splitZone = 0.0 < aedZones(zon-1)%z_env(1)%z_height ! CAB???
               ENDIF
            ELSE
               splitZone = .FALSE.
            ENDIF

            IF (splitZone) THEN
               IF (lev .GT. 1) THEN
            !     scale = (aedZones(zon-1)%z_env(lev-1)%z_height - lheights(lev-1)) / (lheights(lev) - lheights(lev-1)) ! CAB ???
                  scale = (aedZones(zon-1)%z_env(1)%z_height - data(col)%lheights(lev-1)) / &
                                              (data(col)%lheights(lev) - data(col)%lheights(lev-1)) ! CAB ???
               ELSE
                  scale = (aedZones(zon-1)%z_env(1)%z_height - 0.0) / (data(col)%lheights(lev) - 0.0) ! CAB???
               ENDIF
               flux_pel(v_start:v_end,lev) = flux_pel_z(v_start:v_end,zon) * scale

               zon = zon - 1

               flux_pel(v_start:v_end,lev) = flux_pel(v_start:v_end,lev) + &
                                        flux_pel_z(v_start:v_end,zon) * (1.0 - scale)
            ELSE
               flux_pel(v_start:v_end,lev) = flux_pel_z(v_start:v_end,zon)
            ENDIF
         ENDDO
         !# Limit flux out of bottom waters to concentration of that layer
         !# i.e. don't flux out more than is there & distribute
         !# bottom flux into pelagic over bottom box (i.e., divide by layer height).
         !# scaled to proportion of layer area that is "bottom"
         DO lev=1,wlev
            if(lev>1)flux_pel(:, lev) = flux_pel(:, lev) * (data(col)%area(lev)-data(col)%area(lev-1))/data(col)%area(lev)
            DO v=v_start,v_end
              IF ( data(col)%cc(1, v) .GE. 0.0 ) flux_pel(lev, v) = &
                             max(-1.0 * data(col)%cc(v, lev), flux_pel(v, lev)/data(col)%dz(lev))
            END DO
         ENDDO
      ELSE
         !# Sediment zones are not simulated and therefore just operate on the bottom-most
         !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
         !# applied on flanks of the remaining layers, but note this is not suitable for
         !# model configurations where mass balance of benthic variables is required.

         !# Calculate temporal derivatives due to exchanges at the sediment/water interface
         IF ( zone_var .GE. 1 ) column(zone_var)%cell_sheet => aedZones(1)%z_env(1)%z_sed_zones ! CAB ???
         CALL aed_calculate_benthic(column, bot)

         !# Limit flux out of bottom layers to concentration of that layer
         !# i.e. don't flux out more than is there is. Then
         !# distribute bottom flux into pelagic over bottom box (i.e., divide by layer height)
         !# Skip -ve values, as GEO_ubalchg is -ve and doesnt not comply with this logic
         v_start = 1 ; v_end = n_vars
         DO v=v_start,v_end
            IF ( data(col)%cc(v, bot) .GE. 0.0 ) &
               flux_pel(v, bot) = max(-1.0 * data(col)%cc(v, bot), flux_pel(v, bot)/data(col)%dz(bot))
         END DO

         IF ( benthic_mode .EQ. 1 ) THEN
!$OMP DO
            DO lev=2,wlev
               !# Calculate temporal derivatives due to benthic fluxes.
               CALL aed_calculate_benthic(column, lev)

               !# Limit flux out of bottom layers to concentration of that layer
               !# i.e. don't flux out more than is there
               !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
               !# scaled to proportion of area that is "bottom"
               DO v=v_start,v_end
                  IF ( data(col)%cc(v, lev) .GE. 0.0 ) flux_pel(v, lev) = &
                                   max(-1.0 * data(col)%cc(v, lev), flux_pel(v, lev)/data(col)%dz(lev))
               END DO
               flux_pel(:, lev) = flux_pel(:, lev) * (data(col)%area(lev)-data(col)%area(lev-1))/data(col)%area(lev)
            ENDDO
!$OMP END DO
         ENDIF
      ENDIF
   END SUBROUTINE glm_benthics
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE calculate_fluxes(column, wlev)
   !----------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: wlev
   !
   !LOCALS
      INTEGER :: lev
      INTEGER :: layer_map(wlev)
      LOGICAL :: glm_version = .FALSE.
   !
   !----------------------------------------------------------------------------
   !BEGIN
      flux_pel = zero_
      flux_atm = zero_
      flux_ben = zero_
      flux_zon = zero_

      flux_pel_pre = zero_
      flux_pel_z = zero_

      !# Start with updating column diagnostics (currently only used for light)

      !# WATER COLUMN UPDATING
      DO lev=1,wlev
         IF (bottom_one) THEN   !# not TFV_VARIANT
            layer_map(lev) = 1 + wlev-lev
         ELSE
            layer_map(lev) = lev
         ENDIF
      ENDDO
      CALL aed_calculate_column(column, layer_map)

      !# Now do the general calculation all flux terms for rhs in mass/m3/s
      !# Includes (i) benthic flux, (ii) surface exchange and (ii) kinetic updates in each cell
      !# as calculated by glm

      !# BENTHIC FLUXES
      IF ( glm_version ) THEN
         CALL glm_benthics(column, wlev)
      ELSE
! CAB : not yet
!        IF ( do_zone_averaging ) THEN
!           flux_pel(:,wlev) = flux_pel(:,wlev) + flux_pel_z(bot, z) !/h(wlev)
!
!           !# Calculate temporal derivatives due to benthic exchange processes.
!           CALL aed_calculate_benthic(column, bot, .FALSE.)
!        ELSE
            CALL aed_calculate_benthic(column, bot)
!        ENDIF

         !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
         flux_pel(:,bot) = flux_pel(:,bot)/data(col)%lheights(top)
      ENDIF

      !# SURFACE FLUXES
      !# Calculate temporal derivatives due to air-water exchange.
      IF (doSurface) THEN !# no surface exchange under ice cover
         CALL aed_calculate_surface(column, top)

         !# Distribute the fluxes into pelagic surface layer
         flux_pel(:, top) = flux_pel(:, top) + flux_atm(:)/data(col)%dz(top)
      ENDIF

      !# WATER CELL KINETICS
      !# Add pelagic sink and source terms in cells of all depth levels.
      DO lev=1,wlev
         CALL aed_calculate(column, lev)
      ENDDO
   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END SUBROUTINE aed_run_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, col, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: col, nlev
!
!LOCALS
   INTEGER :: i, start, end, step, first
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_; localext_up = zero_

   IF (bottom_one) THEN
      first = nlev ; start = nlev-1 ; end = 1    ; step = -1
   ELSE
      first = 1    ; start = 2      ; end = nlev ; step = 1
   ENDIF

   ! Surface Kd
   CALL aed_light_extinction(column, first, localext)

   ! Surface PAR
   data(col)%par(first) = par_fraction * data(col)%rad(first) * EXP( -(Kw+localext)*1e-6*data(col)%dz(first) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = start,end,step
      localext_up = localext
      CALL aed_light_extinction(column, i, localext)

      data(col)%par(i) = data(col)%par(i+1) * EXP( -(Kw + localext_up) * data(col)%dz(i+1) )

      IF (bioshade_feedback) data(col)%extc(i) = Kw + localext
   ENDDO
END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_clean_model()
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed_delete()
   ! Deallocate internal arrays
   IF (ALLOCATED(ws))   DEALLOCATE(ws)
   IF (ALLOCATED(max_)) DEALLOCATE(max_)
   IF (ALLOCATED(min_)) DEALLOCATE(min_)
END SUBROUTINE aed_clean_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed_var_index(name)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*) :: name
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv
   INTEGER i,v
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0
   DO i=1, n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
            v = v + 1
            IF ( name .EQ. tv%name ) THEN
               aed_var_index = v
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   aed_var_index = -1
END FUNCTION aed_var_index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_api
