!###############################################################################
!#                                                                             #
!# aed_api.F90                                                                 #
!#                                                                             #
!# A generic interface between model and libaed-xxx                            #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2024 - The University of Western Australia                        #
!#                                                                             #
!#  This file is part of GLM (General Lake Model)                              #
!#                                                                             #
!#  GLM is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  GLM is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!###############################################################################

#include "aed_api.h"

#ifndef __GFORTRAN__
#  ifndef isnan
#    define isnan(x) ieee_is_nan(x)
#  endif
#endif


#define GLM_VARIANT   1
#define TFV_VARIANT   0

!-------------------------------------------------------------------------------
MODULE aed_api
!
   USE ISO_C_BINDING
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private
!
   PUBLIC aed_init_model,       &
          aed_run_model,        &
          aed_var_index,        &
          aed_clean_model,      &
          api_config_t,         &
          aed_config_model,     &
          api_env_t,            &
          aed_set_model_env,    &
          api_data_t,           &
          aed_set_model_data,   &
          aed_mobility_t,       &
          aed_set_mobility

   !#===========================================================#!
   TYPE api_config_t
      INTEGER  :: MaxLayers

      LOGICAL  :: mobility_off
      LOGICAL  :: bioshade_feedback
      LOGICAL  :: repair_state
      LOGICAL  :: link_rain_loss
      LOGICAL  :: link_solar_shade
      LOGICAL  :: link_bottom_drag
      LOGICAL  :: ice

      INTEGER  :: split_factor
      INTEGER  :: benthic_mode

      AED_REAL :: rain_factor
      AED_REAL :: sw_factor
      AED_REAL :: friction

      AED_REAL :: Kw

      AED_REAL :: nir_fraction =  0.52   ! 0.51
      AED_REAL :: par_fraction =  0.43   ! 0.45
      AED_REAL :: uva_fraction =  0.048  ! 0.035
      AED_REAL :: uvb_fraction =  0.002  ! 0.005
   END TYPE api_config_t
   !#===========================================================#!

   !#===========================================================#!
   TYPE api_data_t
      AED_REAL,DIMENSION(:,:),POINTER :: cc         => null()
      AED_REAL,DIMENSION(:),  POINTER :: cc_hz      => null()
      AED_REAL,DIMENSION(:,:),POINTER :: cc_diag    => null()
      AED_REAL,DIMENSION(:),  POINTER :: cc_diag_hz => null()
   END TYPE api_data_t
   !#===========================================================#!

   !#===========================================================#!
   TYPE api_env_t
      AED_REAL,DIMENSION(:),POINTER :: height    => null() !# layer height (previously "h")
      AED_REAL,DIMENSION(:),POINTER :: area      => null() !# layer area
      AED_REAL,DIMENSION(:),POINTER :: dz        => null() !# layer thickness
      AED_REAL,DIMENSION(:),POINTER :: depth     => null() !# layer_depth (previously "z")
      AED_REAL,DIMENSION(:),POINTER :: col_depth => null() !# col_depth (sheet - total depth of column)

      AED_REAL,DIMENSION(:),POINTER :: temp      => null() !# temperature
      AED_REAL,DIMENSION(:),POINTER :: salt      => null() !# salinity
      AED_REAL,DIMENSION(:),POINTER :: rho       => null() !# density
      AED_REAL,DIMENSION(:),POINTER :: rad       => null()
      AED_REAL,DIMENSION(:),POINTER :: extc      => null() !# extinction coefficient
      AED_REAL,DIMENSION(:),POINTER :: tss       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss1       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss2       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss3       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss4       => null()
      AED_REAL,DIMENSION(:),POINTER :: cvel      => null() !# cell velocity
      AED_REAL,DIMENSION(:),POINTER :: vvel      => null() !# vertical velocity
      AED_REAL,DIMENSION(:),POINTER :: bio_drag  => null()
      AED_REAL,DIMENSION(:),POINTER :: I_0       => null() !# par_sf
      AED_REAL,DIMENSION(:),POINTER :: wnd       => null()
      AED_REAL,DIMENSION(:),POINTER :: air_temp  => null()
      AED_REAL,DIMENSION(:),POINTER :: air_pres  => null()
      AED_REAL,DIMENSION(:),POINTER :: rain      => null()
      AED_REAL,DIMENSION(:),POINTER :: evap      => null()
      AED_REAL,DIMENSION(:),POINTER :: humidity  => null()
      AED_REAL,DIMENSION(:),POINTER :: longwave  => null()
      AED_REAL,DIMENSION(:),POINTER :: bathy     => null()
      AED_REAL,DIMENSION(:),POINTER :: rainloss  => null()
      AED_REAL,DIMENSION(:),POINTER :: ustar_bed => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_uorb   => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_t      => null()
      AED_REAL,DIMENSION(:),POINTER :: layer_stress => null()
      AED_REAL,DIMENSION(:),POINTER :: sed_zones => null()
      AED_REAL,DIMENSION(:),POINTER :: pres      => null()

      AED_REAL,DIMENSION(:),POINTER :: nir => null()
      AED_REAL,DIMENSION(:),POINTER :: par => null()
      AED_REAL,DIMENSION(:),POINTER :: uva => null()
      AED_REAL,DIMENSION(:),POINTER :: uvb => null()

      INTEGER, DIMENSION(:,:),POINTER :: mat_id => null()
      LOGICAL, DIMENSION(:),  POINTER :: active => null()

      AED_REAL,POINTER :: longitude => null()
      AED_REAL,POINTER :: latitude  => null()

      AED_REAL,POINTER :: yearday  => null()
      AED_REAL,POINTER :: timestep => null()
   END TYPE api_env_t
   !#===========================================================#!

!  !#===========================================================#!
!  TYPE api_water_col_t
!     TYPE(aed_column_t),POINTER :: column(:) => null()
!     TYPE(api_env_t),POINTER :: env          => null()
!  END TYPE api_water_col_t
!  !#===========================================================#!
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   CHARACTER(len=80) :: cfg_fname = "none"

   AED_REAL :: Kw !, Ksed

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
!  LOGICAL :: do_zone_averaging = .FALSE.
   LOGICAL :: link_solar_shade = .TRUE.
   LOGICAL :: link_rain_loss = .FALSE.
   LOGICAL :: link_bottom_drag = .FALSE.
!  LOGICAL :: link_surface_drag = .FALSE.
!  LOGICAL :: link_water_density = .FALSE.
!  LOGICAL :: link_water_clarity = .FALSE.
   LOGICAL :: link_ext_par = .FALSE.

   !-------------------------------------------------------------
   !# External variables
   !-------------------------------------------------------------

   !# Main arrays storing/pointing to the state and diagnostic variables
   AED_REAL,DIMENSION(:,:),POINTER :: cc          => null()
   AED_REAL,DIMENSION(:),  POINTER :: cc_hz       => null()
   AED_REAL,DIMENSION(:,:),POINTER :: cc_diag     => null()
   AED_REAL,DIMENSION(:),  POINTER :: cc_diag_hz  => null()

   AED_REAL,DIMENSION(:),POINTER :: temp      => null()
   AED_REAL,DIMENSION(:),POINTER :: salt      => null()
   AED_REAL,DIMENSION(:),POINTER :: rho       => null()
   AED_REAL,DIMENSION(:),POINTER :: dz        => null()
   AED_REAL,DIMENSION(:),POINTER :: lheights  => null()
   AED_REAL,DIMENSION(:),POINTER :: area      => null()
   AED_REAL,DIMENSION(:),POINTER :: col_depth => null()
   AED_REAL,DIMENSION(:),POINTER :: depth     => null()
   AED_REAL,DIMENSION(:),POINTER :: extc      => null()
   AED_REAL,DIMENSION(:),POINTER :: tss       => null()
   AED_REAL,DIMENSION(:),POINTER :: ss1       => null()
   AED_REAL,DIMENSION(:),POINTER :: ss2       => null()
   AED_REAL,DIMENSION(:),POINTER :: ss3       => null()
   AED_REAL,DIMENSION(:),POINTER :: ss4       => null()
   AED_REAL,DIMENSION(:),POINTER :: cvel      => null()  !# cell velocity
   AED_REAL,DIMENSION(:),POINTER :: bio_drag  => null()
   AED_REAL,DIMENSION(:),POINTER :: rad       => null()
   AED_REAL,DIMENSION(:),POINTER :: I_0       => null()
   AED_REAL,DIMENSION(:),POINTER :: wnd       => null()
   AED_REAL,DIMENSION(:),POINTER :: air_temp  => null()
   AED_REAL,DIMENSION(:),POINTER :: air_pres  => null()
   AED_REAL,DIMENSION(:),POINTER :: rain      => null()
   AED_REAL,DIMENSION(:),POINTER :: evap      => null()
   AED_REAL,DIMENSION(:),POINTER :: humidity  => null()
   AED_REAL,DIMENSION(:),POINTER :: longwave  => null()
   AED_REAL,DIMENSION(:),POINTER :: bathy     => null()
   AED_REAL,DIMENSION(:),POINTER :: rainloss  => null()
   AED_REAL,DIMENSION(:),POINTER :: ustar_bed => null()
   AED_REAL,DIMENSION(:),POINTER :: wv_uorb   => null()
   AED_REAL,DIMENSION(:),POINTER :: wv_t      => null()
   AED_REAL,DIMENSION(:),POINTER :: vvel      => null()  !# vertical velocity
   AED_REAL,DIMENSION(:),POINTER :: layer_stress => null()
   AED_REAL,DIMENSION(:),POINTER :: sed_zones    => null()
   AED_REAL,DIMENSION(:),POINTER :: pres         => null()

   !# Arrays for environmental variables (used if they are not supplied externally)
   AED_REAL,DIMENSION(:),POINTER :: nir => null()
   AED_REAL,DIMENSION(:),POINTER :: par => null()
   AED_REAL,DIMENSION(:),POINTER :: uva => null()
   AED_REAL,DIMENSION(:),POINTER :: uvb => null()

   INTEGER, DIMENSION(:,:),POINTER :: mat_id => null()
   INTEGER, DIMENSION(:,:),POINTER :: colnums => null()
   LOGICAL, DIMENSION(:),  POINTER :: active => null()

!  !# Maps to nearest cell with water (for riparian exchange)
! may be tuflow specific
!  AED_REAL,DIMENSION(:),POINTER :: nearest_active => null()
!  AED_REAL,DIMENSION(:),POINTER :: nearest_depth => null()
!  INTEGER, DIMENSION(:),POINTER :: route_table => null()

   !# Maps of surface, bottom and wet/dry (active) cells
!  INTEGER,DIMENSION(:),POINTER :: surf_map => null()
!  INTEGER,DIMENSION(:),POINTER :: benth_map => null()
!  LOGICAL,DIMENSION(:),POINTER :: active => null()

   !# Arrays for work, vertical movement (ws), and cross-boundary fluxes
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: ws
   AED_REAL,DIMENSION(:),ALLOCATABLE   :: min_, max_

   !# To support light
   AED_REAL,POINTER :: yearday => null()
   AED_REAL,POINTER :: timestep => null()
   AED_REAL,POINTER :: longitude => null()
   AED_REAL,POINTER :: latitude => null()

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

  !-----------------------------------------------------------------------------
  INTERFACE

    SUBROUTINE aed_mobility_t(N,dt,h,A,ww,min_C,cc)
       INTEGER,INTENT(in)     :: N       !# number of vertical layers
       AED_REAL,INTENT(in)    :: dt      !# time step (s)
       AED_REAL,INTENT(in)    :: h(*)    !# layer thickness (m)
       AED_REAL,INTENT(in)    :: A(*)    !# layer areas (m2)
       AED_REAL,INTENT(in)    :: ww(*)   !# vertical speed (m/s)
       AED_REAL,INTENT(in)    :: min_C   !# minimum allowed cell concentration
       AED_REAL,INTENT(inout) :: cc(*)   !# cell concentration
    END SUBROUTINE aed_mobility_t

  END INTERFACE
  !-----------------------------------------------------------------------------

  PROCEDURE(aed_mobility_t),POINTER :: doMobility => null()

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_show_vars
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
INTEGER FUNCTION aed_init_model(fname, NumWQ_Vars, NumWQ_Ben, NumWQ_Diag, NumWQ_DiagS)
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
   print *,'    aed_api built using intel fortran version ', __INTEL_COMPILER
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

   print *,'    libaed enabled.... aed_init_model processing: ', TRIM(fname)
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
   print "(  5X,'AED : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
   print "(  5X,'AED : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

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

   aed_init_model = n_aed_vars
END FUNCTION aed_init_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_mobility(mobl)
!-------------------------------------------------------------------------------
!ARGUMENTS
   PROCEDURE(aed_mobility_t),POINTER :: mobl
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
!
    doMobility => mobl
END SUBROUTINE aed_set_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_config_model(conf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_config_t), INTENT(in) :: conf
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
   benthic_mode = conf%benthic_mode

   rain_factor = conf%rain_factor
   sw_factor = conf%sw_factor
   friction = conf%friction

   Kw = conf%Kw
END SUBROUTINE aed_config_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_data(dat)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_data_t), INTENT(in) :: dat
!
!LOCALS
   INTEGER :: av, v, sv, status
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   cc         => dat%cc
!# Eventually this will be properly separated
!  cc_hz      => dat%cc_hz
   cc_hz      => dat%cc(1,n_vars:n_vars+n_vars_ben)
   cc_diag    => dat%cc_diag
   cc_diag_hz => dat%cc_diag_hz

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
            cc(:, n_vars+sv) = tvar%initial
         ELSE
            v = v + 1
            cc(:, v) = tvar%initial
         ENDIF
      ENDIF
   ENDDO

   !# Allocate array with vertical movement rates (m/s, positive for upwards),
   !# and set these to the values provided by the model.
   !# allocated for all vars even though only state vars entries will be used
   ALLOCATE(ws(MaxLayers, n_aed_vars),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = zero_

   !# Trigger an error if we don't have all we need.
   CALL check_data
END SUBROUTINE aed_set_model_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_env(env)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_env_t),INTENT(in) :: env
!
!LOCALS
   INTEGER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   yearday  => env%yearday
   timestep => env%timestep

   dt_eff = timestep/FLOAT(split_factor)

   longitude => env%longitude
   latitude  => env%latitude

   temp         => env%temp
   salt         => env%salt
   rho          => env%rho
   dz           => env%dz
   lheights     => env%height
   area         => env%area
   depth        => env%depth
   col_depth    => env%col_depth
   extc         => env%extc
   tss          => env%tss
   ss1          => env%ss1
   ss2          => env%ss2
   ss3          => env%ss3
   ss4          => env%ss4
   cvel         => env%cvel
   vvel         => env%vvel
   bio_drag     => env%bio_drag
   rad          => env%rad
   I_0          => env%I_0
   wnd          => env%wnd
   air_temp     => env%air_temp
   air_pres     => env%air_pres
   rain         => env%rain
   evap         => env%evap
   humidity     => env%humidity
   longwave     => env%longwave
   bathy        => env%bathy
   rainloss     => env%rainloss
   ustar_bed    => env%ustar_bed
   wv_uorb      => env%wv_uorb
   wv_t         => env%wv_t
   layer_stress => env%layer_stress
   sed_zones    => env%sed_zones
   pres         => env%pres

   nir => env%nir
   par => env%par
   uva => env%uva
   uvb => env%uvb

   mat_id => env%mat_id
   active => env%active

   IF (ASSOCIATED(temp))           tv=aed_provide_global('temperature','temperature',           'celsius')
   IF (ASSOCIATED(salt))           tv=aed_provide_global('salinity',   'salinity',              'g/kg'   )
   IF (ASSOCIATED(rho))            tv=aed_provide_global('density',    'density',               'kg/m3'  )
   IF (ASSOCIATED(dz))             tv=aed_provide_global('layer_ht',   'layer heights',         'm'      )
   IF (ASSOCIATED(area))           tv=aed_provide_global('layer_area', 'layer area',            'm2'     )
   IF (ASSOCIATED(depth))          tv=aed_provide_global('depth',      'depth',                 'm'      )
   IF (ASSOCIATED(extc))           tv=aed_provide_global('extc_coef',  'extinction coefficient','/m'     )
   IF (ASSOCIATED(tss))            tv=aed_provide_global('tss',        'tss',                   'g/m3'   )
   IF (ASSOCIATED(ss1))            tv=aed_provide_global('ss1',        'ss1',                   'g/m3'   )
   IF (ASSOCIATED(ss2))            tv=aed_provide_global('ss2',        'ss2',                   'g/m3'   )
   IF (ASSOCIATED(ss3))            tv=aed_provide_global('ss3',        'ss3',                   'g/m3'   )
   IF (ASSOCIATED(ss4))            tv=aed_provide_global('ss4',        'ss4',                   'g/m3'   )
   IF (ASSOCIATED(cvel))           tv=aed_provide_global('cell_vel',   'cell velocity',         'm/s'    )
   IF (ASSOCIATED(pres))           tv=aed_provide_global('pressure',   'pressure',              ''       )

   tv=aed_provide_global('nir',        'nir',                   'W/m2'   )
   tv=aed_provide_global('par',        'par',                   'W/m2'   )
   tv=aed_provide_global('uva',        'uva',                   'W/m2'   )
   tv=aed_provide_global('uvb',        'uvb',                   'W/m2'   )

   IF (ASSOCIATED(I_0))            tv=aed_provide_sheet_global('par_sf',        'par_sf',            'W/m2'          )
   IF (ASSOCIATED(wnd))            tv=aed_provide_sheet_global('wind_speed',    'wind speed',        'm/s'           )
   IF (ASSOCIATED(air_temp))       tv=aed_provide_sheet_global('air_temp',      'air temperature',   'celsius'       )
   IF (ASSOCIATED(air_pres))       tv=aed_provide_sheet_global('air_pres',      'air pressure',      'Pa'            )
   IF (ASSOCIATED(rain))           tv=aed_provide_sheet_global('rain',          'rainfall',          'm/s'           )
   IF (ASSOCIATED(evap))           tv=aed_provide_sheet_global('evap',          'evaporation',       'm/s'           )
   IF (ASSOCIATED(humidity))       tv=aed_provide_sheet_global('humidity',      'relative humidity', '-'             )
   IF (ASSOCIATED(longwave))       tv=aed_provide_sheet_global('longwave',      'longwave',          'W/m2'          )
   IF (ASSOCIATED(bathy))          tv=aed_provide_sheet_global('bathy',         'bathy',             'm above datum' )
   IF (ASSOCIATED(rainloss))       tv=aed_provide_sheet_global('rainloss',      'rain loss',         'm/s'           )
   IF (ASSOCIATED(layer_stress))   tv=aed_provide_sheet_global('taub',          'layer stress',      'N/m2'          )
   IF (ASSOCIATED(sed_zones))      tv=aed_provide_sheet_global('sed_zone',      'sediment zone',     '-'             )

 ! IF (ASSOCIATED(col_num))        tv=aed_provide_sheet_global('col_num',       'column number',     '-'             )
   IF (ASSOCIATED(col_depth))      tv=aed_provide_sheet_global('col_depth',     'column water depth','m above bottom')
 ! IF (ASSOCIATED(nearest_active)) tv=aed_provide_sheet_global('nearest_active','nearest active',    '-'             )
 ! IF (ASSOCIATED(nearest_depth))  tv=aed_provide_sheet_global('nearest_depth', 'nearest depth',     'm'             )

   IF (ASSOCIATED(mat_id))    tv=aed_provide_sheet_global('material',      'material',          '-'             )

   IF (ASSOCIATED(longitude)) tv=aed_provide_sheet_global('longitude',     'longitude',         'radians'       )
   IF (ASSOCIATED(latitude))  tv=aed_provide_sheet_global('latitude',      'latitude',          'radians'       )

   IF (ASSOCIATED(yearday))   tv=aed_provide_sheet_global('yearday', 'yearday', 'day'    )
   IF (ASSOCIATED(timestep))  tv=aed_provide_sheet_global('timestep','timestep','seconds')
END SUBROUTINE aed_set_model_env
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_data
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
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = ASSOCIATED(temp)
            CASE ( 'salinity' )    ; tvar%found = ASSOCIATED(salt)
            CASE ( 'density' )     ; tvar%found = ASSOCIATED(rho)
            CASE ( 'layer_ht' )    ; tvar%found = ASSOCIATED(dz)
            CASE ( 'extc_coef' )   ; tvar%found = ASSOCIATED(extc)
            CASE ( 'tss' )         ; tvar%found = ASSOCIATED(tss)
            CASE ( 'cell_vel' )    ; tvar%found = ASSOCIATED(cvel)
            CASE ( 'par' )         ; tvar%found = ASSOCIATED(par)
            CASE ( 'nir' )         ; tvar%found = ASSOCIATED(nir)
            CASE ( 'uva' )         ; tvar%found = ASSOCIATED(uva)
            CASE ( 'uvb' )         ; tvar%found = ASSOCIATED(uvb)
            CASE ( 'pressure' )    ; tvar%found = ASSOCIATED(pres)
            CASE ( 'depth' )       ; tvar%found = ASSOCIATED(depth)
            CASE ( 'sed_zone' )    ; tvar%found = ASSOCIATED(sed_zones)
            CASE ( 'wind_speed' )  ; tvar%found = ASSOCIATED(wnd)
            CASE ( 'par_sf' )      ; tvar%found = ASSOCIATED(I_0)
            CASE ( 'taub' )        ; tvar%found = ASSOCIATED(layer_stress)
            CASE ( 'col_depth' )   ; tvar%found = ASSOCIATED(col_depth)
            CASE ( 'layer_area' )  ; tvar%found = ASSOCIATED(area)
            CASE ( 'rain' )        ; tvar%found = ASSOCIATED(rain)
            CASE ( 'evap' )        ; tvar%found = ASSOCIATED(evap)
            CASE ( 'air_temp' )    ; tvar%found = ASSOCIATED(air_temp)
            CASE ( 'air_pres' )    ; tvar%found = ASSOCIATED(air_pres)
            CASE ( 'humidity' )    ; tvar%found = ASSOCIATED(humidity)
            CASE ( 'longitude' )   ; tvar%found = ASSOCIATED(longitude)
            CASE ( 'latitude' )    ; tvar%found = ASSOCIATED(latitude)
            CASE ( 'yearday' )     ; tvar%found = ASSOCIATED(yearday)
            CASE ( 'timestep' )    ; tvar%found = ASSOCIATED(timestep)

            CASE ( 'rainloss' )    ; tvar%found = ASSOCIATED(rainloss)
            CASE ( 'material' )    ; tvar%found = ASSOCIATED(mat_id)
            CASE ( 'bathy' )       ; tvar%found = ASSOCIATED(bathy)
            CASE ( 'ss1' )         ; tvar%found = ASSOCIATED(ss1)
            CASE ( 'ss2' )         ; tvar%found = ASSOCIATED(ss2)
            CASE ( 'ss3' )         ; tvar%found = ASSOCIATED(ss3)
            CASE ( 'ss4' )         ; tvar%found = ASSOCIATED(ss4)
            CASE ( 'col_num' )     ; tvar%found = ASSOCIATED(colnums)
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
END SUBROUTINE check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if GLM_VARIANT
!###############################################################################
SUBROUTINE define_column(column, top, flux_pel, flux_atm, flux_ben)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: top
   AED_REAL,TARGET,INTENT(inout) :: flux_pel(:,:) !# (n_layers, n_vars)
   AED_REAL,TARGET,INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL,TARGET,INTENT(inout) :: flux_ben(:)   !# (n_vars)
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
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(:)
            CASE ( 'salinity' )    ; column(av)%cell => salt(:)
            CASE ( 'density' )     ; column(av)%cell => rho(:)
            CASE ( 'layer_ht' )    ; column(av)%cell => dz(:)
            CASE ( 'extc_coef' )   ; column(av)%cell => extc(:)
            CASE ( 'tss' )         ; column(av)%cell => tss(:)
            CASE ( 'cell_vel' )    ; column(av)%cell => cvel(:)
            CASE ( 'par' )         ; column(av)%cell => par(:)
            CASE ( 'nir' )         ; column(av)%cell => nir(:)
            CASE ( 'uva' )         ; column(av)%cell => uva(:)
            CASE ( 'uvb' )         ; column(av)%cell => uvb(:)
            CASE ( 'pressure' )    ; column(av)%cell => pres(:)
            CASE ( 'depth' )       ; column(av)%cell => depth(:)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => sed_zones(1)
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd(1)
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0(1)
            CASE ( 'taub' )        ; column(av)%cell_sheet => layer_stress(1) ! CAB?
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => col_depth(1)
            CASE ( 'layer_area' )  ; column(av)%cell => area(:)
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(1)
            CASE ( 'evap' )        ; column(av)%cell_sheet => evap(1)
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp(1)
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => air_pres(1)
            CASE ( 'humidity' )    ; column(av)%cell_sheet => humidity(1)
            CASE ( 'longitude' )   ; column(av)%cell_sheet => longitude
            CASE ( 'latitude' )    ; column(av)%cell_sheet => latitude
            CASE ( 'yearday' )     ; column(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column(av)%cell_sheet => timestep
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column(av)%cell_sheet => cc_diag_hz(sd)
         ELSE
            d = d + 1
            column(av)%cell => cc_diag(:, d)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => cc(1, n_vars+sv)
!              print *,'av',av,sv
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => cc(top, n_vars+sv)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => cc(:, v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(:, v)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


#if TFV_VARIANT
!###############################################################################
SUBROUTINE define_column(column, col, cc, cc_diag, flux_pel, flux_atm, flux_ben, flux_rip)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: col
   AED_REAL, TARGET, INTENT(in) :: cc(:,:)       !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(in) :: cc_diag(:,:)  !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_rip(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av, i, top, bot
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   top = surf_map(col)
   bot = benth_map(col)

   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed_vars

      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "ERROR getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(top:bot)
            CASE ( 'salinity' )    ; column(av)%cell => salt(top:bot)
            CASE ( 'density' )     ; column(av)%cell => rho(top:bot)
            CASE ( 'layer_ht' )    ; column(av)%cell => h(top:bot)
            CASE ( 'layer_area' )  ; column(av)%cell_sheet => area(col)
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(col)
            CASE ( 'rainloss' )    ; column(av)%cell_sheet => rainloss(col)
            CASE ( 'material' )    ; IF ( do_zone_averaging ) THEN
                                        column(av)%cell_sheet => zone(zm(col))
                                     ELSE
                                        column(av)%cell_sheet => mat(col)
                                     ENDIF
            CASE ( 'bathy' )       ; column(av)%cell_sheet => bathy(col)
            CASE ( 'extc_coef' )   ; column(av)%cell => extcoeff(top:bot)
            CASE ( 'tss' )         ; column(av)%cell => tss(top:bot)
            CASE ( 'ss1' )         ; column(av)%cell => tss(top:bot)   !   For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss2' )         ; column(av)%cell => tss(top:bot)   !   For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss3' )         ; column(av)%cell => tss(top:bot)   !   For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss4' )         ; column(av)%cell => tss(top:bot)   !   For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'cell_vel' )    ; column(av)%cell => cvel(top:bot)
            CASE ( 'nir' )         ; column(av)%cell => nir(top:bot)
            CASE ( 'par' )         ; IF (link_ext_par) THEN
                                        column(av)%cell => lpar(top:bot)
                                     ELSE
                                        column(av)%cell => par(top:bot)
                                     ENDIF
            CASE ( 'uva' )         ; column(av)%cell => uva(top:bot)
            CASE ( 'uvb' )         ; column(av)%cell => uvb(top:bot)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => zone(zm(col))
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd(col)
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0(col)
            CASE ( 'taub' )        ; column(av)%cell_sheet => col_taub
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp(col)
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => air_pres(col)
            CASE ( 'humidity' )    ; column(av)%cell_sheet => humidity(col)
            CASE ( 'longwave' )    ; column(av)%cell_sheet => longwave(col)
            CASE ( 'col_num' )     ; column(av)%cell_sheet => colnums(col)
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => z(col)

            CASE ( 'nearest_active' ) ; column(av)%cell_sheet => nearest_active(col)
            CASE ( 'nearest_depth' )  ; column(av)%cell_sheet => nearest_depth(col)

            CASE ( 'longitude' )   ; column(av)%cell_sheet => lon(col)
            CASE ( 'latitude' )    ; column(av)%cell_sheet => lat(col)
            CASE ( 'yearday' )     ; column(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column(av)%cell_sheet => dt

            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         d = d + 1
         IF ( tvar%sheet ) THEN
            column(av)%cell_sheet => cc_diag(d, bot)
         ELSE
            column(av)%cell => cc_diag(d,top:bot)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, bot)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, top)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
            column(av)%flux_rip => flux_rip(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => cc(v,top:bot)
            column(av)%flux_pel => flux_pel(v,top:bot)
            column(av)%flux_ben => flux_ben(v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_rip => flux_rip(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


!###############################################################################
SUBROUTINE check_states(column, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: wlev
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
                  IF ( isnan(cc(lev, v)) ) last_naned = i
#endif
                  IF ( .NOT. isnan(min_(v)) ) THEN
                     IF ( cc(lev, v) < min_(v) ) cc(lev, v) = min_(v)
                  ENDIF
                  IF ( .NOT. isnan(max_(v)) ) THEN
                     IF ( cc(lev, v) > max_(v) ) cc(lev, v) = max_(v)
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

   TYPE (aed_column_t),ALLOCATABLE :: column(:)
   TYPE (aed_column_t),ALLOCATABLE :: column_sed(:)
   AED_REAL,ALLOCATABLE,TARGET :: flux_ben(:)
   AED_REAL,ALLOCATABLE,TARGET :: flux_atm(:)
   AED_REAL,ALLOCATABLE,TARGET :: flux_pel(:,:)
   AED_REAL,ALLOCATABLE,TARGET :: flux_zon(:,:)
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(flux_pel(MAX(wlev,aed_n_zones),n_vars+n_vars_ben))
   ALLOCATE(flux_zon(aed_n_zones, n_vars+n_vars_ben))
   ALLOCATE(flux_ben(n_vars+n_vars_ben))
   ALLOCATE(flux_atm(n_vars+n_vars_ben))
   ALLOCATE(column(n_aed_vars))
   ALLOCATE(column_sed(n_aed_vars))

   cc_diag = 0.
   cc_diag_hz = 0.

   IF ( .NOT. reinited ) CALL re_initialize()

   DO col=1, nCols
      CALL define_column(column, wlev, flux_pel, flux_atm, flux_ben)
      CALL pre_kinetics(col)
   ENDDO

   cc_diag = 0.
   cc_diag_hz = 0.

   DO col=1, nCols
      top = wlev ; bot = 1
      CALL define_column(column, wlev, flux_pel, flux_atm, flux_ben)
      CALL aed_run_column(wlev, doSurface)
   ENDDO

   DEALLOCATE(column_sed)
   DEALLOCATE(column)
   DEALLOCATE(flux_pel)
   DEALLOCATE(flux_zon)
   DEALLOCATE(flux_ben)
   DEALLOCATE(flux_atm)

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
      INTEGER  :: v, lev, zon=1
      AED_REAL :: pa = 0.
   !
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( benthic_mode .GT. 1 ) THEN
         zon = 1
         DO lev=1,wlev
!           print *,'zone_heights',height(lev),zone_height(zon),aedZones(1)%zheight,aedZones(2)%zheight
            IF (lheights(lev) .GT. aedZones(zon)%z_env(1)%z_height) THEN  !CAB ???
               sed_zones(lev) = zon * ( area(lev) - pa ) / area(lev)
               pa = area(lev)
               zon = zon+1
            ELSE
               sed_zones(lev) = zon
            ENDIF
         ENDDO
      ENDIF

!   DO split=1,split_factor
      IF (benthic_mode .GT. 1) THEN
         CALL p_calc_zone_areas(aedZones, aed_n_zones, area, lheights, wlev)
         CALL p_copy_to_zone(aedZones, aed_n_zones, lheights, cc, cc_hz, cc_diag, cc_diag_hz, wlev)
      ENDIF

      !# Time-integrate one biological time step
      CALL calculate_fluxes(wlev)

      !# Update the water column layers
      DO v = 1, n_vars
         DO lev = 1, wlev
            cc(lev, v) = cc(lev, v) + dt_eff*flux_pel(lev, v)
         ENDDO
      ENDDO

      !# Now update benthic variables, depending on whether zones are simulated
      IF ( benthic_mode .GT. 1 ) THEN
         ! Loop through benthic state variables to update their mass
         DO v = n_vars+1, n_vars+n_vars_ben
            ! Loop through each sediment zone
            DO zon = 1, aed_n_zones
               ! Update the main cc_sed data array with the
               aedZones(zon)%z_cc(1, v) = aedZones(zon)%z_cc(1, v) + dt_eff*flux_zon(zon, v)
            ENDDO
         ENDDO

         !# Distribute cc-sed benthic properties back into main cc array
         CALL p_copy_from_zone(aedZones, aed_n_zones, lheights, cc, cc_hz, cc_diag, cc_diag_hz, wlev)
      ELSE
         DO v = n_vars+1, n_vars+n_vars_ben
            cc(1, v) = cc(1, v) + dt_eff*flux_ben(v)
         ENDDO
      ENDIF

      CALL check_states(column, wlev)
!   ENDDO
   END SUBROUTINE aed_run_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if GLM_VARIANT
   !############################################################################
   SUBROUTINE pre_kinetics(col)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in) :: col
   !
   !LOCALS
      TYPE(aed_variable_t),POINTER :: tv
      AED_REAL :: min_C
      INTEGER  :: i, v
   !
   !----------------------------------------------------------------------------
   !BEGIN
!     CALL define_column(column, wlev, flux_pel, flux_atm, flux_ben)

      IF ( .NOT. mobility_off ) THEN
         v = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tv) ) THEN
               IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
                  v = v + 1
                  ws(:,i) = zero_
                  ! only for state_vars that are not sheet
                  IF ( .NOT. isnan(tv%mobility) ) THEN
                     ! default to ws that was set during initialisation
                     ws(1:wlev,i) = tv%mobility
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         DO i = 1, wlev
            ! update ws for modules that use the mobility method
            CALL aed_mobility(column, i, ws(i,:))
         ENDDO

         !# (3) Calculate source/sink terms due to the settling or rising of
         !# state variables in the water column (note that settling into benthos
         !# is done in aed_do_benthos)
         v = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tv) ) THEN
               IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
                  v = v + 1
                  !# only for state_vars that are not sheet, and also non-zero ws
                  IF ( .NOT. isnan(tv%mobility) .AND. SUM(ABS(ws(1:wlev,i)))>zero_ ) THEN
                     min_C = tv%minimum
                     CALL doMobility(wlev, timestep, dz, area, ws(:, i), min_C, cc(:, v))
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      CALL check_states(column,wlev)

!CAB ??
      IF (benthic_mode .GT. 1) THEN
         CALL p_calc_zone_areas(aedZones, aed_n_zones, area, lheights, wlev)
         CALL p_copy_to_zone(aedZones, aed_n_zones, lheights, cc, cc_hz, cc_diag, cc_diag_hz, wlev)
      ENDIF

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer
      IF (.NOT. link_ext_par) &
         CALL update_light(column, wlev)

      !# Fudge
      nir(:) = (par(:)/par_fraction) * nir_fraction
      uva(:) = (par(:)/par_fraction) * uva_fraction
      uvb(:) = (par(:)/par_fraction) * uvb_fraction
   END SUBROUTINE pre_kinetics
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


#if TFV_VARIANT
   !############################################################################
   SUBROUTINE pre_kinetics(col)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in) :: col
   !
   !LOCALS
      TYPE(aed_variable_t),POINTER :: tv
!     AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben),    &
!                 flux_rip(n_vars+n_vars_ben)
      INTEGER  :: i, lev, top, bot, v
      TYPE (aed_column_t) :: column(n_aed_vars)
   !
   !----------------------------------------------------------------------------
   !BEGIN
      ! move to next column if dry
      IF (.NOT. active(col)) RETURN

      ! identify cell indicies within the domain
      top = surf_map(col)
      bot = benth_map(col)

      ! set column data
      CALL define_column(column, col, cc, cc_diag, flux_pel, flux_atm, flux_ben, flux_rip)

      ! compute vertical settling/mobility
      v = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               ! only for state_vars that are not sheet
               IF ( .NOT. ieee_is_nan(tv%mobility) ) THEN
                  ! default to ws that was set during initialisation
                  ws(top:bot,i) = tv%mobility
               ELSE
                  ! zero nan values
                  ws(top:bot,i) = zero_
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      DO lev = top, bot
        ! update ws for modules that use the mobility method
        CALL aed_mobility(column, lev-top+1, ws(lev,:))
      ENDDO
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
         !  IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) .AND. SUM(ABS(ws(top:bot,i)))>zero_ ) THEN
         !     CALL Settling(bot-top+1, dt, h(top:bot), ws(top:bot,i), Fsed_setl(col), column(i)%cell)
         !  ENDIF
            IF ( .NOT. isnan(tv%mobility) .AND. SUM(ABS(ws(top:bot,i)))>zero_ ) THEN
               min_C = tv%minimum
               CALL doMobility(wlev, timestep, dz, area, ws(top:bot, i), min_C, cc(:, v))
            ENDIF
         ENDIF
      ENDDO
      CALL check_states(top, bot)

      !# populate local light/extc arrays one column at a time
      IF (.NOT. link_ext_par) THEN
         CALL Light(column, bot-top+1, I_0(col), extcoeff(top:bot), par(top:bot), h(top:bot))
         ! non PAR bandwidth fractions (set assuming single light extinction)
         nir(top:bot) = (par(top:bot)/par_frac) * nir_frac
         uva(top:bot) = (par(top:bot)/par_frac) * uva_frac
         uvb(top:bot) = (par(top:bot)/par_frac) * uvb_frac
      ELSE
        ! non PAR bandwidth fractions (set assuming single light extinction)
        nir(top:bot) = (lpar(top:bot)/par_frac) * nir_frac
        uva(top:bot) = (lpar(top:bot)/par_frac) * uva_frac
        uvb(top:bot) = (lpar(top:bot)/par_frac) * uvb_frac
      ENDIF
   END SUBROUTINE pre_kinetics
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


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
         SELECT CASE (tvar%name)
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
        !   CASE ( 'sed_zone' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(bot)%z_sed_zones; zone_var = av !CAB ??? (bot)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => aedZones(zon)%z_env(1)%z_sed_zones; zone_var = av
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd(1)
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0(1)
            CASE ( 'taub' )        ; column(av)%cell_sheet => layer_stress(1) !CAB ??? (bot)
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => col_depth(1)
            CASE ( 'layer_area' )  ; column(av)%cell => aedZones(zon)%z_env(:)%z_area
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(1)
            CASE ( 'evap' )        ; column(av)%cell_sheet => evap(1)
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp(1)
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => air_pres(1)
            CASE ( 'humidity' )    ; column(av)%cell_sheet => humidity(1)
            CASE ( 'longitude' )   ; column(av)%cell_sheet => longitude
            CASE ( 'latitude' )    ; column(av)%cell_sheet => latitude
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
            column(av)%cell => aedZones(zon)%z_cc_diag(:, d)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
!CAB ??
#if 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(bot, n_vars+sv)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(top, n_vars+sv)
            ENDIF
#else
            column(av)%cell_sheet => aedZones(zon)%z_cc_hz(sv)
#endif
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => aedZones(zon)%z_cc(:, v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(:, v)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
   END SUBROUTINE define_zone_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if GLM_VARIANT
   !############################################################################
   SUBROUTINE re_initialize
   !----------------------------------------------------------------------------
   !ARGUMENTS
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
!CAB ??
#if 1
          !          IF ( tvar%bot ) THEN
          !             column_sed(av)%cell_sheet => aedZones(zon)%z_cc(1, n_vars+sv)
          !          ENDIF
#else
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc(1, n_vars+sv)
#endif
                  ENDIF
               ENDIF
            ENDDO

            CALL aed_initialize_benthic(column_sed, zon)
         ENDDO
      ENDIF
      reinited = .TRUE.
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


#if TFV_VARIANT
   !############################################################################
   SUBROUTINE re_initialize()
   !----------------------------------------------------------------------------
   !ARGUMENTS
   !
   !LOCALS
   !  INTEGER  :: i, col, lev, top, bot, count, nCols
      INTEGER  :: col, lev, top, bot, count, nCols
      AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben),    &
                  flux_rip(n_vars+n_vars_ben)
      TYPE (aed_column_t) :: column(n_aed_vars)
   !
   !----------------------------------------------------------------------------
   !BEGIN
      nCols = ubound(active, 1)

      DO col=1, nCols
         top = surf_map(col)
         bot = benth_map(col)
         count = bot-top+1
         CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben, flux_rip)
         DO lev=1, count
            CALL aed_initialize(column, lev)
         ENDDO
         IF ( .NOT. do_zone_averaging ) &
            CALL aed_initialize_benthic(column, 1)
      ENDDO

      IF ( do_zone_averaging ) &
         CALL aed_initialize_zone_benthic(nCols, active, n_aed_vars, cc_diag, benth_map)

      reinited = .TRUE.
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


#if GLM_VARIANT
   !############################################################################
   SUBROUTINE calculate_fluxes(wlev)
   !----------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER, INTENT(in) :: wlev
   !
   !LOCALS
      INTEGER :: lev,zon,v,v_start,v_end,av,sv,sd
      INTEGER,ALLOCATABLE :: layer_map(:)
      AED_REAL :: scale
      AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_pre
      AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_z
      AED_REAL :: localrainl, localshade, localdrag
      LOGICAL :: splitZone
      TYPE(aed_variable_t),POINTER :: tvar
   !
   !----------------------------------------------------------------------------
   !BEGIN
   flux_pel = zero_
   flux_atm = zero_
   flux_ben = zero_
   flux_zon = zero_

   ALLOCATE(flux_pel_pre(MAX(wlev, aed_n_zones), n_vars+n_vars_ben))
   flux_pel_pre = zero_
   ALLOCATE(flux_pel_z(MAX(wlev, aed_n_zones), n_vars+n_vars_ben))
   flux_pel_z = zero_
   !# Start with updating column diagnostics (currently only used for light)

   !# (1) WATER COLUMN UPDATING
   ALLOCATE(layer_map(wlev))
   DO lev=1,wlev
      layer_map(lev) = 1 + wlev-lev
   ENDDO
   CALL aed_calculate_column(column, layer_map)
   DEALLOCATE(layer_map)

   !# Now do the general calculation all flux terms for rhs in mass/m3/s
   !# Includes (i) benthic flux, (ii) surface exchange and (ii) kinetic updates in each cell
   !# as calculated by glm

   !# (2) BENTHIC FLUXES
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

!        print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
         CALL aed_calculate_benthic_zone(column_sed, 1, zon)  !CAB???

         !# Record benthic fluxes in the zone array
         flux_zon(zon, :) = flux_ben(:)

         !# Now we have to find out the water column flux that occured and
         !# disaggregate it to relevant layers
         flux_pel_z(zon,:) = flux_pel(zon,:)-flux_pel_pre(zon,:)
      ENDDO
!$OMP END DO

      !# Disaggregation of zone induced fluxes to overlying layers
      v_start = 1 ; v_end = n_vars
      zon = aed_n_zones
      DO lev=wlev,1,-1
         IF ( zon .GT. 1 ) THEN
            IF (lev .GT. 1) THEN
             ! splitZone = lheights(lev-1) < aedZones(zon-1)%z_env(lev-1)%z_height !CAB ???
               splitZone = lheights(lev-1) < aedZones(zon-1)%z_env(1)%z_height !CAB ???
            ELSE
               splitZone = 0.0 < aedZones(zon-1)%z_env(1)%z_height ! CAB???
            ENDIF
         ELSE
            splitZone = .FALSE.
         ENDIF

         IF (splitZone) THEN
            IF (lev .GT. 1) THEN
            !  scale = (aedZones(zon-1)%z_env(lev-1)%z_height - lheights(lev-1)) / (lheights(lev) - lheights(lev-1)) !CAB ???
               scale = (aedZones(zon-1)%z_env(1)%z_height - lheights(lev-1)) / (lheights(lev) - lheights(lev-1)) !CAB ???
            ELSE
               scale = (aedZones(zon-1)%z_env(1)%z_height - 0.0) / (lheights(lev) - 0.0) ! CAB???
            ENDIF
            flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end) * scale

            zon = zon - 1

            flux_pel(lev,v_start:v_end) = flux_pel(lev,v_start:v_end) + &
                                        flux_pel_z(zon,v_start:v_end) * (1.0 - scale)
         ELSE
            flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end)
         ENDIF
      ENDDO
      !# Limit flux out of bottom waters to concentration of that layer
      !# i.e. don't flux out more than is there & distribute
      !# bottom flux into pelagic over bottom box (i.e., divide by layer height).
      !# scaled to proportion of layer area that is "bottom"
      DO lev=1,wlev
         if(lev>1)flux_pel(lev, :) = flux_pel(lev, :) * (area(lev)-area(lev-1))/area(lev)
         DO v=v_start,v_end
           IF ( cc(1, v) .GE. 0.0 ) flux_pel(lev, v) = max(-1.0 * cc(lev, v), flux_pel(lev, v)/dz(lev))
         END DO
      ENDDO
   ELSE
      !# Sediment zones are not simulated and therefore just operate on the bottom-most
      !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
      !# applied on flanks of the remaining layers, but note this is not suitable for
      !# model configurations where mass balance of benthic variables is required.

      !# Calculate temporal derivatives due to exchanges at the sediment/water interface
      IF ( zone_var .GE. 1 ) column(zone_var)%cell_sheet => aedZones(1)%z_env(1)%z_sed_zones ! CAB ???
      CALL aed_calculate_benthic(column, 1)

      !# Limit flux out of bottom layers to concentration of that layer
      !# i.e. don't flux out more than is there is. Then
      !# distribute bottom flux into pelagic over bottom box (i.e., divide by layer height)
      !# Skip -ve values, as GEO_ubalchg is -ve and doesnt not comply with this logic
      v_start = 1 ; v_end = n_vars
      DO v=v_start,v_end
        IF ( cc(1, v) .GE. 0.0 ) flux_pel(1, v) = max(-1.0 * cc(1, v), flux_pel(1, v)/dz(1))
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
              IF ( cc(lev, v) .GE. 0.0 ) flux_pel(lev, v) = &
                                   max(-1.0 * cc(lev, v), flux_pel(lev, v)/dz(lev))
            END DO
            flux_pel(lev, :) = flux_pel(lev, :) * (area(lev)-area(lev-1))/area(lev)
         ENDDO
!$OMP END DO
      ENDIF
   ENDIF

   !# (3) SURFACE FLUXES
   !# Calculate temporal derivatives due to air-water exchange.
   IF (doSurface) THEN !# no surface exchange under ice cover
      CALL aed_calculate_surface(column, wlev)

      !# Distribute the fluxes into pelagic surface layer
      flux_pel(wlev, :) = flux_pel(wlev, :) + flux_atm(:)/dz(wlev)
   ENDIF

   !# (4) WATER CELL KINETICS
   !# Add pelagic sink and source terms in cells of all depth levels.
   DO lev=1,wlev
      CALL aed_calculate(column, lev)
   ENDDO

   DEALLOCATE(flux_pel_pre)
   DEALLOCATE(flux_pel_z)
   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


#if TFV_VARIANT
   !############################################################################
   SUBROUTINE calculate_fluxes(column, count, z, flux_pel, flux_atm, flux_ben, flux_rip, h)
   !----------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE (aed_column_t), INTENT(inout) :: column(:)
      INTEGER, INTENT(in) :: count, z
      AED_REAL, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
      AED_REAL, INTENT(inout) :: flux_atm(:)   !# (n_vars+n_ben)
      AED_REAL, INTENT(inout) :: flux_ben(:)   !# (n_vars+n_ben)
      AED_REAL, INTENT(inout) :: flux_rip(:)   !# (n_vars)
      AED_REAL, INTENT(inout) :: h(:)          !# (n_layers)
   !
      INTEGER :: layer_map(count)
   !
   !LOCALS
      INTEGER :: i
   !
   !----------------------------------------------------------------------------
   !BEGIN
      flux_pel = zero_
      flux_atm = zero_
      flux_ben = zero_
      flux_rip = zero_

      !# Start with updating column diagnostics (currently only used for light)
      DO i=1, count
   !     layer_map(i) = 1 + count-i
         layer_map(i) = i
      ENDDO
      CALL aed_calculate_column(column, layer_map)

      !# Calculate temporal derivatives due to air-water exchange.
      CALL aed_calculate_surface(column, 1)

      !# Distribute the fluxes into pelagic surface layer
      IF ( do_2d_atm_flux .OR. count > 1 ) &
         flux_pel(:,1) = flux_pel(:,1) + flux_atm(:)/h(1)

      IF ( do_zone_averaging ) THEN
         flux_pel(:,count) = flux_pel(:,count) + flux_pelz(:,z) !/h(count)

         !# Calculate temporal derivatives due to benthic exchange processes.
         CALL aed_calculate_benthic(column, count, .FALSE.)
      ELSE
         CALL aed_calculate_benthic(column, count)
      ENDIF

      !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
      flux_pel(:,count) = flux_pel(:,count)/h(count)

      !# Add pelagic sink and source terms for all depth levels.
      DO i=1,count
         CALL aed_calculate(column, i)
      ENDDO
   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


END SUBROUTINE aed_run_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: nlev
!
!LOCALS
   INTEGER :: i
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_; localext_up = zero_

   ! Surface Kd
   CALL aed_light_extinction(column, nlev, localext)

   ! Surface PAR
   par(nlev) = par_fraction * rad(nlev) * EXP( -(Kw+localext)*1e-6*dz(nlev) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = (nlev-1),1,-1
      localext_up = localext
      CALL aed_light_extinction(column, i, localext)

      par(i) = par(i+1) * EXP( -(Kw + localext_up) * dz(i+1) )

      IF (bioshade_feedback) extc(i) = Kw + localext
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
   IF (ALLOCATED(ws))         DEALLOCATE(ws)
   IF (ALLOCATED(max_))       DEALLOCATE(max_)
   IF (ALLOCATED(min_))       DEALLOCATE(min_)
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
