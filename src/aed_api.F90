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
          sub_mobility_t

   !#===========================================================#!
   TYPE api_config_t
      INTEGER  :: MaxLayers

      LOGICAL  :: mobility_off
      LOGICAL  :: bioshade_feedback
      LOGICAL  :: repair_state
      LOGICAL  :: do_plots
      LOGICAL  :: link_rain_loss
      LOGICAL  :: link_solar_shade
      LOGICAL  :: link_bottom_drag
      LOGICAL  :: ice

      INTEGER  :: split_factor
      INTEGER  :: ode_method
      INTEGER  :: benthic_mode

      AED_REAL :: rain_factor
      AED_REAL :: sw_factor
      AED_REAL :: friction

      AED_REAL :: Kw
      AED_REAL :: dt

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

   TYPE api_env_t
      AED_REAL,DIMENSION(:),POINTER :: temp      => null() !# temperature
      AED_REAL,DIMENSION(:),POINTER :: salt      => null() !# salinity
      AED_REAL,DIMENSION(:),POINTER :: rho       => null() !# density
      AED_REAL,DIMENSION(:),POINTER :: dz        => null() !# layer thickness
      AED_REAL,DIMENSION(:),POINTER :: height    => null() !# layer height (previously "h")
      AED_REAL,DIMENSION(:),POINTER :: area      => null() !# layer area
      AED_REAL,DIMENSION(:),POINTER :: depth     => null() !# layer_depth (previously "z")
      AED_REAL,DIMENSION(:),POINTER :: col_depth => null() !# col_depth (sheet - total depth of column)
      AED_REAL,DIMENSION(:),POINTER :: extc      => null() !# extinction coefficient
      AED_REAL,DIMENSION(:),POINTER :: tss       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss1       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss2       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss3       => null()
      AED_REAL,DIMENSION(:),POINTER :: ss4       => null()
      AED_REAL,DIMENSION(:),POINTER :: cvel      => null() !# cell velocity
      AED_REAL,DIMENSION(:),POINTER :: vvel      => null() !# vertical velocity
      AED_REAL,DIMENSION(:),POINTER :: bio_drag  => null()
      AED_REAL,DIMENSION(:),POINTER :: rad       => null()
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
      AED_REAL,DIMENSION(:),POINTER :: pres => null()

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

   !#===========================================================#!
   TYPE api_water_col_t
      TYPE(aed_column_t),POINTER :: column(:) => null()
      TYPE(api_env_t),POINTER :: env          => null()
   END TYPE api_water_col_t
   !#===========================================================#!
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

!  CHARACTER(len=64) :: NULCSTR = ""
   CHARACTER(len=80) :: cfg_fname = "none"

   AED_REAL :: dt
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
   LOGICAL :: do_plots = .TRUE.
!  LOGICAL :: do_zone_averaging = .FALSE.
   LOGICAL :: link_solar_shade = .TRUE.
   LOGICAL :: link_rain_loss = .FALSE.
!  LOGICAL :: depress_clutch = .FALSE.
   LOGICAL :: link_bottom_drag = .FALSE.
!  LOGICAL :: link_surface_drag = .FALSE.
!  LOGICAL :: link_water_density = .FALSE.
!  LOGICAL :: link_water_clarity = .FALSE.
!  LOGICAL :: link_ext_par = .FALSE.

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
   AED_REAL,DIMENSION(:),POINTER :: height    => null()
   AED_REAL,DIMENSION(:),POINTER :: area      => null()
   AED_REAL,DIMENSION(:),POINTER :: depth     => null()
   AED_REAL,DIMENSION(:),POINTER :: col_depth => null()
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
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: ws
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
!  LOGICAL :: request_nearest = .FALSE.
!  LOGICAL :: have_nearest = .FALSE.
!  LOGICAL :: reinited = .FALSE.
!  INTEGER :: ThisStep = 0
!  INTEGER :: n_cellids = 0

   !#------------------------------------------------------------
   !#   internal variables
   !#------------------------------------------------------------

   !# Integers storing number of variables being simulated
   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet

   CHARACTER(len=30),ALLOCATABLE,TARGET :: names(:)
   CHARACTER(len=30),ALLOCATABLE,TARGET :: bennames(:)
!  CHARACTER(len=30),ALLOCATABLE,TARGET :: diagnames(:)

   INTEGER :: MaxLayers = 0
   INTEGER :: zone_var = 0

   AED_REAL :: dt_eff
   LOGICAL :: reinited = .false.

  !-----------------------------------------------------------------------------
  INTERFACE

    SUBROUTINE sub_mobility_t(N,dt,h,A,ww,min_C,cc)
       INTEGER,INTENT(in)     :: N       !# number of vertical layers
       AED_REAL,INTENT(in)    :: dt      !# time step (s)
       AED_REAL,INTENT(in)    :: h(*)    !# layer thickness (m)
       AED_REAL,INTENT(in)    :: A(*)    !# layer areas (m2)
       AED_REAL,INTENT(in)    :: ww(*)   !# vertical speed (m/s)
       AED_REAL,INTENT(in)    :: min_C   !# minimum allowed cell concentration
       AED_REAL,INTENT(inout) :: cc(*)   !# cell concentration
    END SUBROUTINE sub_mobility_t

  END INTERFACE
  !-----------------------------------------------------------------------------


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
INTEGER FUNCTION aed_init_model(fname, NumWQ_Vars, NumWQ_Ben, NumWQ_Diag, NumWQ_DiagL)
!-------------------------------------------------------------------------------
! Initialize the AED driver by reading settings from "fname".
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: fname
   INTEGER,INTENT(out)     :: NumWQ_Vars, NumWQ_Ben
   INTEGER,INTENT(out)     :: NumWQ_Diag, NumWQ_DiagL
!
!LOCALS
   INTEGER :: status, namlst, i

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

   write(*,"(/,5X,'---------- AED API config : start ----------')")
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
   NumWQ_DiagL = n_vars_diag_sheet

   aed_init_model = n_aed_vars
END FUNCTION aed_init_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_config_model(conf)
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_config_t), INTENT(in) :: conf
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
!
   MaxLayers = conf%MaxLayers
   mobility_off = conf%mobility_off
   bioshade_feedback = conf%bioshade_feedback
   repair_state = conf%repair_state
   do_plots = conf%do_plots
   link_rain_loss = conf%link_rain_loss
   link_solar_shade = conf%link_solar_shade
   link_bottom_drag = conf%link_bottom_drag

   split_factor = conf%split_factor
!  ode_method = conf%ode_method
   benthic_mode = conf%benthic_mode

   rain_factor = conf%rain_factor
   sw_factor = conf%sw_factor
   friction = conf%friction

   Kw = conf%Kw
   dt = conf%dt
END SUBROUTINE aed_config_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_data(env)
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_data_t), INTENT(in) :: env
!
!LOCALS
   INTEGER :: av, v, sv, status
   TYPE(aed_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   cc         => env%cc
   cc_hz      => env%cc_hz
   cc_diag    => env%cc_diag
   cc_diag_hz => env%cc_diag_hz

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

   CALL aed_check_data

   write(*,"(/,5X,'----------  AED config : end  ----------',/)")

END SUBROUTINE aed_set_model_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_env(env)
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_env_t), INTENT(in) :: env
!
!LOCALS
   INTEGER :: tv
!-------------------------------------------------------------------------------
!BEGIN
!
   yearday  => env%yearday
   timestep => env%timestep

   longitude  => env%longitude
   latitude   => env%latitude

   temp         => env%temp
   salt         => env%salt
   rho          => env%rho
   dz           => env%dz
   height       => env%height
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

   IF (ASSOCIATED(yearday))   tv=aed_provide_sheet_global('yearday', 'yearday', 'day'    )
   IF (ASSOCIATED(timestep))  tv=aed_provide_sheet_global('timestep','timestep','seconds')

   IF (ASSOCIATED(longitude)) tv=aed_provide_sheet_global('longitude',     'longitude',         'radians'       )
   IF (ASSOCIATED(latitude))  tv=aed_provide_sheet_global('latitude',      'latitude',          'radians'       )

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

   IF (ASSOCIATED(mat_id))         tv=aed_provide_sheet_global('material',      'material',          '-'             )

   tv=aed_provide_global('nir',        'nir',                   'W/m2'   )
   tv=aed_provide_global('par',        'par',                   'W/m2'   )
   tv=aed_provide_global('uva',        'uva',                   'W/m2'   )
   tv=aed_provide_global('uvb',        'uvb',                   'W/m2'   )
END SUBROUTINE aed_set_model_env
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_check_data
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: av
   INTEGER :: v, d, sv, sd, ev, err_count
   TYPE(aed_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = .true.
            CASE ( 'salinity' )    ; tvar%found = .true.
            CASE ( 'density' )     ; tvar%found = .true.
            CASE ( 'layer_ht' )    ; tvar%found = .true.
            CASE ( 'extc_coef' )   ; tvar%found = .true.
            CASE ( 'tss' )         ; tvar%found = .true.
            CASE ( 'cell_vel' )    ; tvar%found = .true.
            CASE ( 'par' )         ; tvar%found = .true.
            CASE ( 'nir' )         ; tvar%found = .true.
            CASE ( 'uva' )         ; tvar%found = .true.
            CASE ( 'uvb' )         ; tvar%found = .true.
            CASE ( 'pressure' )    ; tvar%found = .true.
            CASE ( 'depth' )       ; tvar%found = .true.
            CASE ( 'sed_zone' )    ; tvar%found = .true.
            CASE ( 'wind_speed' )  ; tvar%found = .true.
            CASE ( 'par_sf' )      ; tvar%found = .true.
            CASE ( 'taub' )        ; tvar%found = .true.
            CASE ( 'col_depth' )   ; tvar%found = .true.
            CASE ( 'layer_area' )  ; tvar%found = .true.
            CASE ( 'rain' )        ; tvar%found = .true.
            CASE ( 'evap' )        ; tvar%found = .true.
            CASE ( 'air_temp' )    ; tvar%found = .true.
            CASE ( 'air_pres' )    ; tvar%found = .true.
            CASE ( 'humidity' )    ; tvar%found = .true.
            CASE ( 'longitude' )   ; tvar%found = .true.
            CASE ( 'latitude' )    ; tvar%found = .true.
            CASE ( 'yearday' )     ; tvar%found = .true.
            CASE ( 'timestep' )    ; tvar%found = .true.
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
END SUBROUTINE aed_check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, top, flux_pel, flux_atm, flux_ben)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in)  :: top
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_layers, n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av !, i
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed_variable_t),POINTER :: tvar
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
            CASE ( 'taub' )        ; column(av)%cell_sheet => layer_stress(1)
            CASE ( 'layer_area' )  ; column(av)%cell => area(:)
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(1)
            CASE ( 'evap' )        ; column(av)%cell_sheet => evap(1)
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp(1)
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => air_pres(1)
            CASE ( 'humidity' )    ; column(av)%cell_sheet => humidity(1)
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => col_depth(1)
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
SUBROUTINE aed_run_model(wlev, doMobility, doSurface)
!-------------------------------------------------------------------------------
!                        wlev is the number of levels used;
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: wlev
   PROCEDURE(sub_mobility_t),POINTER :: doMobility
   LOGICAL,INTENT(in) :: doSurface
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv

   AED_REAL :: min_C
   INTEGER  :: i, j, v, lev, split, zon=1

   TYPE (aed_column_t) :: column(n_aed_vars)
   TYPE (aed_column_t) :: column_sed(n_aed_vars)
   AED_REAL,TARGET :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben)
   AED_REAL,TARGET,ALLOCATABLE :: flux_pel(:, :)
   AED_REAL,TARGET :: flux_zon(aed_n_zones, n_vars+n_vars_ben)
   AED_REAL :: pa = 0.
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(flux_pel(MAX(wlev,aed_n_zones),n_vars+n_vars_ben))

   IF ( benthic_mode .GT. 1 ) THEN
      j = 1
      DO i=1,wlev
!        print *,'i =',i," j =",j
!        print *,'zone_heights',height(i),j,aedZones(j)%z_heights(1)
         IF (height(i) .GT. aedZones(j)%z_heights(1)) THEN
            sed_zones(i) = j * ( area(i) - pa ) / area(i)
            pa = area(i)
            j = j+1
         ELSE
            sed_zones(i) = j
         ENDIF
      ENDDO
   ENDIF

   CALL define_column(column, wlev, flux_pel, flux_atm, flux_ben)

   IF ( .NOT. reinited ) CALL re_initialize

   cc_diag = 0.
   cc_diag_hz = 0.

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
                  CALL doMobility(wlev, dt, dz, area, ws(:, i), min_C, cc(:, v))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   CALL check_states(column,wlev)

   DO split=1,split_factor

      IF (benthic_mode .GT. 1) THEN
         CALL p_copy_to_zone(aedZones, aed_n_zones, cc, cc_hz, cc_diag, cc_diag_hz, wlev)
         CALL p_calc_zone_areas(aedZones, aed_n_zones, area, wlev, height(wlev))
      ENDIF

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer
      CALL update_light(column, wlev)

      !# Fudge
      nir(:) = (par(:)/par_fraction) * nir_fraction
      uva(:) = (par(:)/par_fraction) * uva_fraction
      uvb(:) = (par(:)/par_fraction) * uvb_fraction

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
               aedZones(zon)%z_cc(1, v) = aedZones(zon)%z_cc(1, v)+ dt_eff*flux_zon(zon, v)
            ENDDO
         ENDDO

         !# Distribute cc-sed benthic properties back into main cc array
         CALL p_copy_from_zone(aedZones, aed_n_zones, cc, cc_hz, cc_diag, cc_diag_hz, wlev)
      ELSE
         DO v = n_vars+1, n_vars+n_vars_ben
            cc(1, v) = cc(1, v) + dt_eff*flux_ben(v)
         ENDDO
      ENDIF

      CALL check_states(column, wlev)
   ENDDO
   DEALLOCATE(flux_pel)

  ! IF ( display_minmax ) THEN
  !    v = 0; d = 0
  !    DO i=1,n_aed_vars
  !       IF ( aed_get_var(i, tv) ) THEN
  !          IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
  !             v = v + 1
  !             WRITE(*,'(1X,"VarLims: ",I0,1X,"<=> ",f15.8,f15.8," : ",A," (",A,")")') &
  !                                        v,MINVAL(cc(:,v)),MAXVAL(cc(:,v)),TRIM(tv%name),TRIM(tv%units)
  !             !print *,'VarLims',v,TRIM(tv%name),MINVAL(cc(v,:)),MAXVAL(cc(v,:))
  !          ELSE IF ( tv%diag  ) THEN
  !             d = d + 1
  !             WRITE(*,'(1X,"DiagLim: ",I0,1X,"<=> ",f15.8,f15.8," : ",A," (",A,")")') &
  !                                        d,MINVAL(cc_diag(:,d)),MAXVAL(cc_diag(:,d)),TRIM(tv%name),TRIM(tv%units)
  !             !print *,'DiagLim',d,TRIM(tv%name),MINVAL(cc_diag(d,:)),MAXVAL(cc_diag(d,:))
  !          ENDIF
  !       ENDIF
  !    ENDDO
  !  ENDIF

CONTAINS
!-------------------------------------------------------------------------------

   !############################################################################
   SUBROUTINE define_sed_column(column, zon, top, bot)
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
   !----------------------------------------------------------------------------
   !BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0

   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => aedZones(zon)%z_temp(:)
            CASE ( 'salinity' )    ; column(av)%cell => aedZones(zon)%z_salt(:)
            CASE ( 'density' )     ; column(av)%cell => aedZones(zon)%z_rho(:)
            CASE ( 'layer_ht' )    ; column(av)%cell => aedZones(zon)%z_dz(:)
            CASE ( 'extc_coef' )   ; column(av)%cell => aedZones(zon)%z_extc(:)
            CASE ( 'tss' )         ; column(av)%cell => aedZones(zon)%z_tss(:)
            CASE ( 'cell_vel' )    ; column(av)%cell => aedZones(zon)%z_vel(:)
            CASE ( 'par' )         ; column(av)%cell => aedZones(zon)%z_par(:)
            CASE ( 'nir' )         ; column(av)%cell => aedZones(zon)%z_nir(:)
            CASE ( 'uva' )         ; column(av)%cell => aedZones(zon)%z_uva(:)
            CASE ( 'uvb' )         ; column(av)%cell => aedZones(zon)%z_uvb(:)
            CASE ( 'pressure' )    ; column(av)%cell => aedZones(zon)%z_pres(:)
            CASE ( 'depth' )       ; column(av)%cell => aedZones(zon)%z_depth(:)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => aedZones(zon)%z_sed_zones(bot); zone_var = av
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd(1)
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0(1)
            CASE ( 'taub' )        ; column(av)%cell_sheet => layer_stress(1)
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => col_depth(1)
            CASE ( 'layer_area' )  ; column(av)%cell => aedZones(zon)%z_area(:)
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
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(bot, n_vars+sv)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => aedZones(zon)%z_cc(top, n_vars+sv)
            ENDIF
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
   END SUBROUTINE define_sed_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE re_initialize
   !----------------------------------------------------------------------------
   !ARGUMENTS
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
            CALL define_sed_column(column_sed, zon, 1, 1)

            aedZones(zon)%z_sed_zones(1) = zon  !MH TMP
!           print *,'aedZones(zon)%z_sed_zones',aedZones(zon)%z_sed_zones !MH TMP
            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_sed_zones(1)

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

            CALL aed_initialize_benthic(column_sed, zon)
         ENDDO
       ENDIF

      reinited = .TRUE.
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
      INTEGER, ALLOCATABLE :: layer_map(:)
      AED_REAL :: scale
!     AED_REAL, DIMENSION(wlev, n_vars+n_vars_ben)    :: flux_pel_pre
      AED_REAL, DIMENSION(:, :),ALLOCATABLE  :: flux_pel_pre
      AED_REAL, DIMENSION(aed_n_zones, n_vars+n_vars_ben) :: flux_pel_z
      AED_REAL :: localrainl, localshade, localdrag
      LOGICAL :: splitZone
      TYPE(aed_variable_t),POINTER :: tvar
   !----------------------------------------------------------------------------
   !BEGIN
   flux_pel = zero_
   flux_atm = zero_
   flux_ben = zero_
   flux_zon = zero_
   flux_pel_z = zero_

   ALLOCATE(flux_pel_pre(MAX(wlev, aed_n_zones), n_vars+n_vars_ben))
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
         CALL define_sed_column(column_sed, zon, 1, 1)

         !# Reinitialise flux_ben to be repopulated for this zone
         flux_ben = zero_
         flux_pel_pre = zero_

         !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
         IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_sed_zones(1)

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
            CALL aed_calculate_riparian(column_sed, zon, aedZones(zon)%z_pc_wet(1))
            IF (aedZones(zon)%z_pc_wet(1) .EQ. 0. ) CALL aed_calculate_dry(column_sed, zon)

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
         CALL aed_calculate_benthic(column_sed, zon, .TRUE.)

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
               splitZone = height(lev-1) < aedZones(zon-1)%z_heights(1)
            ELSE
               splitZone = ( 0.0 < aedZones(zon-1)%z_heights(1) )
            ENDIF
         ELSE
            splitZone = .FALSE.
         ENDIF

         IF (splitZone) THEN
            IF (lev .GT. 1) THEN
               scale = (aedZones(zon-1)%z_heights(1) - height(lev-1)) / (height(lev) - height(lev-1))
            ELSE
               scale = (aedZones(zon-1)%z_heights(1) - 0.0) / (height(lev) - 0.0)
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
      IF ( zone_var .GE. 1 ) column(zone_var)%cell_sheet => aedZones(1)%z_sed_zones(1)
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
            !# scaled to proportion of layer area that is "bottom"
            DO v=v_start,v_end
              IF ( cc(1, v) .GE. 0.0 ) flux_pel(lev, v) = max(-1.0 * cc(lev, v), flux_pel(lev, v)/dz(lev))
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
   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE aed_run_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: nlev
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
!  IF (ASSOCIATED(par))       DEALLOCATE(par)
!  IF (ASSOCIATED(nir))       DEALLOCATE(nir)
!  IF (ASSOCIATED(uva))       DEALLOCATE(uva)
!  IF (ASSOCIATED(uvb))       DEALLOCATE(uvb)
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
