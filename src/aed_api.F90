!###############################################################################
!#                                                                             #
!# aed_api.F90                                                                 #
!#                                                                             #
!# A generic interface between hist models and libaed-xxx                      #
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
!#  This file is part of libaed (Library for Aquatic Eco Dynamics)             #
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
!# CAB: redefining them.  Fluxes should be OK because data only persists in a
!# CAB: column run anyway

#include "aed_api.h"

#ifndef __GFORTRAN__
#  ifndef isnan
#    define isnan(x) ieee_is_nan(x)
#  endif
#endif

!# Temporary HOST specific flags
!# Ultimately these will go - the code should be sufficiently generic
!# to handle all cases, but while merging is happening we probably need them
#define GLM_VARIANT   1
#define TFV_VARIANT   2
#define SCH_VARIANT   3
#define ELC_VARIANT   4

#define CUR_VARIANT   GLM_VARIANT
!#define CUR_VARIANT   SCH_VARIANT


!-------------------------------------------------------------------------------
MODULE aed_api
!
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE aed_zones
   USE aed_ptm, ONLY : Particles, aed_calculate_particles

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
   !* A structure to pass coupling configuration values to AED  *!
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
      LOGICAL  :: do_particle_bgc

      INTEGER  :: split_factor = 1
      INTEGER  :: benthic_mode = 1

      AED_REAL :: rain_factor  = 1.
      AED_REAL :: sw_factor    = 1.
      AED_REAL :: friction     = 1.

      AED_REAL :: Kw
      AED_REAL :: Ksed

      AED_REAL :: nir_fraction = 0.520  ! 0.51
      AED_REAL :: par_fraction = 0.430  ! 0.45
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
      INTEGER                       :: n_layers

      !# time and location
      AED_REAL,POINTER              :: timestep       => null()
      AED_REAL,POINTER              :: yearday        => null()
      AED_REAL,POINTER              :: longitude      => null()
      AED_REAL,POINTER              :: latitude       => null()
      INTEGER, POINTER              :: col_num        => null()

      !# a columns active (wet or dry) status
      LOGICAL, POINTER              :: active         => null()

      !# above water (meteorological) conditions
      AED_REAL,POINTER              :: longwave       => null() !# longwave
      AED_REAL,POINTER              :: air_temp       => null() !# air_temp
      AED_REAL,POINTER              :: air_pres       => null() !# air_pres
      AED_REAL,POINTER              :: humidity       => null() !# humidity
      AED_REAL,POINTER              :: wind           => null() !# wind_speed
      AED_REAL,POINTER              :: rain           => null() !# rain
      AED_REAL,POINTER              :: evap           => null() !# evap
      AED_REAL,POINTER              :: I_0            => null() !# rad_sf

      !# water column depths and layer information
      AED_REAL,POINTER              :: col_depth      => null() !# col_depth (sheet - total depth of column)
      AED_REAL,POINTER              :: col_area       => null() !# col_area (sheet - area of column)
      AED_REAL,DIMENSION(:),POINTER :: height         => null() !# layer height (previously "h")
      AED_REAL,DIMENSION(:),POINTER :: depth          => null() !# layer_depth (previously "z")
      AED_REAL,DIMENSION(:),POINTER :: area           => null() !# area :layer area
      AED_REAL,DIMENSION(:),POINTER :: dz             => null() !# thick : layer thickness (previously "layer_ht")

      !# water column hydrodynamic information
      AED_REAL,DIMENSION(:),POINTER :: temp           => null() !# temperature
      AED_REAL,DIMENSION(:),POINTER :: salt           => null() !# salinity
      AED_REAL,DIMENSION(:),POINTER :: cvel           => null() !# velocity
      AED_REAL,DIMENSION(:),POINTER :: pres           => null() !# pressure
      AED_REAL,DIMENSION(:),POINTER :: rho            => null() !# density
      AED_REAL,DIMENSION(:),POINTER :: rad            => null() !# radiation

      !# water column light information
      AED_REAL,DIMENSION(:),POINTER :: extc           => null() !# extinction coefficient
      AED_REAL,DIMENSION(:),POINTER :: par            => null() !# if missing the following is calculated from I_0 (rad_sf) above
      AED_REAL,DIMENSION(:),POINTER :: nir            => null()
      AED_REAL,DIMENSION(:),POINTER :: uva            => null()
      AED_REAL,DIMENSION(:),POINTER :: uvb            => null()

      !# water column external host turbidity information
      AED_REAL,DIMENSION(:),POINTER :: tss            => null() !# total suspended solids
      AED_REAL,DIMENSION(:),POINTER :: ss1            => null()
      AED_REAL,DIMENSION(:),POINTER :: ss2            => null()
      AED_REAL,DIMENSION(:),POINTER :: ss3            => null()
      AED_REAL,DIMENSION(:),POINTER :: ss4            => null()

      !# bottom wave and current stress information
      AED_REAL,DIMENSION(:),POINTER :: ustar_bed      => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_uorb        => null()
      AED_REAL,DIMENSION(:),POINTER :: wv_t           => null()
      AED_REAL,POINTER              :: layer_stress   => null()
      AED_REAL,POINTER              :: dz_benthic     => null()

      !# bottom sediment zone classification
      AED_REAL,DIMENSION(:),POINTER :: sed_zones      => null()  !# sedzones are an odd mix - for GLM a zone will be different at
      AED_REAL,POINTER :: sed_zone     => null()  !# different levels - while for others a column would be all one zone
      AED_REAL,POINTER :: mat_id       => null()

      !# inforation for dry column / groundwater calculations
      AED_REAL,POINTER              :: bathy          => null()
      AED_REAL,POINTER              :: datum          => null()
      AED_REAL,POINTER              :: col_height     => null()
      AED_REAL,POINTER              :: nearest_height => null()
      INTEGER, POINTER              :: nearest_active => null()

      !# arrays for AED to feedback to hosts environment
      AED_REAL,DIMENSION(:),POINTER :: biodrag        => null()
      AED_REAL,DIMENSION(:),POINTER :: bioextc        => null()
      AED_REAL,POINTER              :: solarshade     => null()
      AED_REAL,POINTER              :: windshade      => null()
      AED_REAL,POINTER              :: rainloss       => null()
   END TYPE aed_env_t
   !#===========================================================#!

   !#===========================================================#!
   !* The data structure defining an AED water column           *!
   !*-----------------------------------------------------------*!
   TYPE api_col_data_t
      INTEGER :: n_layers           = 0 !# number of layers in this column
                                        !# in cases like GLM this may vary each timestep

      !# Arrays storing/pointing to the AED module state and diagnostic variables
      AED_REAL,DIMENSION(:,:),POINTER :: cc             => null()  !# (n_vars, n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_hz          => null()  !# (n_ben_vars)
      AED_REAL,DIMENSION(:,:),POINTER :: cc_diag        => null()  !# (n_diag_vars, n_layers)
      AED_REAL,DIMENSION(:),  POINTER :: cc_diag_hz     => null()  !# (n_diag_ben_vars)

      !# Column index, location and status
      AED_REAL,POINTER                :: yearday        => null()
      AED_REAL,POINTER                :: longitude      => null()
      AED_REAL,POINTER                :: latitude       => null()
      AED_REAL                        :: col_num
     !INTEGER, POINTER                :: col_num        => null()
      LOGICAL,POINTER                 :: active         => null()

      !# Arrays storing/pointing to surface (above water) environment
      AED_REAL,POINTER                :: longwave       => null()
      AED_REAL,POINTER                :: air_temp       => null()
      AED_REAL,POINTER                :: air_pres       => null()
      AED_REAL,POINTER                :: humidity       => null()
      AED_REAL,POINTER                :: wind           => null()
      AED_REAL,POINTER                :: rain           => null()
      AED_REAL,POINTER                :: evap           => null()
      AED_REAL,POINTER                :: I_0            => null()

      !# Arrays storing/pointing to the grid/domain information
      AED_REAL,POINTER                :: col_depth      => null() !# total depth of column (sheet)
      AED_REAL,POINTER                :: col_area       => null() !# area of column (sheet)
      AED_REAL,DIMENSION(:),POINTER   :: lheights       => null() !# top height of layers (n_layers)
      AED_REAL,DIMENSION(:),POINTER   :: depth          => null() !# depth at top of layers (n_layers)
      AED_REAL,DIMENSION(:),POINTER   :: area           => null() !# area of layers (n_layers)
      AED_REAL,DIMENSION(:),POINTER   :: dz             => null() !# vertical thickness of layers (n_layers)

      !# Arrays storing/pointing to water column environment
      AED_REAL,DIMENSION(:),POINTER   :: temp           => null()
      AED_REAL,DIMENSION(:),POINTER   :: salt           => null()
      AED_REAL,DIMENSION(:),POINTER   :: cvel           => null()
      AED_REAL,DIMENSION(:),POINTER   :: pres           => null()
      AED_REAL,DIMENSION(:),POINTER   :: rho            => null()
      AED_REAL,DIMENSION(:),POINTER   :: rad            => null()

      AED_REAL,DIMENSION(:),POINTER   :: extc           => null()
      AED_REAL,DIMENSION(:),POINTER   :: par            => null()
      AED_REAL,DIMENSION(:),POINTER   :: nir            => null()
      AED_REAL,DIMENSION(:),POINTER   :: uva            => null()
      AED_REAL,DIMENSION(:),POINTER   :: uvb            => null()

      AED_REAL,DIMENSION(:),POINTER   :: tss            => null()
      AED_REAL,DIMENSION(:),POINTER   :: ss1            => null()
      AED_REAL,DIMENSION(:),POINTER   :: ss2            => null()
      AED_REAL,DIMENSION(:),POINTER   :: ss3            => null()
      AED_REAL,DIMENSION(:),POINTER   :: ss4            => null()

      !# Arrays storing/pointing to benthic (bottom) environment
      AED_REAL,DIMENSION(:),POINTER   :: ustar_bed      => null()
      AED_REAL,DIMENSION(:),POINTER   :: wv_uorb        => null()
      AED_REAL,DIMENSION(:),POINTER   :: wv_t           => null()
      AED_REAL,POINTER                :: layer_stress   => null()
      AED_REAL,POINTER                :: dz_benthic     => null()

      AED_REAL,DIMENSION(:),POINTER   :: sed_zones      => null()
      AED_REAL,POINTER                :: sed_zone       => null()
      AED_REAL,POINTER                :: mat_id         => null()

      !# Arrays storing/pointing to riparian environment
      AED_REAL,POINTER                :: bathy          => null()
      AED_REAL,POINTER                :: datum          => null()
      AED_REAL,POINTER                :: col_height     => null()
      AED_REAL,POINTER                :: nearest_height => null()
      INTEGER, POINTER                :: nearest_active => null()

      !# Arrays storing/pointing information for feedback to host
      AED_REAL,DIMENSION(:),POINTER   :: biodrag        => null()
      AED_REAL,DIMENSION(:),POINTER   :: bioextc        => null()
      AED_REAL,POINTER                :: solarshade     => null()
      AED_REAL,POINTER                :: windshade      => null()
      AED_REAL,POINTER                :: rainloss       => null()
   END TYPE api_col_data_t
   !#===========================================================#!

   !#===========================================================#!
   !* A structure defining anviromental variable.               *!
   !*-----------------------------------------------------------*!
   TYPE api_env_def_t
      TYPE(aed_variable_t),POINTER :: tvar
      AED_REAL,POINTER :: data      => null()
      AED_REAL,DIMENSION(:),POINTER :: datac => null()
      AED_REAL,POINTER :: zdata      => null()
      AED_REAL,DIMENSION(:),POINTER :: zdatac => null()
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
   LOGICAL :: do_particle_bgc = .FALSE.

   LOGICAL :: bottom_one = .TRUE.

   !-------------------------------------------------------------
   !# External variables
   !-------------------------------------------------------------
   TYPE(api_col_data_t),DIMENSION(:),ALLOCATABLE,TARGET :: data
   !# Arrays for environmental variables (used if they are not supplied externally)
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: lpar
!  AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: matz

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

   AED_REAL,DIMENSION(:),  POINTER :: bathy      => null()

   AED_REAL,DIMENSION(:,:),POINTER :: biodrag    => null()
   AED_REAL,DIMENSION(:,:),POINTER :: bioextc    => null()
   AED_REAL,DIMENSION(:),  POINTER :: solarshade => null()
   AED_REAL,DIMENSION(:),  POINTER :: windshade  => null()
   AED_REAL,DIMENSION(:),  POINTER :: rainloss   => null()

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
   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_ptm_vars
   INTEGER :: n_vars_diag, n_vars_diag_sheet

   CHARACTER(len=48),ALLOCATABLE :: names(:)
   CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)

   INTEGER :: MaxLayers = 0
   INTEGER :: zone_var = 0
   INTEGER :: call_count = 0

   AED_REAL :: dt_eff
   LOGICAL :: reinited = .FALSE.

   TYPE(aed_column_t),DIMENSION(:,:),ALLOCATABLE,TARGET :: all_cols !# (n_aed_vars, ncols)
   TYPE(aed_column_t),DIMENSION(:,:),ALLOCATABLE,TARGET :: zon_cols !# (n_aed_vars, nzones)
   AED_REAL,ALLOCATABLE,TARGET :: flux_ben(:)          !# (n_vars+n_vars_ben)
   AED_REAL,ALLOCATABLE,TARGET :: flux_atm(:)          !# (n_vars+n_vars_ben)
   AED_REAL,ALLOCATABLE,TARGET :: flux_rip(:)          !# (n_vars+n_vars_ben)
   AED_REAL,ALLOCATABLE,TARGET :: flux_zon(:,:)        !# (n_vars+n_vars_ben, aed_n_zones)
   AED_REAL,ALLOCATABLE,TARGET :: flux_pel(:,:)        !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_pre !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_z   !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))

!  !# 2D env vars (sheet)
!  ( 'timestep',      'timestep',               'seconds'        )
!  ( 'yearday',       'yearday',                'day'            )
!  ( 'longitude',     'longitude',              'radians'        )
!  ( 'latitude',      'latitude',               'radians'        )
!  ( 'col_num',       'column number',          ''               )
!
!  ( 'longwave',      'longwave',               'W/m2'           )
!  ( 'air_temp',      'air temperature',        'celsius'        )
!  ( 'air_pres',      'air pressure',           'Pa'             )
!  ( 'humidity',      'relative humidity',      ''               )
!  ( 'wind_speed',    'wind speed',             'm/s'            )
!  ( 'rain',          'rainfall',               'm/s'            )
!  ( 'evap',          'evaporation',            'm/s'            )
!
!  ( 'par_sf',        'par_sf',                 'W/m2'           )
!  ( 'col_depth',     'column water depth',     'm above bottom' )
!  ( 'taub',          'layer stress',           'N/m2'           )
!  ( 'nearest_active','nearest active',         ''               )
!  ( 'nearest_depth', 'nearest depth',          'm'              )
!  ( 'sed_zone',      'current sediment zone',  ''               )
!  ( 'material',      'material',               ''               )
!  ( 'bathy',         'bathy',                  'm above datum'  )
!  ( 'biodrag',       'biodrag',                ''               )
!  ( 'bioextc',       'bioextc',                ''               )
!  ( 'solarshade',    'solarshade',             ''               )
!  ( 'windshade',     'windshade',              ''               )
!  ( 'rainloss',      'rain loss',              'm/s'            )
!
!  !# 3D env vars
!  ( 'temperature',   'temperature',            'celsius'        )
!  ( 'salinity',      'salinity',               'g/kg'           )
!  ( 'density',       'density',                'kg/m3'          )
!  ( 'layer_ht',      'layer heights',          'm'              )
!  ( 'layer_area',    'layer area',             'm2'             )
!  ( 'depth',         'depth',                  'm'              )
!  ( 'extc_coef',     'extinction coefficient', '/m'             )
!  ( 'tss',           'tss',                    'g/m3'           )
!  ( 'ss1',           'ss1',                    'g/m3'           )
!  ( 'ss2',           'ss2',                    'g/m3'           )
!  ( 'ss3',           'ss3',                    'g/m3'           )
!  ( 'ss4',           'ss4',                    'g/m3'           )
!  ( 'cell_vel',      'cell velocity',          'm/s'            )
!  ( 'pressure',      'pressure',               ''               )
!  ( 'sed_zones',     'sediment zones',         ''               )
!  ( 'nir',           'nir',                    'W/m2'           )
!  ( 'par',           'par',                    'W/m2'           )
!  ( 'uva',           'uva',                    'W/m2'           )
!  ( 'uvb',           'uvb',                    'W/m2'           )

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
! Print summary (names and attributes) of configured AED variables
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
         IF ( .NOT. tvar%sheet .AND. tvar%var_type == V_STATE ) THEN
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
         IF ( tvar%sheet .AND. tvar%var_type == V_STATE ) THEN
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
         IF ( tvar%var_type == V_DIAGNOSTIC .AND. .NOT. tvar%sheet ) THEN
            j = j + 1
            print "(7X,'D(',I4,') water column diagnostic   > ',A)",j , TRIM(tvar%name)
            !print *,"     D(",j,") AED diagnostic 3Dvariable: ", TRIM(tvar%name)
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_DIAGNOSTIC .AND. tvar%sheet ) THEN
            j = j + 1
            !print *,"     D(",j,") AED diagnostic 2Dvariable: ", TRIM(tvar%name)
            print "(7X,'D(',I4,') bottom/surface diagnostic ~ ',A)",j , TRIM(tvar%name)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE aed_show_vars
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed_configure_models(fname, NumWQ_Vars, NumWQ_Ben, NumWQ_Diag, NumWQ_DiagS, NumPTM_Vars)
!-------------------------------------------------------------------------------
! Initializes the AED-API driver by reading settings from "fname"
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: fname
   INTEGER,INTENT(out) :: NumWQ_Vars, NumWQ_Ben, NumPTM_Vars
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

   n_aed_vars = aed_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet, n_ptm_vars)

#if DEBUG
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         print *,"AED var ", i, tvar%sheet, tvar%var_type, TRIM(tvar%name)
      ELSE
         print *,"AED var ", i, " is empty"
      ENDIF
   ENDDO
#endif

   print "(/,5X,'AED : n_aed_vars  = ',I3,' ; n_vars            = ',I4)",n_aed_vars,n_vars
   print "(  5X,'AED : n_vars_ben  = ',I3,' ; n_ptm_vars        = ',I4)",n_vars_ben,n_ptm_vars
   print "(  5X,'AED : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I4,/)",n_vars_diag,n_vars_diag_sheet

   print "(/,5X,'AED : MaxLayers         = ',I4)",MaxLayers

   !# names = grab the names from info
   ALLOCATE(names(n_vars),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (names)'
   ALLOCATE(bennames(n_vars_ben),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (bennames)'

   NumWQ_Vars  = n_vars
   NumWQ_Ben   = n_vars_ben
   NumWQ_Diag  = n_vars_diag
   NumWQ_DiagS = n_vars_diag_sheet
   NumPTM_Vars = n_ptm_vars

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
! Routine to set the coupling information provided by the host
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
   do_particle_bgc = conf%do_particle_bgc

   split_factor = conf%split_factor
   IF (split_factor == 0) split_factor = 1
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
! Routine to set and then initialise AED (state) variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: ncols, nlevs
   TYPE(aed_data_t),INTENT(in) :: dat(ncols)
!
!LOCALS
   INTEGER :: av, v, sv, status, col, zon
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (nlevs > MaxLayers) MaxLayers = nlevs

   IF (.NOT. ALLOCATED(all_cols) ) THEN
      ALLOCATE(all_cols(n_aed_vars, ncols),stat=status)
      IF (status /= 0) STOP 'allocate_memory(): Error allocating "all_cols"'
   ENDIF
   IF (.NOT. ALLOCATED(zon_cols) ) THEN
      ALLOCATE(zon_cols(n_aed_vars, aed_n_zones),stat=status)
      IF (status /= 0) STOP 'allocate_memory(): Error allocating "zon_cols"'
   ENDIF

   IF (.NOT.ALLOCATED(flux_zon))     ALLOCATE(flux_zon(n_aed_vars, aed_n_zones))
   IF (.NOT.ALLOCATED(flux_ben))     ALLOCATE(flux_ben(n_aed_vars))
   IF (.NOT.ALLOCATED(flux_atm))     ALLOCATE(flux_atm(n_aed_vars))
   IF (.NOT.ALLOCATED(flux_rip))     ALLOCATE(flux_rip(n_aed_vars))

   IF (.NOT.ALLOCATED(flux_pel))     ALLOCATE(flux_pel(    n_aed_vars, MAX(MaxLayers, aed_n_zones)))
   IF (.NOT.ALLOCATED(flux_pel_pre)) ALLOCATE(flux_pel_pre(n_aed_vars, MAX(MaxLayers, aed_n_zones)))
   IF (.NOT.ALLOCATED(flux_pel_z))   ALLOCATE(flux_pel_z(  n_aed_vars, MAX(MaxLayers, aed_n_zones)))

   !# Should have been done in the env setup
   IF (.NOT. ALLOCATED(data) ) ALLOCATE(data(ncols))
   !----------------------------------------------------------------------------
   !# Set pointers to state variable and diagnist variables arrays
   DO col=1,ncols
      data(col)%cc         => dat(col)%cc
!     data(col)%cc_hz      => dat(col)%cc_hz          !# Eventually this will be properly separated
      data(col)%cc_hz      => dat(col)%cc(n_vars:n_vars+n_vars_ben,1)
      data(col)%cc_diag    => dat(col)%cc_diag
      data(col)%cc_diag_hz => dat(col)%cc_diag_hz
   ENDDO

   ALLOCATE(min_((n_vars + n_vars_ben))) ; ALLOCATE(max_((n_vars + n_vars_ben)))

   CALL aed_show_vars

   !----------------------------------------------------------------------------
   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "     ERROR getting variable info"
      IF ( tvar%var_type == V_STATE ) THEN  !# state variable
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

   DO col=1, nCols
      CALL define_column(all_cols(:,col), col)
   ENDDO
   DO zon=1,aed_n_zones
      CALL define_zone_column(zon_cols(:,zon), zon)
   ENDDO
END SUBROUTINE aed_set_model_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_model_env(env, nlayrs, ncols)
!-------------------------------------------------------------------------------
! Routine to set environment (external) variables to main column "data" array
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nlayrs, ncols
   TYPE(aed_env_t),INTENT(in) :: env(ncols)
!
!LOCALS
   INTEGER :: tv, col, status = 0
   LOGICAL :: no_biodg = .FALSE., no_bioex = .FALSE.
   LOGICAL :: no_sshad = .FALSE., no_wshad = .FALSE.
   LOGICAL :: no_rianl = .FALSE., no_bathy = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (MaxLayers < nlayrs) MaxLayers = nlayrs

   !# Allocate main AED column arrays
   IF (.NOT. ALLOCATED(data) ) THEN
      ALLOCATE(data(ncols),stat=status)
      IF (status /= 0) STOP 'allocate_memory(): Error allocating "data"'
   ENDIF

   !# Check whether external light field is to be used, or allocate locally
   IF (.NOT. link_ext_par) ALLOCATE(lpar(MaxLayers,ncols))
!  ALLOCATE(matz(ncols)) ; matz(:) = zero_

   !# Set effective time/step
   dt_eff = env(1)%timestep/FLOAT(split_factor)  ! was timestep/FLOAT(split_factor)

   !# Check for "optional" environment/feedback vars that were not provided
   !  and allocate locally
   no_bathy = (.NOT.ASSOCIATED(env(1)%bathy))
   no_biodg = (.NOT.ASSOCIATED(env(1)%biodrag))
   no_bioex = (.NOT.ASSOCIATED(env(1)%bioextc))
   no_sshad = (.NOT.ASSOCIATED(env(1)%solarshade))
   no_wshad = (.NOT.ASSOCIATED(env(1)%windshade))
   no_rianl = (.NOT.ASSOCIATED(env(1)%rainloss))

   IF (no_bathy) THEN ; ALLOCATE(bathy(ncols))             ; bathy = zero_      ; ENDIF
   IF (no_biodg) THEN ; ALLOCATE(biodrag(MaxLayers,ncols)) ; biodrag = zero_    ; ENDIF
   IF (no_bioex) THEN ; ALLOCATE(bioextc(MaxLayers,ncols)) ; bioextc = zero_    ; ENDIF
   IF (no_sshad) THEN ; ALLOCATE(solarshade(ncols))        ; solarshade = zero_ ; ENDIF
   IF (no_wshad) THEN ; ALLOCATE(windshade(ncols))         ; windshade = zero_  ; ENDIF
   IF (no_rianl) THEN ; ALLOCATE(rainloss(ncols))          ; rainloss = zero_   ; ENDIF

   !# Set local AED column data structure to point to environment vars from host
   timestep => env(1)%timestep
   yearday  => env(1)%yearday

   DO col=1,ncols

      data(col)%n_layers     =  env(col)%n_layers

      data(col)%yearday      => env(col)%yearday
      data(col)%longitude    => env(col)%longitude
      data(col)%latitude     => env(col)%latitude
      data(col)%col_num      =  FLOAT(col)  !  => env(col)%col_num

      data(col)%active       => env(col)%active

      data(col)%longwave     => env(col)%longwave
      data(col)%air_temp     => env(col)%air_temp
      data(col)%air_pres     => env(col)%air_pres
      data(col)%humidity     => env(col)%humidity
      data(col)%wind         => env(col)%wind
      data(col)%rain         => env(col)%rain
      data(col)%evap         => env(col)%evap
      data(col)%I_0          => env(col)%I_0

      data(col)%col_depth    => env(col)%col_depth
      data(col)%col_area     => env(col)%col_area
      data(col)%lheights     => env(col)%height
      data(col)%depth        => env(col)%depth
      data(col)%area         => env(col)%area
      data(col)%dz           => env(col)%dz

      data(col)%temp         => env(col)%temp
      data(col)%salt         => env(col)%salt
      data(col)%cvel         => env(col)%cvel
      data(col)%pres         => env(col)%pres
      data(col)%rho          => env(col)%rho
      data(col)%rad          => env(col)%rad

      data(col)%extc         => env(col)%extc
      IF (link_ext_par) THEN
         data(col)%par       => env(col)%par
      ELSE
         data(col)%par       => lpar(:,col)
      ENDIF
      data(col)%nir          => env(col)%nir
      data(col)%uva          => env(col)%uva
      data(col)%uvb          => env(col)%uvb

      data(col)%tss => env(col)%tss
      data(col)%ss1 => env(col)%ss1
      data(col)%ss2 => env(col)%ss2
      data(col)%ss3 => env(col)%ss3
      data(col)%ss4 => env(col)%ss4

      data(col)%ustar_bed    => env(col)%ustar_bed
      data(col)%wv_uorb      => env(col)%wv_uorb
      data(col)%wv_t         => env(col)%wv_t
      data(col)%layer_stress => env(col)%layer_stress
      data(col)%dz_benthic   => env(col)%dz_benthic

      data(col)%sed_zones    => env(col)%sed_zones
      data(col)%sed_zone     => env(col)%sed_zone
      data(col)%mat_id       => env(col)%mat_id

      IF (no_bathy) THEN ; data(col)%bathy      => bathy(col)
      ELSE ;               data(col)%bathy      => env(col)%bathy      ; ENDIF

      data(col)%datum => env(col)%datum
      data(col)%col_height => env(col)%col_height
      data(col)%nearest_height => env(col)%nearest_height
      data(col)%nearest_active => env(col)%nearest_active

      IF (no_biodg) THEN ; data(col)%biodrag    => biodrag(:,col)
      ELSE ;               data(col)%biodrag    => env(col)%biodrag    ; ENDIF
      IF (no_bioex) THEN ; data(col)%bioextc    => bioextc(:,col)
      ELSE ;               data(col)%bioextc    => env(col)%bioextc    ; ENDIF
      IF (no_sshad) THEN ; data(col)%solarshade => solarshade(col)
      ELSE ;               data(col)%solarshade => env(col)%solarshade ; ENDIF
      IF (no_wshad) THEN ; data(col)%windshade  => windshade(col)
      ELSE ;               data(col)%windshade  => env(col)%windshade  ; ENDIF
      IF (no_rianl) THEN ; data(col)%rainloss   => rainloss(col)
      ELSE ;               data(col)%rainloss   => env(col)%rainloss   ; ENDIF
   ENDDO

   !# Register module accesible environment variables to AED core variable list
   IF (ASSOCIATED(timestep))       tv=aed_provide_sheet_global('timestep',      'timestep',              'seconds' )
   IF (ASSOCIATED(yearday))        tv=aed_provide_sheet_global('yearday',       'yearday',               'day'     )
   IF (BSSOCIATED(longitude))      tv=aed_provide_sheet_global('longitude',     'longitude',             'radians' )
   IF (BSSOCIATED(latitude))       tv=aed_provide_sheet_global('latitude',      'latitude',              'radians' )
                                   tv=aed_provide_sheet_global('col_num',       'column number',         '-'       )

   IF (BSSOCIATED(longwave))       tv=aed_provide_sheet_global('longwave',      'longwave',              'W/m2'    )
   IF (BSSOCIATED(air_temp))       tv=aed_provide_sheet_global('air_temp',      'air temperature',       'celsius' )
   IF (BSSOCIATED(air_pres))       tv=aed_provide_sheet_global('air_pres',      'air pressure',          'Pa'      )
   IF (BSSOCIATED(humidity))       tv=aed_provide_sheet_global('humidity',      'relative humidity',     '-'       )
   IF (BSSOCIATED(wind))           tv=aed_provide_sheet_global('wind_speed',    'wind speed',            'm/s'     )
   IF (BSSOCIATED(rain))           tv=aed_provide_sheet_global('rain',          'rainfall',              'm/s'     )
   IF (BSSOCIATED(evap))           tv=aed_provide_sheet_global('evap',          'evaporation',           'm/s'     )
   IF (BSSOCIATED(I_0))            tv=aed_provide_sheet_global('par_sf',        'par_sf',                'W/m2'    )

   IF (BSSOCIATED(col_depth))      tv=aed_provide_sheet_global('col_depth',     'column water depth',    'm above bottom')
   IF (BSSOCIATED(depth))          tv=aed_provide_global      ('depth',         'depth',                 'm'       )
   IF (BSSOCIATED(area))           tv=aed_provide_global      ('layer_area',    'layer area',            'm2'      )
   IF (BSSOCIATED(dz))             tv=aed_provide_global      ('layer_ht',      'layer heights',         'm'       )

   IF (BSSOCIATED(temp))           tv=aed_provide_global      ('temperature',   'temperature',           'celsius' )
   IF (BSSOCIATED(salt))           tv=aed_provide_global      ('salinity',      'salinity',              'g/kg'    )
   IF (BSSOCIATED(cvel))           tv=aed_provide_global      ('cell_vel',      'cell velocity',         'm/s'     )
   IF (BSSOCIATED(pres))           tv=aed_provide_global      ('pressure',      'pressure',              ''        )
   IF (BSSOCIATED(rho))            tv=aed_provide_global      ('density',       'density',               'kg/m3'   )
   IF (BSSOCIATED(rad))            tv=aed_provide_global      ('rad',           'radiation',             'W/m2'    )

   IF (BSSOCIATED(extc))           tv=aed_provide_global      ('extc_coef',     'extinction coefficient','/m'      )
                                   tv=aed_provide_global      ('nir',           'nir',                   'W/m2'    )
                                   tv=aed_provide_global      ('par',           'par',                   'W/m2'    )
                                   tv=aed_provide_global      ('uva',           'uva',                   'W/m2'    )
                                   tv=aed_provide_global      ('uvb',           'uvb',                   'W/m2'    )

   IF (BSSOCIATED(tss))            tv=aed_provide_global      ('tss',           'tss',                   'g/m3'    )
   IF (BSSOCIATED(ss1))            tv=aed_provide_global      ('ss1',           'ss1',                   'g/m3'    )
   IF (BSSOCIATED(ss2))            tv=aed_provide_global      ('ss2',           'ss2',                   'g/m3'    )
   IF (BSSOCIATED(ss3))            tv=aed_provide_global      ('ss3',           'ss3',                   'g/m3'    )
   IF (BSSOCIATED(ss4))            tv=aed_provide_global      ('ss4',           'ss4',                   'g/m3'    )

   IF (BSSOCIATED(layer_stress))   tv=aed_provide_sheet_global('taub',          'layer stress',          'N/m2'    )

   IF (BSSOCIATED(sed_zones))      tv=aed_provide_global      ('sed_zones',     'sediment zones',        '-'       )
   IF (BSSOCIATED(sed_zone))       tv=aed_provide_sheet_global('sed_zone',      'current sediment zone', '-'       )
   IF (BSSOCIATED(mat_id))         tv=aed_provide_sheet_global('material',      'material',              '-'       )

                                   tv=aed_provide_sheet_global('bathy',         'bathy',                 'm above datum' )
   IF (ASSOCIATED(nearest_active)) tv=aed_provide_sheet_global('nearest_active','nearest active',        '-'       )
   IF (ASSOCIATED(nearest_depth))  tv=aed_provide_sheet_global('nearest_depth', 'nearest depth',         'm'       )

                                   tv=aed_provide_global      ('biodrag',       'biodrag',               ''        )
                                   tv=aed_provide_global      ('bioextc',       'bioextc',               ''        )
                                   tv=aed_provide_sheet_global('solarshade',    'solarshade',            ''        )
                                   tv=aed_provide_sheet_global('windshade',     'windshade',             ''        )
                                   tv=aed_provide_sheet_global('rainloss',      'rain loss',             'm/s'     )

   !# env vars currently not made available
    !active
    !col_area
    !height
    !ustar_bed
    !wv_uorb
    !wv_t
    !dz_benthic
    !bathy
    !datum
    !col_height

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
   INTEGER :: av, v, d, sv, sd, ev, err_count
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed_vars
      IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# Environment (external/global) variable
         ev = ev + 1
         SELECT CASE (tvar%name)

            CASE ( 'timestep' )    ; tvar%found = ASSOCIATED(timestep)
            CASE ( 'yearday' )     ; tvar%found = ASSOCIATED(yearday)
            CASE ( 'longitude' )   ; tvar%found = BSSOCIATED(longitude)
            CASE ( 'latitude' )    ; tvar%found = BSSOCIATED(latitude)
            CASE ( 'col_num' )     ; tvar%found = .TRUE.
           !CASE ( 'col_num' )     ; tvar%found = BSSOCIATED(col_num)

            CASE ( 'longwave' )    ; tvar%found = BSSOCIATED(longwave)
            CASE ( 'air_temp' )    ; tvar%found = BSSOCIATED(air_temp)
            CASE ( 'air_pres' )    ; tvar%found = BSSOCIATED(air_pres)
            CASE ( 'humidity' )    ; tvar%found = BSSOCIATED(humidity)
            CASE ( 'wind_speed' )  ; tvar%found = BSSOCIATED(wind)
            CASE ( 'rain' )        ; tvar%found = BSSOCIATED(rain)
            CASE ( 'evap' )        ; tvar%found = BSSOCIATED(evap)
            CASE ( 'par_sf' )      ; tvar%found = BSSOCIATED(I_0)

            CASE ( 'col_depth' )   ; tvar%found = BSSOCIATED(col_depth)
            CASE ( 'depth' )       ; tvar%found = BSSOCIATED(depth)
            CASE ( 'layer_area' )  ; tvar%found = BSSOCIATED(area)
            CASE ( 'layer_ht' )    ; tvar%found = BSSOCIATED(dz)

            CASE ( 'temperature' ) ; tvar%found = BSSOCIATED(temp)
            CASE ( 'salinity' )    ; tvar%found = BSSOCIATED(salt)
            CASE ( 'cell_vel' )    ; tvar%found = BSSOCIATED(cvel)
            CASE ( 'pressure' )    ; tvar%found = BSSOCIATED(pres)
            CASE ( 'density' )     ; tvar%found = BSSOCIATED(rho)
            CASE ( 'rad' )         ; tvar%found = BSSOCIATED(rad)

            CASE ( 'extc_coef' )   ; tvar%found = BSSOCIATED(extc)
            CASE ( 'par' )         ; tvar%found = BSSOCIATED(par)
            CASE ( 'nir' )         ; tvar%found = BSSOCIATED(nir)
            CASE ( 'uva' )         ; tvar%found = BSSOCIATED(uva)
            CASE ( 'uvb' )         ; tvar%found = BSSOCIATED(uvb)

            CASE ( 'tss' )         ; tvar%found = BSSOCIATED(tss)
            CASE ( 'ss1' )         ; tvar%found = BSSOCIATED(ss1)
            CASE ( 'ss2' )         ; tvar%found = BSSOCIATED(ss2)
            CASE ( 'ss3' )         ; tvar%found = BSSOCIATED(ss3)
            CASE ( 'ss4' )         ; tvar%found = BSSOCIATED(ss4)

            CASE ( 'taub' )        ; tvar%found = BSSOCIATED(layer_stress)

            CASE ( 'sed_zones' )   ; tvar%found = BSSOCIATED(sed_zones)
            CASE ( 'sed_zone' )    ; tvar%found = BSSOCIATED(sed_zone)
            CASE ( 'material' )    ; tvar%found = BSSOCIATED(mat_id)

            CASE ( 'bathy' )       ; tvar%found = BSSOCIATED(bathy)
            CASE ( 'nearest_depth' )  ; tvar%found = have_nearest ; request_nearest = have_nearest
            CASE ( 'nearest_active' ) ; tvar%found = have_nearest ; request_nearest = have_nearest

            CASE ( 'biodrag' )     ; tvar%found = BSSOCIATED(biodrag)
            CASE ( 'bioextc' )     ; tvar%found = BSSOCIATED(bioextc)
            CASE ( 'solarshade' )  ; tvar%found = BSSOCIATED(solarshade)
            CASE ( 'windshade' )   ; tvar%found = BSSOCIATED(windshade)
            CASE ( 'rainloss' )    ; tvar%found = BSSOCIATED(rainloss)

            CASE DEFAULT ; print*,"ERROR: external variable "//TRIM(tvar%name)//" not found. Check host compatibility."
                           err_count = err_count + 1
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
         ELSE
            d = d + 1
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN       !# State variable
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

   IF ( err_count > 0 ) CALL STOPIT("*** ERRORS IN CONFIGURATION")
END SUBROUTINE aed_check_model_setup
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(icolm, col)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: icolm(:)
   INTEGER, INTENT(in) :: col
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

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'timestep' )    ; icolm(av)%cell_sheet => timestep
            CASE ( 'yearday' )     ; icolm(av)%cell_sheet => data(col)%yearday   !yearday
            CASE ( 'longitude' )   ; icolm(av)%cell_sheet => data(col)%longitude
            CASE ( 'latitude' )    ; icolm(av)%cell_sheet => data(col)%latitude
            CASE ( 'col_num' )     ; icolm(av)%cell_sheet => data(col)%col_num

            CASE ( 'longwave' )    ; icolm(av)%cell_sheet => data(col)%longwave
            CASE ( 'air_temp' )    ; icolm(av)%cell_sheet => data(col)%air_temp
            CASE ( 'air_pres' )    ; icolm(av)%cell_sheet => data(col)%air_pres
            CASE ( 'humidity' )    ; icolm(av)%cell_sheet => data(col)%humidity
            CASE ( 'wind_speed' )  ; icolm(av)%cell_sheet => data(col)%wind
            CASE ( 'rain' )        ; icolm(av)%cell_sheet => data(col)%rain
            CASE ( 'evap' )        ; icolm(av)%cell_sheet => data(col)%evap
            CASE ( 'par_sf' )      ; icolm(av)%cell_sheet => data(col)%I_0

            CASE ( 'col_depth' )   ; icolm(av)%cell_sheet => data(col)%col_depth
            CASE ( 'depth' )       ; icolm(av)%cell => data(col)%depth
            CASE ( 'layer_area' )  ; icolm(av)%cell => data(col)%area
            CASE ( 'layer_ht' )    ; icolm(av)%cell => data(col)%dz

            CASE ( 'temperature' ) ; icolm(av)%cell => data(col)%temp
            CASE ( 'salinity' )    ; icolm(av)%cell => data(col)%salt
            CASE ( 'cell_vel' )    ; icolm(av)%cell => data(col)%cvel
            CASE ( 'pressure' )    ; icolm(av)%cell => data(col)%pres
            CASE ( 'density' )     ; icolm(av)%cell => data(col)%rho
            CASE ( 'rad' )         ; icolm(av)%cell => data(col)%rad

            CASE ( 'extc_coef' )   ; icolm(av)%cell => data(col)%extc
            CASE ( 'par' )         ; icolm(av)%cell => data(col)%par
            CASE ( 'nir' )         ; icolm(av)%cell => data(col)%nir
            CASE ( 'uva' )         ; icolm(av)%cell => data(col)%uva
            CASE ( 'uvb' )         ; icolm(av)%cell => data(col)%uvb

            CASE ( 'tss' )         ; icolm(av)%cell => data(col)%tss
            CASE ( 'ss1' )         ; icolm(av)%cell => data(col)%ss1   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss2' )         ; icolm(av)%cell => data(col)%ss2   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss3' )         ; icolm(av)%cell => data(col)%ss3   ! For FV API 2.0 (To be connected to sed_conc)
            CASE ( 'ss4' )         ; icolm(av)%cell => data(col)%ss4   ! For FV API 2.0 (To be connected to sed_conc)

            CASE ( 'taub' )        ; icolm(av)%cell_sheet => data(col)%layer_stress ! CAB? col_taub

            CASE ( 'sed_zones' )   ; icolm(av)%cell => data(col)%sed_zones
            CASE ( 'sed_zone' )    ; icolm(av)%cell_sheet => data(col)%sed_zone

     ! CAB: Not handled yet
     !      CASE ( 'material' )    ; IF ( do_zone_averaging ) THEN
     !                                  icolm(av)%cell_sheet => zone(zm(col))
     !                               ELSE
     !                                  icolm(av)%cell_sheet => mat(col)
     !                               ENDIF
            CASE ( 'material' )    ; icolm(av)%cell_sheet => data(col)%mat_id

            CASE ( 'bathy' )       ; icolm(av)%cell_sheet => data(col)%bathy
            CASE ( 'nearest_active' ) ; icolm(av)%cell_sheet => nearest_active(col)
            CASE ( 'nearest_depth' )  ; icolm(av)%cell_sheet => nearest_depth(col)

            CASE ( 'biodrag' )     ; icolm(av)%cell => data(col)%biodrag
            CASE ( 'bioextc' )     ; icolm(av)%cell => data(col)%bioextc
            CASE ( 'solarshade' )  ; icolm(av)%cell_sheet => data(col)%solarshade
            CASE ( 'windshade' )   ; icolm(av)%cell_sheet => data(col)%windshade
            CASE ( 'rainloss' )    ; icolm(av)%cell_sheet => data(col)%rainloss

            CASE DEFAULT ; CALL STOPIT("ERROR: environment variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd =    sd + 1
            icolm(av)%cell_sheet => data(col)%cc_diag_hz(sd)
         ELSE
            d = d + 1
            icolm(av)%cell => data(col)%cc_diag(d,:)
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            icolm(av)%cell_sheet => data(col)%cc_hz(sv)
            icolm(av)%flux_ben => flux_ben(n_vars+sv)
            icolm(av)%flux_atm => flux_atm(n_vars+sv)
            icolm(av)%flux_rip => flux_rip(n_vars+sv)
         ELSE
            v = v + 1
            icolm(av)%cell => data(col)%cc(v,:)
            icolm(av)%flux_pel => flux_pel(v,:)
            icolm(av)%flux_ben => flux_ben(v)
            icolm(av)%flux_atm => flux_atm(v)
            icolm(av)%flux_rip => flux_rip(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_zone_column(zcolm, zon)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_column_t),INTENT(inout) :: zcolm(n_aed_vars)
   INTEGER, INTENT(in)  :: zon
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
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'timestep' )    ; zcolm(av)%cell_sheet => timestep
            CASE ( 'yearday' )     ; zcolm(av)%cell_sheet => yearday
            CASE ( 'longitude' )   ; zcolm(av)%cell_sheet => aedZones(zon)%longitude
            CASE ( 'latitude' )    ; zcolm(av)%cell_sheet => aedZones(zon)%latitude
            CASE ( 'col_num' )     ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_col_num

            CASE ( 'longwave' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_longwave
            CASE ( 'air_temp' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_air_temp
            CASE ( 'air_pres' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_air_pres
            CASE ( 'humidity' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_humidity
            CASE ( 'wind_speed' )  ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_wind
            CASE ( 'rain' )        ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_rain
            CASE ( 'evap' )        ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_evap
            CASE ( 'par_sf' )      ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_I_0

            CASE ( 'col_depth' )   ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_col_depth
            CASE ( 'depth' )       ; zcolm(av)%cell => aedZones(:)%z_env%z_depth
            CASE ( 'layer_area' )  ; zcolm(av)%cell => aedZones(:)%z_env%z_area
            CASE ( 'layer_ht' )    ; zcolm(av)%cell => aedZones(:)%z_env%z_dz

            CASE ( 'temperature' ) ; zcolm(av)%cell => aedZones(:)%z_env%z_temp
            CASE ( 'salinity' )    ; zcolm(av)%cell => aedZones(:)%z_env%z_salt
            CASE ( 'cell_vel' )    ; zcolm(av)%cell => aedZones(:)%z_env%z_vel
            CASE ( 'pressure' )    ; zcolm(av)%cell => aedZones(:)%z_env%z_pres
            CASE ( 'density' )     ; zcolm(av)%cell => aedZones(:)%z_env%z_rho
            CASE ( 'rad' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_rad

            CASE ( 'extc_coef' )   ; zcolm(av)%cell => aedZones(:)%z_env%z_extc
            CASE ( 'par' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_par
            CASE ( 'nir' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_nir
            CASE ( 'uva' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_uva
            CASE ( 'uvb' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_uvb

            CASE ( 'tss' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_tss
            CASE ( 'ss1' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_ss1
            CASE ( 'ss2' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_ss2
            CASE ( 'ss3' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_ss3
            CASE ( 'ss4' )         ; zcolm(av)%cell => aedZones(:)%z_env%z_ss4

            CASE ( 'taub' )        ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_layer_stress !CAB ??? (bot)

            CASE ( 'sed_zones' )   ; zcolm(av)%cell => aedZones(:)%z_env%z_sed_zones; zone_var = av
            CASE ( 'sed_zone' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_sed_zone; zone_var = av
            CASE ( 'material' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_mat_id

            CASE ( 'bathy' )       ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_bathy
        !   CASE ( 'nearest_active' ) ; zcolm(av)%cell_sheet => nearest_active(zon)
        !   CASE ( 'nearest_depth' )  ; zcolm(av)%cell_sheet => nearest_depth(zon)

            CASE ( 'biodrag' )     ; zcolm(av)%cell => aedZones(:)%z_env%z_biodrag
            CASE ( 'bioextc' )     ; zcolm(av)%cell => aedZones(:)%z_env%z_bioextc
            CASE ( 'solarshade' )  ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_solarshade
            CASE ( 'windshade' )   ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_windshade
            CASE ( 'rainloss' )    ; zcolm(av)%cell_sheet => aedZones(zon)%z_env%z_rainloss

            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            zcolm(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
         ELSE
            d = d + 1
            zcolm(av)%cell => aedZones(zon)%z_cc_diag(d, :)
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            zcolm(av)%cell_sheet => aedZones(zon)%z_cc_hz(sv)
            zcolm(av)%flux_ben => flux_ben(n_vars+sv)
            zcolm(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            zcolm(av)%cell => aedZones(zon)%z_cc(v, :)
            zcolm(av)%flux_atm => flux_atm(v)
            zcolm(av)%flux_pel => flux_pel(v, :)
            zcolm(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_zone_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_states(icolm, col, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: icolm(:)
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
      CALL aed_equilibrate(icolm, lev)    !MH this should be in the main do_glm routine ????!!!
      v = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( tv%var_type == V_STATE ) THEN
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
SUBROUTINE aed_run_model(nCols, nLevs, doSurface)
!-------------------------------------------------------------------------------
!                        nLevs is the number of levels used;
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nCols, nLevs
   LOGICAL,INTENT(in) :: doSurface
!
!LOCALS
   INTEGER :: col, top, bot
!
!-------------------------------------------------------------------------------
!BEGIN
   !----------------------------------------------------------------------------
   !# index and time-step updates
   IF ( bottom_one ) THEN
      top = nLevs ; bot = 1
   ELSE
      top = 1 ; bot = nLevs
   ENDIF

   !# reset effective time/step
   dt_eff = timestep/FLOAT(split_factor)
   call_count = call_count + 1

   !----------------------------------------------------------------------------
   !# Resetting and re-initialisation tasks
   DO col=1, nCols
      data(col)%cc_diag = 0.
      data(col)%cc_diag_hz = 0.
   ENDDO

   IF ( .NOT. reinited ) THEN
      DO col=1, nCols
      !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
         CALL re_initialize(all_cols(:,col), nLevs)
      ENDDO
      reinited = .TRUE.
   ENDIF

   !----------------------------------------------------------------------------
   !# Pre flux integration tasks
   DO col=1, nCols
   !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
      CALL pre_kinetics(all_cols(:,col), col, nLevs)
   ENDDO

   !----------------------------------------------------------------------------
   !# Main time-step tasks
   DO col=1, nCols
   !  IF (.NOT. active(col)) CYCLE  !# skip this column if dry
      CALL aed_run_column(all_cols(:,col), col, nLevs, doSurface)
   ENDDO

   !----------------------------------------------------------------------------
   !# Particle tracking tasks
   IF (do_particle_bgc) THEN
     print *,'Particle BGC', call_count, nLevs
     CALL Particles(nLevs)
     DO col=1, nCols
       !IF (.NOT. active(col)) CYCLE  !# skip this column if dry
       CALL aed_calculate_particles(all_cols(:,col), col, nLevs)
     ENDDO
   ENDIF 

!-------------------------------------------------------------------------------
CONTAINS


   !############################################################################
   SUBROUTINE aed_run_column(icolm, col, nlev, doSurface)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: icolm(:)
      INTEGER,INTENT(in) :: col, nlev
      LOGICAL,INTENT(in) :: doSurface
   !
   !LOCALS
      INTEGER  :: v, lev, zon, split
   !
   !----------------------------------------------------------------------------
   !BEGIN
      DO split=1,split_factor
         IF (benthic_mode .GT. 1) THEN
            CALL p_calc_zone_areas(aedZones, aed_n_zones, data(col)%area, data(col)%lheights, nlev)
            CALL p_copy_to_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, nlev)
         ENDIF

         !# Update local light field (self-shading may have changed through
         !# changes in biological state variables). Update_light is set to
         !# be inline with current aed_phyoplankton, which requires only
         !# surface par, then integrates over depth of a layer
         CALL update_light(icolm, col, nlev)

         ! non PAR bandwidth fractions (set assuming single light extinction)
         data(col)%nir(:) = (data(col)%par(:)/par_fraction) * nir_fraction
         data(col)%uva(:) = (data(col)%par(:)/par_fraction) * uva_fraction
         data(col)%uvb(:) = (data(col)%par(:)/par_fraction) * uvb_fraction

         !# Time-integrate one biological time step
         CALL calculate_fluxes(icolm, col, nlev, doSurface)

         !# Update the water column layers
         DO v = 1, n_vars
            DO lev = 1, nlev
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
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, nlev)
         ELSE
            DO v = n_vars+1, n_vars+n_vars_ben
               data(col)%cc(v, bot) = data(col)%cc(v, bot) + dt_eff*flux_ben(v)
            ENDDO
         ENDIF

         CALL check_states(icolm, col, nlev)
      ENDDO
   END SUBROUTINE aed_run_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE pre_kinetics(icolm, col, nlev)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: icolm(:)
      INTEGER,INTENT(in) :: col, nlev
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
               IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
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
         DO lev = 1, nlev
            ! update ws for modules that use the mobility method
            !# direction doesn't seem to matter ? CAB
            CALL aed_mobility(icolm, lev, ws(:,lev))
         ENDDO

         !# Calculate source/sink terms due to the settling or rising of
         !# state variables in the water column (note that settling into benthos
         !# is done in aed_do_benthos)
         IF ( ASSOCIATED(doMobility) ) THEN
            v = 0
            DO i=1,n_aed_vars
               IF ( aed_get_var(i, tv) ) THEN
                  IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
                     v = v + 1
                     !# only for state_vars that are not sheet, and also non-zero ws
                     IF ( .NOT. isnan(tv%mobility) .AND. SUM(ABS(ws(i,:)))>zero_ ) THEN
                        min_C = tv%minimum
                        CALL doMobility(nlev, timestep, data(col)%dz, data(col)%area, ws(i,:), min_C, data(col)%cc(v, :))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      CALL check_states(icolm, col, nlev)

      IF (benthic_mode .GT. 1) THEN
         CALL p_calc_zone_areas(aedZones, aed_n_zones, data(col)%area, data(col)%lheights, nlev)
         CALL p_copy_to_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                                 data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, nlev)
      ENDIF

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer

      !# populate local light/extc arrays one column at a time
      IF (.NOT. link_ext_par) CALL update_light(icolm, col, nlev)

      ! non PAR bandwidth fractions (set assuming single light extinction)
      data(col)%nir = (data(col)%par/par_fraction) * nir_fraction
      data(col)%uva = (data(col)%par/par_fraction) * uva_fraction
      data(col)%uvb = (data(col)%par/par_fraction) * uvb_fraction
   END SUBROUTINE pre_kinetics
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE aed_initialize_zone_benthic(nCols, nlev, n_aed_vars, cc_diag)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in)   :: nCols, nlev
      INTEGER,INTENT(in)   :: n_aed_vars
      AED_REAL,INTENT(out) :: cc_diag(:,:)
   !
   !LOCALS
      INTEGER :: col, zon!, bot
   !
   !----------------------------------------------------------------------------
   !BEGIN
      DO zon=1, aed_n_zones
         aedZones(zon)%z_cc_diag(:, zon)  = zero_

         CALL aed_initialize_benthic(zon_cols(:,zon), 1)
      ENDDO

      CALL p_copy_from_zone(aedZones, aed_n_zones, data(col)%lheights, data(col)%cc,  &
                     data(col)%cc_hz, data(col)%cc_diag, data(col)%cc_diag_hz, nlev)
   END SUBROUTINE aed_initialize_zone_benthic
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE re_initialize(icolm, nlev)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: icolm(:)
      INTEGER,INTENT(in) :: nlev
   !
   !LOCALS
      INTEGER :: lev,zon,av,sv,sd
      TYPE(aed_variable_t),POINTER :: tvar
      TYPE(aed_column_t),DIMENSION(:),POINTER :: column_sed    !# (n_aed_vars)
   !
   !----------------------------------------------------------------------------
   !BEGIN
      DO lev=1, nlev
         CALL aed_initialize(icolm, lev)
      ENDDO

#if CUR_VARIANT == GLM_VARIANT
      !# (1) BENTHIC INITIALISATION
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers

         DO zon=1,aed_n_zones
            column_sed => zon_cols(:,zon)

        !   aedZones(zon)%z_env%z_sed_zones = zon  !MH TMP !CAB ???
!           print *,'aedZones(zon)%z_sed_zones',aedZones(zon)%z_sed_zones !MH TMP
            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
!           IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_env%z_sed_zones !CAB ???

            sv = 0 ; sd = 0

            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
                  IF ( tvar%sheet ) THEN
                     sd = sd + 1
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
                  ENDIF
               ELSEIF ( tvar%var_type == V_STATE ) THEN !# State variable
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
         CALL aed_initialize_benthic(icolm, 1)
      ELSE
         CALL aed_initialize_zone_benthic(nCols, nlev, n_aed_vars, data(col)%cc_diag)
      ENDIF
#endif
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE glm_benthics(icolm, col, wlev, bot)
   !----------------------------------------------------------------------------
   ! Calculate the benthic fluxes for GLM. This has been seperated from the
   ! general calc because its more complicated.
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: icolm(:)
      INTEGER, INTENT(in) :: col, wlev, bot
   !
   !LOCALS
      INTEGER :: lev,zon,v,v_start,v_end,av,sv,sd
      AED_REAL :: scale
      AED_REAL :: localrainl, localshade, localdrag
      LOGICAL :: splitZone
      TYPE(aed_variable_t),POINTER :: tvar
      TYPE(aed_column_t),DIMENSION(:),POINTER :: column_sed    !# (n_aed_vars)
      INTEGER :: layer_map(wlev), zlev
   !
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers, and output flux
         !# needs disaggregating from each zone back to the overlying cells/layers

!$OMP DO
         DO zon=1,aed_n_zones
            column_sed => zon_cols(:,zon)

            !# Reinitialise flux_ben to be repopulated for this zone
            flux_ben = zero_
            flux_pel_pre = zero_

            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => aedZones(zon)%z_env%z_sed_zones ! CAB???

            sv = 0 ; sd = 0
            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%sheet ) THEN
                  IF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
                     sd = sd + 1
                     column_sed(av)%cell_sheet => aedZones(zon)%z_cc_diag_hz(sd)
                  ELSEIF ( tvar%var_type == V_STATE ) THEN !# State variable
                     sv = sv + 1
                     IF ( tvar%bot ) THEN
                        column_sed(av)%cell_sheet => aedZones(zon)%z_cc(n_vars+sv, bot)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            !print*,"Calling ben for zone ",zon,zone_var,z_sed_zones(zon)

            !# (1) ZONE COLUMN UPDATING
            zlev = 0
            DO lev=1,wlev
               IF (aedZones(zon)%z_env%z_height < data(col)%lheights(lev)) THEN
                  zlev = lev
                  EXIT
               ENDIF
            ENDDO
            ! The upper height of the last zone may have no water above it.
            ! In this case, set the wet-layer vector to just be the top water layer
            IF(zlev == 0) zlev = wlev

            DO lev=zlev,wlev
              layer_map(lev) = zlev + wlev-lev
            ENDDO
            CALL aed_calculate_column(column_sed, layer_map)
       !# This is unused. Only macrophyte has a routine for this, but it's link is commented out

            IF ( benthic_mode .EQ. 3 ) THEN
               !# Zone is able to operated on by riparian and dry methods
               CALL aed_calculate_riparian(column_sed, zlev, aedZones(zon)%z_env%z_pc_wet) ! CAB???
               IF (aedZones(zon)%z_env%z_pc_wet < 0.01 ) CALL aed_calculate_dry(column_sed, zlev) !CAB ???

               !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
               CALL aed_rain_loss(icolm, bot, localrainl);
               IF (link_rain_loss) rain_factor = localrainl

               !# update feedback arrays to shade the water (ie reduce incoming light, Io)
               CALL aed_light_shading(icolm, bot, localshade)
               IF (link_solar_shade) sw_factor = localshade

               !# now the bgc updates are complete, update links to host model
               IF (aedZones(zon)%z_env%z_pc_wet > 0.99 ) THEN
                  CALL aed_bio_drag(icolm, bot, localdrag)
                  IF (link_bottom_drag) friction = localdrag
               ENDIF
            ENDIF

            !# Calculate temporal derivatives due to benthic processes.
            !# They are stored in flux_ben (benthic vars) and flux_pel (water vars)
            flux_pel_pre = flux_pel

            !print*,"Calling ben for zone ",zon,zone_var,aedZones(zon)%z_env%z_sed_zones
            !# pass zone instead of "bot" so aed puts it in a zone area
            CALL aed_calculate_benthic(column_sed, zon, .TRUE.)

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
                  splitZone = data(col)%lheights(lev-1) < aedZones(zon-1)%z_env%z_height
               ELSE
                  splitZone = 0.0 < aedZones(zon-1)%z_env%z_height ! CAB???
               ENDIF
            ELSE
               splitZone = .FALSE.
            ENDIF

            IF (splitZone) THEN
               IF (lev .GT. 1) THEN
                  scale = (aedZones(zon-1)%z_env%z_height - data(col)%lheights(lev-1)) / &
                                              (data(col)%lheights(lev) - data(col)%lheights(lev-1))
               ELSE
                  scale = (aedZones(zon-1)%z_env%z_height - 0.0) / (data(col)%lheights(lev) - 0.0)
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
            IF (lev > 1) flux_pel(:, lev) = flux_pel(:, lev) * (data(col)%area(lev)-data(col)%area(lev-1))/data(col)%area(lev)
            DO v=v_start,v_end
              IF ( data(col)%cc(v, 1) .GE. 0.0 ) flux_pel(v, lev) = &
                             max(-1.0 * data(col)%cc(v, lev), flux_pel(v, lev)/data(col)%dz(lev))
            END DO
         ENDDO
      ELSE
         !# Sediment zones are not simulated and therefore just operate on the bottom-most
         !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
         !# applied on flanks of the remaining layers, but note this is not suitable for
         !# model configurations where mass balance of benthic variables is required.

         !# Calculate temporal derivatives due to exchanges at the sediment/water interface
         IF ( zone_var .GE. 1 ) icolm(zone_var)%cell_sheet => aedZones(1)%z_env%z_sed_zones
         CALL aed_calculate_benthic(icolm, bot)

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
               CALL aed_calculate_benthic(icolm, lev)

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
   SUBROUTINE calculate_fluxes(icolm, col, nlev, doSurface)
   !----------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !----------------------------------------------------------------------------
   !ARGUMENTS
      TYPE(aed_column_t),INTENT(inout) :: icolm(:)
      INTEGER,INTENT(in) :: col, nlev
      LOGICAL,INTENT(in) :: doSurface
   !
   !LOCALS
      INTEGER :: lev
      INTEGER :: layer_map(nlev)
#if CUR_VARIANT == GLM_VARIANT
      LOGICAL :: glm_version = .TRUE.
#else
      LOGICAL :: glm_version = .FALSE.
#endif
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
      !# Now do the general calculation all flux terms for rhs in mass/m3/s
      !# Includes (i) benthic flux, (ii) surface exchange and (ii) kinetic updates in each cell
      !# as calculated by glm

      !# BENTHIC FLUXES
      IF ( glm_version ) THEN
         CALL glm_benthics(icolm, col, nlev, bot)
      ELSE
         IF ( do_zone_averaging ) THEN
            flux_pel(:,nlev) = flux_pel(:,nlev) + flux_pel_z(:, bot) !/h(nlev)

            !# Calculate temporal derivatives due to benthic exchange processes.
            CALL aed_calculate_benthic(icolm, bot, .FALSE.)
         ELSE
            CALL aed_calculate_benthic(icolm, bot)
         ENDIF

         !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
         flux_pel(:,bot) = flux_pel(:,bot)/data(col)%lheights(top)
      ENDIF

      !# SURFACE FLUXES
      !# Calculate temporal derivatives due to air-water exchange.
      IF (doSurface) THEN !# no surface exchange under ice cover
         CALL aed_calculate_surface(icolm, top)

         !# Distribute the fluxes into pelagic surface layer
         flux_pel(:, top) = flux_pel(:, top) + flux_atm(:)/data(col)%dz(top)
      ENDIF

      !# WATER CELL KINETICS
      !# Add pelagic sink and source terms in cells of all depth levels.
      DO lev=1,nlev
         CALL aed_calculate(icolm, lev)
      ENDDO
   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE aed_run_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(icolm, col, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: icolm(:)
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
   CALL aed_light_extinction(icolm, first, localext)

   ! Surface PAR
   data(col)%par(first) = par_fraction * data(col)%rad(first) * EXP( -(Kw+localext)*1e-6*data(col)%dz(first) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = start,end,step
      localext_up = localext
      CALL aed_light_extinction(icolm, i, localext)

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
         IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
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
