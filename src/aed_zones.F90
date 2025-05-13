!###############################################################################
!#                                                                             #
!# aed_zones.F90                                                               #
!#                                                                             #
!# The sediment zone processing bit for WQ                                     #
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

#include "aed_api.h"

MODULE aed_zones

   USE ISO_C_BINDING

   USE aed_common

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   !#-----------------------------------------------------------#!

   PUBLIC api_zone_t, aed_init_zones, api_set_zone_funcs
   PUBLIC aedZones, aed_n_zones
   PUBLIC p_calc_zone_areas, p_copy_to_zone, p_copy_from_zone
   PUBLIC calc_zone_areas_t, copy_to_zone_t, copy_from_zone_t

   !#===========================================================#!
   !# Structured type for Zones environment
   TYPE :: api_zone_env_t
      AED_REAL :: z_height
      AED_REAL :: z_dz
      AED_REAL :: z_area
      AED_REAL :: z_depth
      AED_REAL :: z_col_depth

      AED_REAL :: z_pc_wet
      AED_REAL :: z_col_num
      AED_REAL :: z_mat_id

      AED_REAL :: z_temp
      AED_REAL :: z_salt
      AED_REAL :: z_rho
      AED_REAL :: z_rad
      AED_REAL :: z_extc
      AED_REAL :: z_vel
      AED_REAL :: z_pres
      AED_REAL :: z_heatflux

      AED_REAL :: z_tss
      AED_REAL :: z_ss1
      AED_REAL :: z_ss2
      AED_REAL :: z_ss3
      AED_REAL :: z_ss4

      AED_REAL :: z_layer_stress
      AED_REAL :: z_sed_zone
      AED_REAL :: z_sed_zones
      AED_REAL :: z_wind
      AED_REAL :: z_air_temp
      AED_REAL :: z_air_pres
      AED_REAL :: z_rain
      AED_REAL :: z_evap
      AED_REAL :: z_humidity
      AED_REAL :: z_I_0
      AED_REAL :: z_longwave
      AED_REAL :: z_nir
      AED_REAL :: z_par
      AED_REAL :: z_uva
      AED_REAL :: z_uvb

      AED_REAL :: z_bathy
      AED_REAL :: z_biodrag
      AED_REAL :: z_bioextc
      AED_REAL :: z_solarshade
      AED_REAL :: z_windshade
      AED_REAL :: z_rainloss
   END TYPE api_zone_env_t
   !#===========================================================#!

   !#===========================================================#!
   !# Structured type for Zones
   TYPE :: api_zone_t
      INTEGER  :: n_levs

      TYPE(api_zone_env_t) :: z_env

      AED_REAL,DIMENSION(:,:),POINTER :: z_cc         !(n_vars, n_levs)
      AED_REAL,DIMENSION(:),  POINTER :: z_cc_hz      !(n_vars_ben, 2)
      AED_REAL,DIMENSION(:,:),POINTER :: z_cc_diag    !(n_diag_vars, n_levs)
      AED_REAL,DIMENSION(:),  POINTER :: z_cc_diag_hz !(n_diag_vars_hz, 2)

      AED_REAL,POINTER :: longitude   => null()
      AED_REAL,POINTER :: latitude    => null()
   END TYPE api_zone_t
   !#===========================================================#!

   INTEGER :: aed_n_zones
   TYPE(api_zone_t),DIMENSION(:),ALLOCATABLE,TARGET :: aedZones

   !----------------------------------------------------------------------------
   INTERFACE

     !#########################################################
     SUBROUTINE calc_zone_areas_t(theZones, n_zones, areas, heights, wlev)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: areas
        AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: heights
        INTEGER,INTENT(in) :: wlev
     END SUBROUTINE calc_zone_areas_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     !#########################################################
     SUBROUTINE copy_to_zone_t(theZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: heights
        AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_cc
        AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_cc_hz
        AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_diag
        AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_diag_hz
        INTEGER,INTENT(in) :: wlev
     END SUBROUTINE copy_to_zone_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     !#########################################################
     SUBROUTINE copy_from_zone_t(theZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: heights
        AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_cc
        AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_cc_hz
        AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_diag
        AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_diag_hz
        INTEGER,INTENT(in) :: wlev
     END SUBROUTINE copy_from_zone_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   END INTERFACE
   !----------------------------------------------------------------------------

   PROCEDURE(calc_zone_areas_t),POINTER :: p_calc_zone_areas
   PROCEDURE(copy_to_zone_t),   POINTER :: p_copy_to_zone
   PROCEDURE(copy_from_zone_t), POINTER :: p_copy_from_zone

CONTAINS

!###############################################################################
SUBROUTINE aed_init_zones(n_zones, n_levs, z_cc, z_cc_hz, z_diag, z_diag_hz)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: n_zones, n_levs
   AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_cc      !(n_vars, n_levs, n_zones)
   AED_REAL,DIMENSION(:,:),  POINTER,INTENT(in) :: z_cc_hz   !(n_vars, n_zones+1)
   AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_diag    !(n_vars, n_levs, n_zones)
   AED_REAL,DIMENSION(:,:)  ,POINTER,INTENT(in) :: z_diag_hz !(n_vars, n_zones+1)
!
!LOCALS
   INTEGER :: zon !, s, e
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_n_zones = n_zones
! print*,"z_cc(",size(z_cc, 1),",",size(z_cc,2),",",size(z_cc,3),")"

   ALLOCATE(aedZones(aed_n_zones+1))

   DO zon=1,n_zones
      aedZones(zon)%n_levs = n_levs

      aedZones(zon)%z_env%z_sed_zone = zon
      aedZones(zon)%z_env%z_sed_zones = zon

      aedZones(zon)%z_cc => z_cc(:, :, zon)
      aedZones(zon)%z_cc_hz => z_cc_hz(:, zon)
      aedZones(zon)%z_cc_diag => z_diag(:, :, zon)
      aedZones(zon)%z_cc_diag_hz => z_diag_hz(:, zon)
! print*,"aedZones%z_cc(",size(aedZones(zon)%z_cc, 1),",",size(aedZones(zon)%z_cc,2),")"
   ENDDO
END SUBROUTINE aed_init_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_zone_funcs(copy_to, copy_from, calc_areas)
!-------------------------------------------------------------------------------
!ARGUMENTS
   PROCEDURE(copy_to_zone_t),POINTER    :: copy_to
   PROCEDURE(copy_from_zone_t),POINTER  :: copy_from
   PROCEDURE(calc_zone_areas_t),POINTER :: calc_areas
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   p_calc_zone_areas => calc_areas
   p_copy_to_zone    => copy_to
   p_copy_from_zone  => copy_from
END SUBROUTINE api_set_zone_funcs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_zones
