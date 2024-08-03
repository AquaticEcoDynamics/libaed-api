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

MODULE aed_zones

   USE ISO_C_BINDING

   USE aed_common

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   !#===========================================================#!
   !# Structured type for Zones
   TYPE :: api_zone_t
      INTEGER :: n_levs
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_temp
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_salt
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_heights
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_rad
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_rho
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_area
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_extc
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_layer_stress
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_tss
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_dz
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_vel
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_par
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_nir
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_uva
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_uvb
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_pres
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_depth
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_sed_zones
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_pc_wet
      AED_REAL,DIMENSION(:),ALLOCATABLE :: z_heatflux

      AED_REAL,DIMENSION(:,:),POINTER :: z_cc         !(n_levs, n_vars)
      AED_REAL,DIMENSION(:),  POINTER :: z_cc_hz      !(2, n_vars_ben)
      AED_REAL,DIMENSION(:,:),POINTER :: z_cc_diag    !(n_levs, n_diag_vars)
      AED_REAL,DIMENSION(:),  POINTER :: z_cc_diag_hz !(2, n_diag_vars_hz)
   END TYPE api_zone_t
   !#===========================================================#!

   INTEGER :: aed_n_zones
   TYPE(api_zone_t),DIMENSION(:),ALLOCATABLE,TARGET :: aedZones

   !----------------------------------------------------------------------------
   INTERFACE

     !#########################################################
     SUBROUTINE calc_zone_areas_t(theZones, n_zones, areas, wlev, surf)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:),INTENT(in) :: areas
        INTEGER,INTENT(in) :: wlev
        AED_REAL :: surf
     END SUBROUTINE calc_zone_areas_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     !#########################################################
     SUBROUTINE copy_to_zone_t(theZones, n_zones, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:,:),INTENT(in) :: x_cc
        AED_REAL,DIMENSION(:),INTENT(in) :: x_cc_hz
        AED_REAL,DIMENSION(:,:),INTENT(in) :: x_diag
        AED_REAL,DIMENSION(:),INTENT(in) :: x_diag_hz
        INTEGER,INTENT(in) :: wlev
     END SUBROUTINE copy_to_zone_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     !#########################################################
     SUBROUTINE copy_from_zone_t(theZones, n_zones, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
        IMPORT :: api_zone_t
        TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: theZones
        INTEGER,INTENT(in) :: n_zones
        AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_cc
        AED_REAL,DIMENSION(:),INTENT(inout) :: x_cc_hz
        AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_diag
        AED_REAL,DIMENSION(:),INTENT(inout) :: x_diag_hz
        INTEGER,INTENT(in) :: wlev
     END SUBROUTINE copy_from_zone_t
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   END INTERFACE
   !----------------------------------------------------------------------------

   PROCEDURE(calc_zone_areas_t),POINTER :: p_calc_zone_areas
   PROCEDURE(copy_to_zone_t),POINTER    :: p_copy_to_zone
   PROCEDURE(copy_from_zone_t),POINTER  :: p_copy_from_zone

   PUBLIC api_zone_t, aed_init_zones, api_set_zone_funcs
   PUBLIC aedZones, aed_n_zones
   PUBLIC p_calc_zone_areas, p_copy_to_zone, p_copy_from_zone
   PUBLIC calc_zone_areas_t, copy_to_zone_t, copy_from_zone_t

CONTAINS

!###############################################################################
SUBROUTINE aed_init_zones(n_zones, n_levs, z_cc, z_cc_hz, z_diag, z_diag_hz)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: n_zones, n_levs
   AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_cc      !(n_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:),  POINTER,INTENT(in) :: z_cc_hz   !(n_zones+1, n_vars)
   AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_diag    !(n_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:)  ,POINTER,INTENT(in) :: z_diag_hz !(n_zones+1, n_vars)
!
!LOCALS
   INTEGER :: zon !, s, e
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_n_zones = n_zones
! print*,"z_cc(",size(z_cc, 1),",",size(z_cc,2),",",size(z_cc,3),")"

   ALLOCATE(aedZones(aed_n_zones))

   DO zon=1,n_zones
      aedZones(zon)%n_levs = n_levs
      ALLOCATE(aedZones(zon)%z_temp(n_levs))
      ALLOCATE(aedZones(zon)%z_salt(n_levs))
      ALLOCATE(aedZones(zon)%z_heights(n_levs))
      ALLOCATE(aedZones(zon)%z_rad(n_levs))
      ALLOCATE(aedZones(zon)%z_rho(n_levs))
      ALLOCATE(aedZones(zon)%z_area(n_levs))
      ALLOCATE(aedZones(zon)%z_extc(n_levs))
      ALLOCATE(aedZones(zon)%z_layer_stress(n_levs))
      ALLOCATE(aedZones(zon)%z_tss(n_levs))
      ALLOCATE(aedZones(zon)%z_dz(n_levs))
      ALLOCATE(aedZones(zon)%z_vel(n_levs))
      ALLOCATE(aedZones(zon)%z_par(n_levs))
      ALLOCATE(aedZones(zon)%z_nir(n_levs))
      ALLOCATE(aedZones(zon)%z_uva(n_levs))
      ALLOCATE(aedZones(zon)%z_uvb(n_levs))
      ALLOCATE(aedZones(zon)%z_pres(n_levs))
      ALLOCATE(aedZones(zon)%z_depth(n_levs))
      ALLOCATE(aedZones(zon)%z_sed_zones(n_levs))
      ALLOCATE(aedZones(zon)%z_pc_wet(n_levs))
      ALLOCATE(aedZones(zon)%z_heatflux(n_levs))

      aedZones(zon)%z_cc => z_cc(zon, :, :)
      aedZones(zon)%z_cc_hz => z_cc_hz(zon, :)
      aedZones(zon)%z_cc_diag => z_diag(zon, :, :)
      aedZones(zon)%z_cc_diag_hz => z_diag_hz(zon, :)
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
