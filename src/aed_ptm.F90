!###############################################################################
!#                                                                             #
!# aed_ptm.F90                                                                 #
!#                                                                             #
!# A generic interface between model and libaed-xxx for particle tracking      #
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

!#define CINTEGER INTEGER(kind=C_INT32_T)
!#define CSIZET   INTEGER(kind=C_SIZE_T)
!#define CLOGICAL LOGICAL(kind=C_BOOL)
!#define CCHARACTER CHARACTER(C_CHAR)

MODULE aed_ptm

   USE ISO_C_BINDING

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   PUBLIC aed_part_group_t

   !#--------------------------------------------------------------------------#
   !# Module Types

   !#===========================================================#!
   TYPE :: aed_part_group_t
      CINTEGER :: num_particles                      !# Number of particles
      CINTEGER :: idx_stat, idx_2, idx_3, idx_layer  !# Particle ISTAT Index Values
      CINTEGER :: idx_bed_layer, idx_motility        !# Particle ISTAT Index Values
      CINTEGER :: idx_uvw0, idx_uvw, idx_nu, idx_wnd !# Particle PROP Index Values
      CINTEGER :: idx_wsel, idx_watd, idx_partd      !# Particle PROP Index Values
      CINTEGER :: idx_age, idx_state                 !# Particle TSTAT Index Values
      CINTEGER :: next                               !# next particle index
      CINTEGER :: stat                               !# particle status
      AED_REAL,DIMENSION(:,:),POINTER :: age         !# particle time/age vector (2,Npart)
      AED_REAL,DIMENSION(:,:),POINTER :: posn        !# particle position vector
      AED_REAL,DIMENSION(:,:),POINTER :: prop        !# particle property vector (12,Npart)
      AED_REAL,DIMENSION(:,:),POINTER :: U           !# particle conserved variable vector (NU,NP)
   ENDTYPE aed_part_group_t

!  TYPE :: partgroup_p
!     INTEGER :: idx, grp
!  ENDTYPE
!  TYPE :: partgroup_cell
!      INTEGER :: count, n
!      TYPE(partgroup_p),ALLOCATABLE,DIMENSION(:) :: prt
!  END TYPE partgroup_cell
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   INTEGER :: aed_n_groups


!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed_ptm_init()
!-------------------------------------------------------------------------------
! Initialise the particle tracker.
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
END SUBROUTINE aed_ptm_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_ptm
