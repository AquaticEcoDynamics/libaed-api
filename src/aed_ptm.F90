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

!#define CINTEGER INTEGER(kind=C_INT32_T)
!#define CSIZET   INTEGER(kind=C_SIZE_T)
!#define CLOGICAL LOGICAL(kind=C_BOOL)
!#define CCHARACTER CHARACTER(C_CHAR)

MODULE aed_ptm

   USE ISO_C_BINDING
   USE aed_common
   USE aed_core, ONLY:  aed_ptm_t

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   PUBLIC aed_part_group_t, aed_ptm_init, ptm_istat, ptm_env


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
      CINTEGER,DIMENSION(:,:),POINTER :: status      !# particle status (4,Npart)
      AED_REAL,DIMENSION(:,:),POINTER :: age         !# particle time/age vector (2,Npart)
      AED_REAL,DIMENSION(:,:),POINTER :: posn        !# particle position vector
      AED_REAL,DIMENSION(:,:),POINTER :: prop        !# particle property vector (12,Npart)
      AED_REAL,DIMENSION(:,:),POINTER :: vars        !# particle conserved variable vector (NU,NPart)
   ENDTYPE aed_part_group_t

   TYPE :: partgroup_p
      INTEGER :: idx, grp
   ENDTYPE
   TYPE :: partgroup_cell
      INTEGER :: count, n
      TYPE(partgroup_p),ALLOCATABLE,DIMENSION(:) :: prt
   END TYPE partgroup_cell


!  TYPE :: partgroup_cell
!      INTEGER :: count, n
!      TYPE(partgroup_p),ALLOCATABLE,DIMENSION(:) :: prt
!  END TYPE partgroup_cell


!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   INTEGER :: aed_n_groups
   TYPE(aed_part_group_t),DIMENSION(:),POINTER :: particle_groups

   CINTEGER,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_istat       !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_env         !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_state       !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_diag        !# AED particle data structure (NGroups,NParticles,NAttributes)

   INTEGER, PARAMETER :: n_ptm_istat = 5
   INTEGER, PARAMETER :: n_ptm_env = 10
   INTEGER :: aed_n_particles


!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed_ptm_init(ng,np,parts,  n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)
!-------------------------------------------------------------------------------
! Initialise the particle tracker.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER, INTENT(in) :: ng,np
   INTEGER, INTENT(in) :: n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
   TYPE(aed_part_group_t),DIMENSION(:),TARGET,INTENT(in) :: parts
!LOCALS
   INTEGER :: rc
!
!-------------------------------------------------------------------------------
!BEGIN

   print *,'sss',ng,np
   aed_n_groups = ng
   aed_n_particles = np

   IF (.NOT. ASSOCIATED(particle_groups)) THEN
      particle_groups => parts
   ENDIF

   ! Allocating AED PTM arrays : status, environment, state and diagnostic
   ALLOCATE(ptm_istat(aed_n_groups,1:aed_n_particles,1:n_ptm_istat),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_istat)'
   ALLOCATE(ptm_env(aed_n_groups,1:aed_n_particles,1:n_ptm_env),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_env)'
   ALLOCATE(ptm_state(aed_n_groups,1:aed_n_particles,1:(n_vars+n_vars_ben)),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_state)'
   ALLOCATE(ptm_diag(aed_n_groups,1:aed_n_particles,1:(n_vars_diag+n_vars_diag_sheet)),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_diag)'

   ptm_istat(:,:,:) = -9999
   ptm_env(:,:,:)   = -9999.
   ptm_state(:,:,:) = -9999.
   ptm_diag(:,:,:)  = -9999.

    !TESTS
    ptm_istat(1,5000,1) = 42
    ptm_istat(1,1,1) = 101
    ptm_istat(1,1,2) = 102
    ptm_istat(1,1,3) = 103
    ptm_istat(1,1,4) = 104

    ptm_istat(1,2,1) = 201
    ptm_istat(1,2,2) = 202
    ptm_istat(1,2,3) = 203
    ptm_istat(1,2,4) = 204

    ptm_env(1,3,3) = 1.33
    ptm_env(1,1,1) = 111.1


END SUBROUTINE aed_ptm_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Particles(column, count, layer_particles)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles 
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   TYPE(partgroup_cell), INTENT(inout) :: layer_particles(:)
   INTEGER,  INTENT(in)    :: count
!
!LOCAL VARIABLES:
   INTEGER :: lev, grp, prt, n, pt, NU
   INTEGER :: ppid
   AED_REAL,DIMENSION(20) :: zz
   INTEGER :: stat, idxi3
   AED_REAL :: dt = 3600

   TYPE (aed_ptm_t) :: ptm

!
!-------------------------------------------------------------------------------
 
!BEGIN
   IF (aed_n_groups == 0 .OR. aed_n_particles == 0) RETURN
   zz = zero_

   DO lev=1,count

      ppid = 0          ! new cell identifier, to allow cumulation of prts

      DO pt=1,layer_particles(lev)%count

         grp = layer_particles(lev)%prt(pt)%grp ; prt = layer_particles(lev)%prt(pt)%idx
         !stat  = ptm_istat(grp,part,idx_stat) ! particle_groups(grp)%idx_stat   ! should be 1
         !idxi3 = ptm_istat(grp,part,idx_3)    ! should be 3

         ptm%ptm_istat => ptm_istat(grp,prt,:)
         ptm%ptm_env => ptm_env(grp,prt,:)
         ptm%ptm_state => ptm_state(grp,prt,:)
         ptm%ptm_diag => ptm_diag(grp,prt,:)

         IF ( ptm_istat(grp,prt,1) >= 0 ) THEN

            !CALL aed_particle_bgc(column,lev,ppid,zz)     !ppid getting incremeted in here

            CALL aed_particle_bgc(column,lev,ppid,p=ptm)     

         ENDIF

      ENDDO
   ENDDO
END SUBROUTINE Particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE Particles_zz(column, count, parts)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles 
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   TYPE(partgroup_cell), INTENT(inout) :: parts(:)
   INTEGER,  INTENT(in)    :: count
!
!LOCAL VARIABLES:
   INTEGER :: lev, grp, prt, n, pt, NU
   INTEGER :: ppid
   AED_REAL,DIMENSION(20) :: zz
   INTEGER :: stat, idxi3
   AED_REAL :: dt = 3600

!
!-------------------------------------------------------------------------------
 
!BEGIN
   IF (.NOT. ASSOCIATED(particle_groups) .OR. aed_n_groups == 0) RETURN
   zz = zero_

   DO lev=1,count

      ppid = 0          ! new cell identifier, to allow cumulation of prts

      DO pt=1,parts(lev)%count

         grp = parts(lev)%prt(pt)%grp ; prt = parts(lev)%prt(pt)%idx
         stat = particle_groups(grp)%idx_stat   ! should be 1
         idxi3 =  particle_groups(grp)%idx_3   ! should be 3

         IF ( particle_groups(grp)%status(stat, prt) >= 0 ) THEN


            NU = ubound(particle_groups(grp)%vars, 1)
            n = min(16, size(particle_groups(grp)%prop(:,prt)))

          ! zz(1:n) = particle_groups(grp)%prop(1:n,prt)
            zz(1)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0, prt)
            zz(2)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0+1, prt)
            zz(3)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0+2, prt)
            zz(4)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw, prt)
            zz(5)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw+1, prt)
            zz(6)  = particle_groups(grp)%prop(particle_groups(grp)%idx_uvw+2, prt)
            zz(7)  = particle_groups(grp)%prop(particle_groups(grp)%idx_nu, prt)
            zz(8)  = particle_groups(grp)%prop(particle_groups(grp)%idx_nu+1, prt)
            zz(9)  = particle_groups(grp)%prop(particle_groups(grp)%idx_nu+2, prt)
            zz(10) = particle_groups(grp)%prop(particle_groups(grp)%idx_nu+3, prt)
            zz(11) = particle_groups(grp)%prop(particle_groups(grp)%idx_wsel, prt)
            zz(12) = particle_groups(grp)%prop(particle_groups(grp)%idx_watd, prt)
            zz(13) = particle_groups(grp)%prop(particle_groups(grp)%idx_partd, prt)
            zz(14) = particle_groups(grp)%prop(particle_groups(grp)%idx_wnd, prt) !Vvel

            IF (NU > 0) zz(15) = particle_groups(grp)%vars(1, prt)  !Mass
            IF (NU > 1) zz(16) = particle_groups(grp)%vars(2, prt)

            zz(17:18) = particle_groups(grp)%age(1:2,prt)   !Birth and Age
            zz(19) = particle_groups(grp)%status(stat, prt)    !Status

        !    CALL aed_particle_bgc(column,lev,ppid,zz)     !ppid getting incremeted in here

           !particle_groups(grp)%prop(1:n,prt) = zz(1:n)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0, prt)   = zz(1)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0+1, prt) = zz(2)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw0+2, prt) = zz(3)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw, prt)    = zz(4)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw+1, prt)  = zz(5)
            particle_groups(grp)%prop(particle_groups(grp)%idx_uvw+2, prt)  = zz(6)
            particle_groups(grp)%prop(particle_groups(grp)%idx_nu, prt)     = zz(7)
            particle_groups(grp)%prop(particle_groups(grp)%idx_nu+1, prt)   = zz(8)
            particle_groups(grp)%prop(particle_groups(grp)%idx_nu+2, prt)   = zz(9)
            particle_groups(grp)%prop(particle_groups(grp)%idx_nu+3, prt)   = zz(10)
            particle_groups(grp)%prop(particle_groups(grp)%idx_wsel, prt)   = zz(11)
            particle_groups(grp)%prop(particle_groups(grp)%idx_watd, prt)   = zz(12)
            particle_groups(grp)%prop(particle_groups(grp)%idx_partd, prt)  = zz(13)
            particle_groups(grp)%prop(particle_groups(grp)%idx_wnd, prt)    = zz(14)

            IF (NU > 0) particle_groups(grp)%vars(1, prt) = zz(15)
            IF (NU > 1) particle_groups(grp)%vars(2, prt) = zz(16)
            particle_groups(grp)%status(stat, prt) = zz(19)
         ENDIF
         particle_groups(grp)%age(2,prt) = particle_groups(grp)%age(2,prt) + dt
      ENDDO
   ENDDO
END SUBROUTINE Particles_zz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_ptm
