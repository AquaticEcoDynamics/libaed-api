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
!# Copyright 2024-2025 - The University of Western Australia                   #
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

MODULE aed_ptm

   USE ISO_C_BINDING
   USE aed_common
   USE aed_core, ONLY:  aed_ptm_t

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   PUBLIC aed_part_group_t, aed_ptm_init, ptm_istat, ptm_env, Particles, aed_calculate_particles, set_ptm_aed_var_num

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

   INTEGER, PARAMETER :: n_ptm_istat = 6
   INTEGER, PARAMETER :: n_ptm_env   = 5
   INTEGER            :: n_ptm_vars  = 0
   INTEGER            :: aed_n_particles
   INTEGER            :: n_aed_vars_

   !# Particle groups
!   INTEGER :: num_groups
!   TYPE(partgroup),DIMENSION(:),POINTER :: particle_groups
   TYPE(partgroup_cell),DIMENSION(:),ALLOCATABLE, TARGET :: all_particles

   INTEGER, PARAMETER :: STAT = 1
   INTEGER, PARAMETER :: IDX2 = 2
   INTEGER, PARAMETER :: IDX3 = 3
   INTEGER, PARAMETER :: LAYR = 4
   INTEGER, PARAMETER :: FLAG = 5
   INTEGER, PARAMETER :: PTID = 6
!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_ptm_init(ng,np,parts,n_ptm_vars_,n_cells)
!-------------------------------------------------------------------------------
! Initialise the particle tracker.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER, INTENT(in) :: ng,np
   INTEGER,  INTENT(in) :: n_ptm_vars_, n_cells
   TYPE(aed_part_group_t),DIMENSION(:),TARGET,INTENT(in) :: parts
!LOCALS
   TYPE(aed_variable_t),POINTER :: tvar
   TYPE(aed_ptm_t), DIMENSION(:), ALLOCATABLE :: ptm
   INTEGER :: rc, pv, av, grp, prt, ppid
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_n_groups = ng
   aed_n_particles = np

   n_ptm_vars = n_ptm_vars_

   IF (.NOT. ASSOCIATED(particle_groups)) THEN
      particle_groups => parts
   ENDIF

   !# Allocating AED PTM arrays : status, environment, state and diagnostic
   ALLOCATE(ptm_istat(aed_n_groups,1:aed_n_particles,1:n_ptm_istat),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_istat)'
   ALLOCATE(ptm_env(aed_n_groups,1:aed_n_particles,1:n_ptm_env+n_ptm_vars),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_env)'
   ALLOCATE(ptm_state(aed_n_groups,1:aed_n_particles,1:n_ptm_vars),stat=rc)
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_state)'
   ALLOCATE(ptm_diag(aed_n_groups,1:aed_n_particles,1:n_ptm_vars),stat=rc) ! Not yet used
     IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ptm_diag)'

   !#---------------------------------------------------------------------------
   !# Now set initial values
   ptm_istat(:,:,:) = -9999
   ptm_env(:,:,:)   = -9999.
   ptm_state(:,:,:) = -9999.
   ptm_diag(:,:,:)  = -9999.

   pv = 0;
   DO av=1,n_aed_vars_
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "     ERROR getting variable info"
      IF ( tvar%var_type == V_PARTICLE ) THEN  !# ptm variable
          pv = pv + 1

         !print *,'PTM',pv,n_ptm_env,tvar%initial
          ptm_state(:,:,pv) = tvar%initial !# Note this is all particles, regardless of status
                                           !#    (ptm_env(:,:,n_ptm_env+1:n_ptm_env+n_ptm_vars))
          ptm_env(:,:,n_ptm_env+pv) = tvar%initial
      ENDIF
   ENDDO

   ALLOCATE(ptm(aed_n_particles*aed_n_groups))

   DO grp=1,aed_n_groups
     DO prt=1,aed_n_particles
         ! Point single particle object to the global particle data structure
         ptm(prt)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(prt)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(prt)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(prt)%ptm_diag  => ptm_diag(grp,prt,:)
     ENDDO
   ENDDO !end particle loop

   ppid = aed_n_particles*aed_n_groups

   CALL aed_initialize_particle(ppid,ptm)
   DEALLOCATE(ptm)

   ! !TESTS
   ! ptm_istat(1,5000,1) = 42
   ! ptm_istat(1,1,1) = 101
   ! ptm_istat(1,1,2) = 102
   ! ptm_istat(1,1,3) = 103
   ! ptm_istat(1,1,4) = 104
   !
   ! ptm_istat(1,2,1) = 201
   ! ptm_istat(1,2,2) = 202
   ! ptm_istat(1,2,3) = 203
   ! ptm_istat(1,2,4) = 204
   !
   ! ptm_env(1,3,3) = 1.33
   ! ptm_env(1,1,1) = 111.1
   !
   !
   !!print*,"allocating all_parts with ", ubound(temp,1), " cells"
   ALLOCATE(all_particles(n_cells))
END SUBROUTINE aed_ptm_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Particles(n_cells)
!-------------------------------------------------------------------------------
! Calculate biogeochemical transformations on particles
!-------------------------------------------------------------------------------
!ARGUMENTS
!   TYPE (aed_column_t), INTENT(inout) :: column(:)
   !TYPE(partgroup_cell), INTENT(inout) :: layer_particles(:)
   INTEGER,  INTENT(in)    :: n_cells
!
!LOCAL VARIABLES:
   INTEGER :: lev, grp, prt, n, pt, NU
   INTEGER :: ppid
   AED_REAL,DIMENSION(20) :: zz
   INTEGER :: cell, j, ii, count_check
   AED_REAL :: dt = 3600

   TYPE (aed_ptm_t) :: ptm
!
!-------------------------------------------------------------------------------

!BEGIN
   IF (aed_n_groups == 0 .OR. aed_n_particles == 0) RETURN
   zz = zero_

!------------
   !print*,"PTM START", aed_n_groups, aed_n_particles

   DO cell=1, size(all_particles)
      IF (ALLOCATED(all_particles(cell)%prt)) DEALLOCATE(all_particles(cell)%prt)
      all_particles(cell)%count = 0
   ENDDO
   count_check = 0
   DO grp=1,aed_n_groups
      !print*,"PTM GRP", grp

      ! First, loop through all particles, and count how mnay are in each cell
      DO prt=1,aed_n_particles
        !IF (prt<30) print *,'STAT', prt,ptm_istat(grp,prt,STAT),ptm_istat(grp,prt,IDX3)
         IF ( ptm_istat(grp,prt,STAT) >= 0 ) THEN
            cell = ptm_istat(grp,prt,IDX3)
            IF ( cell >= 1 .AND. cell <= size(all_particles) ) THEN
               all_particles(cell)%count = all_particles(cell)%count + 1
           !ELSE
           !   print*,"idx out of range", i, size(all_particles)
           !   stop
            ENDIF
         ENDIF
      ENDDO
      DO ii=1,n_cells
          count_check = count_check+all_particles(ii)%count
          !print*,"PTM CELL", ii, all_particles(ii)%count
      ENDDO
      !print*,"PTM CHK", count_check

!  ENDDO
!  DO grp=1,num_groups
      ! Now, loop through all particles, and populate the particle-cell object
      DO prt=1,aed_n_particles
         IF ( ptm_istat(grp,prt,STAT) < 0 ) CYCLE  !# ignore these

         cell = ptm_istat(grp,prt,IDX3)
         IF ( cell >= 1 .AND. cell <= size(all_particles) ) THEN
            ! If the particle is in a cell, first check if the particle-cell
            ! object is already allocated. If not, allocate it and init n to 0
            IF (.NOT. ALLOCATED(all_particles(cell)%prt)) THEN
               ALLOCATE(all_particles(cell)%prt(all_particles(cell)%count))
               all_particles(cell)%n = 0
            ENDIF
            ! Increment the new particle that was found in the particle-cell object
            ! and if this is less than the total count, add the new particle to the list
            j = all_particles(cell)%n + 1              ! add new particle
            IF (j <= all_particles(cell)%count ) THEN
               all_particles(cell)%prt(j)%grp = grp
               all_particles(cell)%prt(j)%idx = prt
               all_particles(cell)%n = j
               !print*,"PTM", grp, prt, cell, j, all_particles(cell)%count
           !ELSE
           !   print*,"Ooops, error in PTM", j, all_particles(i)%count
            ENDIF
!        ELSE
!           print*,"idx out of range", i, size(all_particles)
!           print*,"grp", grp, " prt ",prt
!           print*,"istat 1", particle_groups(grp)%istat(1,prt)
!           print*,"istat 2", particle_groups(grp)%istat(2,prt)
!           print*,"istat 3", particle_groups(grp)%istat(3,prt)
!           print*,"istat 4", particle_groups(grp)%istat(4,prt)
!           stop
         ENDIF
      ENDDO ! particles
   ENDDO    ! groups
!-------

END SUBROUTINE Particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_particles(icolm, col, nlev)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_column_t),INTENT(inout) :: icolm(:)
   INTEGER,INTENT(in)               :: col, nlev
!
!LOCAL VARIABLES:
   INTEGER :: lev, grp, prt, n, pt, NU
   INTEGER :: ppid, i, pid, new_prt
   AED_REAL :: dt = 3600
   AED_REAL :: div
   LOGICAL :: pass

   TYPE (aed_ptm_t), DIMENSION(:), ALLOCATABLE :: ptm

   TYPE(partgroup_cell), POINTER :: layer_particles
!
!-------------------------------------------------------------------------------

!BEGIN
   IF (aed_n_groups == 0 .OR. aed_n_particles == 0) RETURN

!  DO lev=1,nlev
!  ENDDO

   DO lev=1,nlev
      layer_particles => all_particles(lev)

      !print *, "ptm", lev, layer_particles%count
      IF (layer_particles%count == 0) CYCLE

      ALLOCATE(ptm(layer_particles%count))

      ppid = 0          ! new cell identifier, to allow cumulation of prts
      DO pt=1,layer_particles%count

         ! Retrieve particle properties, from the particle-cell object
         grp = layer_particles%prt(pt)%grp ; prt = layer_particles%prt(pt)%idx

         !print *, "ppp", lev, pt, grp, prt

         ! Point single particle object to the global particle data structure
         ptm(pt)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(pt)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(pt)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(pt)%ptm_diag  => ptm_diag(grp,prt,:)

         !print *,'ptm_istat(grp,prt,STAT)',ptm_istat(grp,prt,STAT), ptm%ptm_istat
      ENDDO !end particle loop
      ppid = layer_particles%count

      ! Pass through the particle to AED modules, if its active
      !IF ( ptm_istat(grp,prt,STAT) >= 0 ) THEN
         CALL aed_particle_bgc(icolm,lev,ppid,p=ptm) ! Note: ppid getting incremeted in here
      !ENDIF
      DEALLOCATE(ptm)
   ENDDO !end layer loop

   ALLOCATE(ptm(aed_n_particles*aed_n_groups))

   DO grp=1,aed_n_groups
     DO prt=1,aed_n_particles
         ! Point single particle object to the global particle data structure
         ptm(prt)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(prt)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(prt)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(prt)%ptm_diag  => ptm_diag(grp,prt,:)
     ENDDO
   ENDDO !end particle loop

   ppid = aed_n_particles*aed_n_groups

   CALL aed_split_particle(ppid,ptm)
   DEALLOCATE(ptm)

END SUBROUTINE aed_calculate_particles
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

        !   CALL aed_particle_bgc(column,lev,ppid,zz)     !ppid getting incremeted in here

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
            particle_groups(grp)%status(stat, prt) = INT(zz(19))
         ENDIF
         particle_groups(grp)%age(2,prt) = particle_groups(grp)%age(2,prt) + dt
      ENDDO
   ENDDO
END SUBROUTINE Particles_zz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE set_ptm_aed_var_num(n_aed_vars)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: n_aed_vars
!
!LOCALS
  !INTEGER n_aed_vars_
!
!-------------------------------------------------------------------------------
!BEGIN

  n_aed_vars_ = n_aed_vars
END SUBROUTINE set_ptm_aed_var_num
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_ptm
