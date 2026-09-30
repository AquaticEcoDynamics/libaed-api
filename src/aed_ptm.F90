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
   USE aed_core, ONLY:  aed_ptm_t, aed_ptm_dt_days, aed_ptm_dt_days_set

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   PUBLIC aed_part_group_t, aed_ptm_init, ptm_istat, ptm_env, Particles, aed_calculate_particles, &
          aed_split_particles, aed_ptm_set_cell_map, set_ptm_aed_var_num, &
          Particles_impl, aed_calculate_particles_impl, aed_split_particles_impl, &
          set_ptm_aed_var_num_impl   !# S5: the unguarded bodies, for aed_api's internal calls only

   !# S5 2026-09-15: detection-only overlap guard for the public entry points (libaed-api/src/aed_api_guard.c)
   INTERFACE
      SUBROUTINE aed_api_guard_enter(id) BIND(C, name="aed_api_guard_enter")
         USE ISO_C_BINDING, ONLY : C_INT
         INTEGER(C_INT), VALUE :: id
      END SUBROUTINE aed_api_guard_enter
      SUBROUTINE aed_api_guard_leave(id) BIND(C, name="aed_api_guard_leave")
         USE ISO_C_BINDING, ONLY : C_INT
         INTEGER(C_INT), VALUE :: id
      END SUBROUTINE aed_api_guard_leave
   END INTERFACE

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

   !# 2026-09-09: initialised. These are assigned ONLY in aed_ptm_init, and the guard that
   !# protects every entry point - `IF (aed_n_groups == 0 .OR. aed_n_particles == 0 .OR.
   !# n_ptm_vars == 0) RETURN` - READS them. On the one path where aed_ptm_init is not reached
   !# (glm_ptm.c: .NOT.do_particle_bgc .AND. Num_PTM_Vars == 0) that read was of an
   !# uninitialised variable: undefined behaviour that happened to work only because gfortran
   !# places module data in .bss and zeroes it. Nothing currently reaches those entry points on
   !# that path, so no behaviour changes - but the guard is now well-defined rather than
   !# relying on the storage model.
   INTEGER :: aed_n_groups = 0
   TYPE(aed_part_group_t),DIMENSION(:),POINTER :: particle_groups

   CINTEGER,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_istat       !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_env         !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_state       !# AED particle data structure (NGroups,NParticles,NAttributes)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE, TARGET :: ptm_diag        !# AED particle data structure (NGroups,NParticles,NAttributes)

   ! Was 6 (STAT,IDX2,IDX3,LAYR,FLAG,PTID). Bumped to 7 for aed_phyto_abm's GRP
   ! attribute (which phyto species/group a particle belongs to, added to support
   ! multiple phytoplankton groups in the ABM) - GRP=7 is written/read via
   ! p(i)%ptm_istat(GRP) in aed_phyto_abm.F90, so the underlying array must be at
   ! least that large or it's an out-of-bounds write into whatever memory follows.
   INTEGER, PARAMETER :: n_ptm_istat = 7
   INTEGER, PARAMETER :: n_ptm_env   = 5
   INTEGER            :: n_ptm_vars  = 0
   INTEGER            :: aed_n_particles = 0
   INTEGER            :: n_aed_vars_

   !# Particle groups
!   INTEGER :: num_groups
!   TYPE(partgroup),DIMENSION(:),POINTER :: particle_groups
   TYPE(partgroup_cell),DIMENSION(:),ALLOCATABLE, TARGET :: all_particles

   !# Optional host map (layer,col) -> global cell index, used to translate a
   !# column's local levels into the global cell space the particles are binned
   !# in (IDX3). For single-column hosts (GLM) this is left unset and the global
   !# cell defaults to the (layer) index - preserving existing behaviour.
   INTEGER,DIMENSION(:,:),ALLOCATABLE :: ptm_cell_id
   LOGICAL :: have_cell_id = .FALSE.

   INTEGER, PARAMETER :: STAT = 1
   INTEGER, PARAMETER :: IDX2 = 2
   INTEGER, PARAMETER :: IDX3 = 3
   INTEGER, PARAMETER :: LAYR = 4
   INTEGER, PARAMETER :: FLAG = 5
   INTEGER, PARAMETER :: PTID = 6
!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE ptm_alloc_check(rc, n, what)
!-------------------------------------------------------------------------------
! Fatal report for a failed or implausible ALLOCATE.
!
! Several ALLOCATEs in this module run EVERY TIMESTEP and carried no stat=, so a
! failure was undefined behaviour rather than a reported error: two runs were
! silently truncated (45 of 78 steps) with the log still saying they completed.
!
! stat= alone is not sufficient. A size computed from an overflowed default
! INTEGER product can come out small or negative, and ALLOCATE of a negative
! extent SUCCEEDS in Fortran - it yields a zero-size array - so the run would
! continue on an array that cannot hold the particles it is about to be indexed
! with. The size is therefore validated separately from the status.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: rc, n
   CHARACTER(len=*), INTENT(in) :: what
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (n < 0) THEN
      WRITE(*,*) 'aed_ptm: ERROR allocating ', TRIM(what), &
                 ' : negative element count (', n, ') - integer overflow in the size expression'
      ERROR STOP 1
   ENDIF
   IF (rc /= 0) THEN
      WRITE(*,*) 'aed_ptm: ERROR allocating ', TRIM(what), &
                 ' : requested ', n, ' elements, stat = ', rc
      ERROR STOP 1
   ENDIF
END SUBROUTINE ptm_alloc_check
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION ptm_total_slots() RESULT(n)
!-------------------------------------------------------------------------------
! aed_n_particles * aed_n_groups, checked for overflow BEFORE it is formed.
! Both are default INTEGERs, so the product silently wraps past ~2.1e9 - and the
! wrapped value is then used as an array extent and as the upper bound of the
! loop that fills it.
!-------------------------------------------------------------------------------
!BEGIN
   IF (aed_n_groups > 0 .AND. aed_n_particles > HUGE(n)/aed_n_groups) THEN
      WRITE(*,*) 'aed_ptm: ERROR aed_n_particles * aed_n_groups overflows a default INTEGER: ', &
                 aed_n_particles, ' x ', aed_n_groups
      ERROR STOP 1
   ENDIF
   n = aed_n_particles*aed_n_groups
END FUNCTION ptm_total_slots
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_ptm_init_impl(ng,np,parts,n_ptm_vars_,n_cells)
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
   INTEGER :: n_slots
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
     IF (rc /= 0) ERROR STOP 'allocate_memory(): ERROR allocating (ptm_istat)'
   ALLOCATE(ptm_env(aed_n_groups,1:aed_n_particles,1:n_ptm_env+n_ptm_vars),stat=rc)
     IF (rc /= 0) ERROR STOP 'allocate_memory(): ERROR allocating (ptm_env)'
   ALLOCATE(ptm_state(aed_n_groups,1:aed_n_particles,1:n_ptm_vars),stat=rc)
     IF (rc /= 0) ERROR STOP 'allocate_memory(): ERROR allocating (ptm_state)'
   ALLOCATE(ptm_diag(aed_n_groups,1:aed_n_particles,1:n_ptm_vars),stat=rc) ! Not yet used
     IF (rc /= 0) ERROR STOP 'allocate_memory(): ERROR allocating (ptm_diag)'

   !#---------------------------------------------------------------------------
   !# Now set initial values
   ptm_istat(:,:,:) = -9999
   ptm_env(:,:,:)   = -9999.
   IF (n_ptm_vars > 0) THEN
      ptm_state(:,:,:) = -9999.
      ptm_diag(:,:,:)  = -9999.
   ENDIF

   IF (n_ptm_vars == 0) THEN
      ALLOCATE(all_particles(n_cells),stat=rc)
        CALL ptm_alloc_check(rc, n_cells, 'all_particles (no ptm vars)')
      RETURN
   ENDIF

   pv = 0;
   DO av=1,n_aed_vars_
      IF ( .NOT.  aed_get_var(av, tvar) ) ERROR STOP "     ERROR getting variable info"
      IF ( tvar%var_type == V_PARTICLE ) THEN  !# ptm variable
          pv = pv + 1

         !print *,'PTM',pv,n_ptm_env,tvar%initial
          ptm_state(:,:,pv) = tvar%initial !# Note this is all particles, regardless of status
                                           !#    (ptm_env(:,:,n_ptm_env+1:n_ptm_env+n_ptm_vars))
          ptm_env(:,:,n_ptm_env+pv) = tvar%initial
      ENDIF
   ENDDO

   n_slots = ptm_total_slots()
   ALLOCATE(ptm(n_slots),stat=rc)
     CALL ptm_alloc_check(rc, n_slots, 'ptm (aed_ptm_init)')

   DO grp=1,aed_n_groups
     DO prt=1,aed_n_particles
         ppid = (grp-1)*aed_n_particles + prt
         ! Point single particle object to the global particle data structure
         ptm(ppid)%ptm_group = grp
         ptm(ppid)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(ppid)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(ppid)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_env+n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(ppid)%ptm_diag  => ptm_diag(grp,prt,:)
     ENDDO
   ENDDO !end particle loop

   ppid = n_slots

   IF (n_ptm_vars > 0) CALL aed_initialize_particle(ppid,ptm)
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
   ALLOCATE(all_particles(n_cells),stat=rc)
     CALL ptm_alloc_check(rc, n_cells, 'all_particles (aed_ptm_init)')
END SUBROUTINE aed_ptm_init_impl

!###############################################################################
SUBROUTINE aed_ptm_init(ng,np,parts,n_ptm_vars_,n_cells)   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 10 (aed_api_guard.c); the body is aed_ptm_init_impl, unchanged
   CINTEGER, INTENT(in) :: ng,np
   INTEGER,  INTENT(in) :: n_ptm_vars_, n_cells
   TYPE(aed_part_group_t),DIMENSION(:),TARGET,INTENT(in) :: parts
   CALL aed_api_guard_enter(10)
   CALL aed_ptm_init_impl(ng,np,parts,n_ptm_vars_,n_cells)
   CALL aed_api_guard_leave(10)
END SUBROUTINE aed_ptm_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Particles_impl(n_cells)
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
   INTEGER :: cell, j
   INTEGER :: rc
   AED_REAL :: dt = 3600

   TYPE (aed_ptm_t) :: ptm
!
!-------------------------------------------------------------------------------

!BEGIN
   IF (aed_n_groups == 0 .OR. aed_n_particles == 0 .OR. n_ptm_vars == 0) RETURN
   zz = zero_

!------------
   ! Bin every active particle into all_particles(cell), the per-cell lists that
   ! aed_calculate_particles consumes. Two bugs fixed here (2026-08):
   !
   ! 1. STAT predicate. glm_ptm.c writes STAT=0 for an inactive/dead slot and STAT=1
   !    for active (see e.g. ptm_init_glm, ptm_removeparticles); STAT is never negative
   !    at runtime (see the -9999 fill in aed_ptm_init - it exists to mark an
   !    unassigned PTID, not STAT, and every STAT slot is overwritten with 0 before the
   !    first timestep). >= 0 / < 0 therefore tested the wrong sentinel: >= 0 admitted
   !    BOTH active and inactive particles, and the `< 0 CYCLE` below never fired.
   !    Corrected to > 0 / <= 0.
   !
   ! 2. Counting and populating must be TWO FULL PASSES OVER ALL GROUPS, not nested
   !    per-group. The previous structure allocated all_particles(cell)%prt sized at
   !    whichever group's count reached that cell first (guarded by
   !    .NOT. ALLOCATED(...%prt)), then never re-sized it as later groups added to
   !    %count - so with more than one group, particles from every group after the
   !    first were silently dropped by the `j <= %count` bound below. Harmless with a
   !    single group (the only configuration this has ever been run with until
   !    num_phytos > 1), which is why it went unnoticed.
   DO cell=1, size(all_particles)
      IF (ALLOCATED(all_particles(cell)%prt)) DEALLOCATE(all_particles(cell)%prt)
      all_particles(cell)%count = 0
      all_particles(cell)%n = 0   !# reset unconditionally - was left stale for any
                                   !# cell that emptied out since the previous step,
                                   !# since it was previously reset only inside the
                                   !# "not yet allocated" branch of the populate pass.
   ENDDO

   ! Pass 1: count every active particle, across ALL groups, before allocating anything.
   DO grp=1,aed_n_groups
      DO prt=1,aed_n_particles
         IF ( ptm_istat(grp,prt,STAT) > 0 ) THEN
            cell = ptm_istat(grp,prt,IDX3)
            IF ( cell >= 1 .AND. cell <= n_cells ) THEN
               all_particles(cell)%count = all_particles(cell)%count + 1
            ENDIF
         ENDIF
      ENDDO
   ENDDO

   DO cell=1, size(all_particles)
      IF (all_particles(cell)%count > 0) THEN
         ALLOCATE(all_particles(cell)%prt(all_particles(cell)%count),stat=rc)
           CALL ptm_alloc_check(rc, all_particles(cell)%count, 'all_particles(:)%prt')
      ENDIF
   ENDDO

   ! Pass 2: now that every cell's list is sized for its TOTAL across all groups,
   ! populate it - a second full pass over all groups, not resumed from pass 1.
   DO grp=1,aed_n_groups
      DO prt=1,aed_n_particles
         IF ( ptm_istat(grp,prt,STAT) <= 0 ) CYCLE  !# ignore inactive particles (must match pass 1)

         cell = ptm_istat(grp,prt,IDX3)
         IF ( cell >= 1 .AND. cell <= size(all_particles) ) THEN
            j = all_particles(cell)%n + 1              ! add new particle
            IF (j <= all_particles(cell)%count ) THEN
               all_particles(cell)%prt(j)%grp = grp
               all_particles(cell)%prt(j)%idx = prt
               all_particles(cell)%n = j
            ENDIF
         ENDIF
      ENDDO
   ENDDO
!-------

END SUBROUTINE Particles_impl

!###############################################################################
SUBROUTINE Particles(n_cells)   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 11 (aed_api_guard.c); the body is Particles_impl, unchanged
   INTEGER,  INTENT(in)    :: n_cells
   CALL aed_api_guard_enter(11)
   CALL Particles_impl(n_cells)
   CALL aed_api_guard_leave(11)
END SUBROUTINE Particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_ptm_set_cell_map_impl(cell_id)
!-------------------------------------------------------------------------------
! Provide the host map (layer,col) -> global cell index so multi-column hosts
! (e.g. ELCOM) can translate a column's local levels into the global cell space
! the particles are binned in (IDX3). Single-column hosts (GLM) need not call
! this; the global cell then defaults to the layer index.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,DIMENSION(:,:),INTENT(in) :: cell_id   !# (layer, col)
!LOCALS
   INTEGER :: rc
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (ALLOCATED(ptm_cell_id)) DEALLOCATE(ptm_cell_id)
   ALLOCATE(ptm_cell_id(size(cell_id,1),size(cell_id,2)),stat=rc)
   IF (rc /= 0) ERROR STOP 'aed_ptm_set_cell_map(): ERROR allocating (ptm_cell_id)'
   ptm_cell_id = cell_id
   have_cell_id = .TRUE.
END SUBROUTINE aed_ptm_set_cell_map_impl

!###############################################################################
SUBROUTINE aed_ptm_set_cell_map(cell_id)   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 14 (aed_api_guard.c); the body is aed_ptm_set_cell_map_impl, unchanged
   INTEGER,DIMENSION(:,:),INTENT(in) :: cell_id
   CALL aed_api_guard_enter(14)
   CALL aed_ptm_set_cell_map_impl(cell_id)
   CALL aed_api_guard_leave(14)
END SUBROUTINE aed_ptm_set_cell_map
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_particles_impl(icolm, nlev, idx_lo, idx_hi, col_no, dt_days)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles for a single column.
!
! Particles are binned (in Particles()) into all_particles() by their GLOBAL
! cell index (IDX3). idx_lo:idx_hi is the layer range of THIS column; local level
! `lev` maps to the host layer `idx_lo+lev-1`. For multi-column hosts the host
! supplies ptm_cell_id(layer,col) (via aed_ptm_set_cell_map) to convert that to
! the global cell. For single-column hosts (GLM) no map is set and the global
! cell defaults to the layer index - behaviour is unchanged. Particle splitting
! is handled separately, once per step, by aed_split_particles().
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_column_t),INTENT(inout) :: icolm(:)
   INTEGER,INTENT(in)               :: nlev          !# layers in this column (= idx_hi-idx_lo+1)
   INTEGER,INTENT(in)               :: idx_lo,idx_hi !# host layer index range of this column
   INTEGER,INTENT(in)               :: col_no        !# host column number (for the cell map)
   !# Biological timestep for THIS call, in days. The host must pass the
   !# interval it is actually advancing over - the full host step if this is
   !# called once per step, or dt_eff/secs_per_day if called once per
   !# split_factor sub-step. Omitting it leaves the previous hourly default.
   AED_REAL,INTENT(in),OPTIONAL     :: dt_days
!
!LOCAL VARIABLES:
   INTEGER :: lev, lyr, gc, grp, prt, pt
   INTEGER :: ppid
   INTEGER :: rc
   TYPE (aed_ptm_t), DIMENSION(:), ALLOCATABLE :: ptm
   TYPE(partgroup_cell), POINTER :: layer_particles
!
!-------------------------------------------------------------------------------

!BEGIN
   !# Publish this call's biological timestep for the particle BGC models to
   !# read (aed_core%aed_ptm_dt_days). Set before the models run, not after; the
   !# flag lets consumers tell a published value from the module default.
   IF ( PRESENT(dt_days) ) THEN
      aed_ptm_dt_days = dt_days
      aed_ptm_dt_days_set = .TRUE.
   ENDIF

   IF (aed_n_groups == 0 .OR. aed_n_particles == 0 .OR. n_ptm_vars == 0) RETURN

   DO lev=1,nlev
      lyr = idx_lo + lev - 1           !# host layer index for this local level
      IF ( have_cell_id ) THEN
         IF ( lyr < 1 .OR. lyr > size(ptm_cell_id,1) .OR. &
              col_no < 1 .OR. col_no > size(ptm_cell_id,2) ) CYCLE
         gc = ptm_cell_id(lyr, col_no) !# global cell index
      ELSE
         gc = lyr                      !# single-column fallback (GLM)
      ENDIF
      IF ( gc < 1 .OR. gc > size(all_particles) ) CYCLE
      layer_particles => all_particles(gc)

      ! Iterate %n (the number of slots pass 2 of Particles() actually FILLED), not
      ! %count (the number pass 1 allocated %prt for). The two now always agree since
      ! Particles() was fixed to a genuine two-pass count-then-populate structure, but
      ! %n is the correct one to depend on: %prt is a freshly-ALLOCATEd INTEGER array,
      ! and any future divergence between the passes would otherwise feed an
      ! uninitialised subscript straight into the pointer associations below.
      IF (layer_particles%n == 0) CYCLE

      IF (layer_particles%n == 0) CYCLE
      ALLOCATE(ptm(layer_particles%n),stat=rc)
        CALL ptm_alloc_check(rc, layer_particles%n, 'ptm (aed_calculate_particles)')

      DO pt=1,layer_particles%n

         ! Retrieve particle properties, from the particle-cell object
         grp = layer_particles%prt(pt)%grp ; prt = layer_particles%prt(pt)%idx

         ! Point single particle object to the global particle data structure
         ptm(pt)%ptm_group = grp
         ptm(pt)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(pt)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(pt)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_env+n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(pt)%ptm_diag  => ptm_diag(grp,prt,:)

      ENDDO !end particle loop
      ppid = layer_particles%n

      ! Pass the particles in this cell to AED modules
      CALL aed_particle_bgc(icolm,lev,ppid,p=ptm) ! Note: ppid getting incremeted in here

      DEALLOCATE(ptm)
   ENDDO !end cell loop

END SUBROUTINE aed_calculate_particles_impl

!###############################################################################
SUBROUTINE aed_calculate_particles(icolm, nlev, idx_lo, idx_hi, col_no, dt_days)   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 12 (aed_api_guard.c); the body is aed_calculate_particles_impl, unchanged
   TYPE(aed_column_t),INTENT(inout) :: icolm(:)
   INTEGER,INTENT(in)               :: nlev, idx_lo, idx_hi, col_no
   AED_REAL,INTENT(in),OPTIONAL     :: dt_days
   CALL aed_api_guard_enter(12)
   CALL aed_calculate_particles_impl(icolm, nlev, idx_lo, idx_hi, col_no, dt_days)
   CALL aed_api_guard_leave(12)
END SUBROUTINE aed_calculate_particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_split_particles_impl()
!-------------------------------------------------------------------------------
!
! Apply particle splitting (ABM) across ALL particle groups, once per step.
!
! Previously this was done inside aed_calculate_particles, which is called once
! per column - for multi-column hosts that would split particles N_cols times.
! This routine is host-column independent and must be called exactly once after
! all columns have been processed.
!-------------------------------------------------------------------------------
!LOCAL VARIABLES:
   INTEGER :: grp, prt, ppid
   INTEGER :: rc, n_slots
   TYPE (aed_ptm_t), DIMENSION(:), ALLOCATABLE :: ptm
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (aed_n_groups == 0 .OR. aed_n_particles == 0 .OR. n_ptm_vars == 0) RETURN

   n_slots = ptm_total_slots()
   ALLOCATE(ptm(n_slots),stat=rc)
     CALL ptm_alloc_check(rc, n_slots, 'ptm (aed_split_particles)')

   DO grp=1,aed_n_groups
     DO prt=1,aed_n_particles
         ppid = (grp-1)*aed_n_particles + prt
         ! Point single particle object to the global particle data structure
         ptm(ppid)%ptm_group = grp
         ptm(ppid)%ptm_istat => ptm_istat(grp,prt,:)
         ptm(ppid)%ptm_env   => ptm_env(grp,prt,1:n_ptm_env)
         ptm(ppid)%ptm_state => ptm_env(grp,prt,n_ptm_env+1:n_ptm_env+n_ptm_vars)    !ptm_state(grp,prt,:)
         ptm(ppid)%ptm_diag  => ptm_diag(grp,prt,:)
     ENDDO
   ENDDO !end particle loop

   ppid = n_slots

   CALL aed_split_particle(ppid,ptm)
   DEALLOCATE(ptm)

END SUBROUTINE aed_split_particles_impl

!###############################################################################
SUBROUTINE aed_split_particles()   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 13 (aed_api_guard.c); the body is aed_split_particles_impl, unchanged
   CALL aed_api_guard_enter(13)
   CALL aed_split_particles_impl()
   CALL aed_api_guard_leave(13)
END SUBROUTINE aed_split_particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Particles_zz(column, count, parts)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles
!
! DEAD CODE (2026-08). Not in this module's PUBLIC list and called from nowhere in
! the tree - its only real work call (aed_particle_bgc, below) is commented out, so
! it currently does nothing but copy values out of particle_groups and back again.
! Retained for reference only.
!
! It carried the same two bugs fixed in Particles() above, and they are corrected
! here too so this routine cannot reintroduce them if it is ever revived:
!   1. the status predicate was `>= 0`, admitting inactive particles (status 0);
!      every live consumer tests `> 0` (see the STAT predicate note in Particles()).
!   2. the per-cell loop ran to parts(lev)%count (the ALLOCATED size) rather than
!      parts(lev)%n (the FILLED count) - the same %count-vs-%n hazard fixed in
!      aed_calculate_particles above.
! If this is revived, re-verify both against the current Particles() /
! aed_calculate_particles first.
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

      DO pt=1,parts(lev)%n            !# was %count - see the DEAD CODE note above

         grp = parts(lev)%prt(pt)%grp ; prt = parts(lev)%prt(pt)%idx
         stat = particle_groups(grp)%idx_stat   ! should be 1
         idxi3 =  particle_groups(grp)%idx_3   ! should be 3

         IF ( particle_groups(grp)%status(stat, prt) > 0 ) THEN   !# was >= 0 - see above
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
SUBROUTINE set_ptm_aed_var_num_impl(n_aed_vars)
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
END SUBROUTINE set_ptm_aed_var_num_impl

!###############################################################################
SUBROUTINE set_ptm_aed_var_num(n_aed_vars)   !# S5 wrapper
!# S5 2026-09-15: guarded entry point 15 (aed_api_guard.c); the body is set_ptm_aed_var_num_impl, unchanged
   INTEGER, INTENT(in) :: n_aed_vars
   CALL aed_api_guard_enter(15)
   CALL set_ptm_aed_var_num_impl(n_aed_vars)
   CALL aed_api_guard_leave(15)
END SUBROUTINE set_ptm_aed_var_num
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_ptm
