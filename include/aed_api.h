!###############################################################################
!#                                                                             #
!# aed_api.h                                                                   #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2024 - The University of Western Australia                       #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#                                                                             #
!# NB - This is currently a "work-in-progress" do not rely on it's function    #
!#      declarations/constants/whatever because it will almost certainly be    #
!#      changed.                                                               #
!#                                                                             #
!###############################################################################
#ifndef _AED_API_H_
#define _AED_API_H_

#include "aed.h"
#include "aed_api_env.h"

#define AED_API_VERSION  0.9.1"

#ifndef AED_REAL
#  define AED_REAL REAL(kind=C_DOUBLE)
#endif
! #define CAED_REAL REAL(kind=C_DOUBLE)
! #define NF90_REALTYPE NF90_DOUBLE
! #define NC_FILLER NC_FILL_DOUBLE
! #define IFIX IDINT
! #define AMOD DMOD
! #define ALOG10 DLOG10
! #define EXP DEXP
! #define AINT DINT
! #define FLOAT
#define DOUBLETYPE double precision
#define CINTEGER INTEGER(kind=C_INT)
#define CSIZET   INTEGER(kind=C_SIZE_T)
#define CLOGICAL LOGICAL(kind=C_BOOL)
#define CCHARACTER CHARACTER(C_CHAR)


#ifndef isnan
#  define isnan(x) ieee_is_nan(x)
#  define HAVE_IEEE_ARITH
#endif

#if 0

INTERFACE

  !#############################################################################
  SUBROUTINE aed_config_model(conf)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     TYPE(api_config_t), INTENT(in) :: conf
  END SUBROUTINE aed_config_model
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE aed_set_model_env(env)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     TYPE(api_env_t),INTENT(in) :: env
  END SUBROUTINE aed_set_model_env
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE aed_set_model_data(dat)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     TYPE(api_data_t), INTENT(in) :: dat
  END SUBROUTINE aed_set_model_data
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE aed_set_mobility(mobl) 
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     PROCEDURE(aed_mobility_t),POINTER :: mobl
  END SUBROUTINE aed_set_mobility
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !*****************************************************************************
  !* initialise the aed libraries: provide an aed.nml file name - returns the
  !* numbers of variables of different types
  !*****************************************************************************

  !#############################################################################
  SUBROUTINE aed_init_model(fname, NumWQ_Vars, NumWQ_Ben, NumWQ_Diag, NumWQ_DiagL)
  !-----------------------------------------------------------------------------
  ! Initialize the AED driver by reading settings from "fname".
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     CHARACTER(*),INTENT(in) :: fname
     INTEGER,INTENT(out)     :: NumWQ_Vars, NumWQ_Ben
     INTEGER,INTENT(out)     :: NumWQ_Diag, NumWQ_DiagL
  !
  END SUBROUTINE aed_init_model
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !*****************************************************************************
  !* run the model for one timestep over "wlev" levels
  !* the "doSurface" parameter was added to alow for ice coverage where there
  !* would be no interaction with the atmosphere
  !*****************************************************************************

  !#############################################################################
  SUBROUTINE aed_run_model(wlev, doSurface)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     INTEGER,INTENT(in) :: wlev
     LOGICAL,INTENT(in) :: doSurface
  !
  END SUBROUTINE aed_run_model
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !*****************************************************************************
  !* For zones.
  !*****************************************************************************

  !#############################################################################
  SUBROUTINE aed_init_zones(n_zones, n_levs, z_cc, z_cc_hz, z_diag, z_diag_hz)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     INTEGER,INTENT(in) :: n_zones, n_levs
     AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_cc      !(n_zones, n_levs, n_vars)
     AED_REAL,DIMENSION(:,:),  POINTER,INTENT(in) :: z_cc_hz   !(n_zones+1, n_vars)
     AED_REAL,DIMENSION(:,:,:),POINTER,INTENT(in) :: z_diag    !(n_zones, n_levs, n_vars)
     AED_REAL,DIMENSION(:,:)  ,POINTER,INTENT(in) :: z_diag_hz !(n_zones+1, n_vars)
  END SUBROUTINE aed_init_zones
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE api_set_zone_funcs(copy_to, copy_from, calc_areas)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     PROCEDURE(copy_to_zone_t),POINTER    :: copy_to
     PROCEDURE(copy_from_zone_t),POINTER  :: copy_from
     PROCEDURE(calc_zone_areas_t),POINTER :: calc_areas
  END SUBROUTINE api_set_zone_funcs
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !*****************************************************************************
  !* Clean up any data allocated.
  !*****************************************************************************

  !#############################################################################
  SUBROUTINE aed_clean_model()
  !-----------------------------------------------------------------------------
  END SUBROUTINE aed_clean_model
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END INTERFACE

#endif

#endif
