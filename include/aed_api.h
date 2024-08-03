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
!#  Copyright 2024 -  The University of Western Australia                      #
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

#include <aed.h>
#include <aed_api_env.h>

#define AED_API_VERSION  0.9.0"

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

  !*****************************************************************************
  !* The 3 environment configuration routines below declare which variables are
  !* to be made available. The "which" parameter is a constant defined in the
  !* header file  aed_api_env.h
  !* The data parameter is an optional pointer to the data, if it is not suplied
  !* the api will create one.
  !*****************************************************************************

  !#############################################################################
  SUBROUTINE aed_environ_3d_var(which, data)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     INTEGER,INTENT(in) :: which
   ! AED_REAL,DIMENSION(:,:),POINTER,INTENT(in),OPTIONAL :: data
     AED_REAL,DIMENSION(:),POINTER,INTENT(in),OPTIONAL :: data
  !
  END SUBROUTINE aed_environ_3d_var
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE aed_environ_2d_var(which, data)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     INTEGER,INTENT(in) :: which
     AED_REAL,DIMENSION(:),POINTER,INTENT(in),OPTIONAL   :: data
  !
  END SUBROUTINE aed_environ_2d_var
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  SUBROUTINE aed_environ_0d_var(which, data)
  !-----------------------------------------------------------------------------
  !ARGUMENTS
     INTEGER,INTENT(in) :: which
     AED_REAL,POINTER,INTENT(in),OPTIONAL :: data
  !
  END SUBROUTINE aed_environ_0d_var
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
  !* run the model for ne timestep over "wlev" levels
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
