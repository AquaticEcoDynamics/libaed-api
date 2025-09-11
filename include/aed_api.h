/*#############################################################################*
 #                                                                             #
 # aed_api.h                                                                   #
 #                                                                             #
 #  Developed by :                                                             #
 #      AquaticEcoDynamics (AED) Group                                         #
 #      School of Agriculture and Environment                                  #
 #      The University of Western Australia                                    #
 #                                                                             #
 #      http://aquatic.science.uwa.edu.au/                                     #
 #                                                                             #
 #  Copyright 2024-2025 - The University of Western Australia                  #
 #                                                                             #
 #  This file is part of libaed (Library for AquaticEco Dynamics)              #
 #                                                                             #
 #   AED is free software: you can redistribute it and/or modify               #
 #   it under the terms of the GNU General Public License as published by      #
 #   the Free Software Foundation, either version 3 of the License, or         #
 #   (at your option) any later version.                                       #
 #                                                                             #
 #   AED is distributed in the hope that it will be useful,                    #
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
 #   GNU General Public License for more details.                              #
 #                                                                             #
 #   You should have received a copy of the GNU General Public License         #
 #   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
 #                                                                             #
 #                                                                             #
 # NB - This is currently a "work-in-progress" do not rely on it's function    #
 #      declarations/constants/whatever because it will almost certainly be    #
 #      changed.                                                               #
 #                                                                             #
 *#############################################################################*/
#ifndef _AED_API_H_
#define _AED_API_H_

#include "aed.h"
#include "aed_api_env.h"

#define AED_API_VERSION  "0.9.5"

#ifndef __STDC__
#ifndef AED_REAL
#  define AED_REAL REAL(kind=C_DOUBLE)
#endif

#ifndef FLOAT
#define FLOAT(x) (x)
#endif

#define CINTEGER INTEGER(kind=C_INT_T)
#define CSIZET   INTEGER(kind=C_SIZE_T)
#define CCHARACTER CHARACTER(C_CHAR)
#endif

#endif
