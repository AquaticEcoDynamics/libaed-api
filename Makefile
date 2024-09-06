###############################################################################
#                                                                             #
# Makefile to build libaed-api                                                #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2024 -  The University of Western Australia                      #
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
###############################################################################

VERSION=$(shell grep AED_API_VERSION include/aed_api.h | head -1 | cut -f2 -d\")
SOVERS=$(shell echo $(VERSION) | cut -f1 -d\.)
VERS=$(shell echo $(VERSION) | cut -f2- -d\.)

LIBAEDAPI=aed-api
OUTLIB=lib$(LIBAEDAPI)

INCLUDES=-I../libaed-water/${incdir}  -I../libaed-water/${moddir}

ifeq ($(AEDWATDIR),)
  AEDWATDIR=../libaed-water
endif
INCLUDES+=-I${AEDWATDIR}/include -I${incdir} -I${AEDWATDIR}/mod

ifeq ($(OSTYPE),Darwin)
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else
  SHARED=-shared -Wl,-soname,$(OUTLIB).so.$(SOVERS)
  so_ext=so
endif


include ../libaed-water/make_defs.inc

OBJS=${objdir}/aed_bivalve.o \
     ${objdir}/aed_habitat_benthic.o \
     ${objdir}/aed_habitat_ruppia.o \
     ${objdir}/aed_habitat_chara.o \
     ${objdir}/aed_habitat_galaxiid.o \
     ${objdir}/aed_habitat_seagrass.o \
     ${objdir}/aed_macroalgae.o \
     ${objdir}/aed_macroalgae2.o \
     ${objdir}/aed_macrophyte.o \
     ${objdir}/aed_benthic.o


SOFLAGS = ${libdir}/${OUTLIB}.a

EXTFLAG=
ifneq ($(AEDBENDIR),)
  LIBAEDBBEN=aed-benthic
  SOFLAGS+=${AEDBENDIR}/lib/lib${LIBAEDBEN}.a
endif
ifneq ($(AEDRIPDIR),)
  LIBAEDBRIP=aed-riparian
  SOFLAGS+=${AEDRIPDIR}/lib/lib${LIBAEDRIP}.a
else
  EXTFLAG+=-DNO_RIPARIAN
endif
ifneq ($(AEDDMODIR),)
  LIBAEDBDMO=aed-demo
  SOFLAGS+=${AEDDMODIR}/lib/lib${LIBAEDDMO}.a
endif
ifneq ($(AEDLGTDIR),)
  LIBAEDBLGT=aed-lighting
  SOFLAGS+=${AEDLGTDIR}/lib/lib${LIBAEDLGT}.a
else
  EXTFLAG+=-DNO_LGT
endif
ifneq ($(AEDDEVDIR),)
  LIBAEDBDEV=aed-dev
  SOFLAGS+=${AEDDEVDIR}/lib/lib${LIBAEDDEV}.a
else
  EXTFLAG+=-DNO_DEV
endif

OBJS=${objdir}/aed_zones.o \
     ${objdir}/aed_api.o

ifeq ($(EXTERNAL_LIBS),shared)
  FFLAGS+=-fPIC
  TARGET = ${libdir}/$(OUTLIB).${so_ext}
else
  FFLAGS+=-fPIE
  TARGET = ${libdir}/$(OUTLIB).a
endif

include ../libaed-water/make_rules.inc

${libdir}/${OUTLIB}.${so_ext}: ${libdir}/${OUTLIB}.a ${OBJS}
	$(F90) ${SHARED} -o $@.${SOVERS}.${VERS} ${OBJS} ${LDFLAGS} ${SOFLAGS}
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@.${SOVERS}
