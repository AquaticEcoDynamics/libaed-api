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
#  Copyright 2024-2026 : The University of Western Australia                  #
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

LIBAEDAPI=aed-api
APILIB=lib$(LIBAEDAPI)

VERSION=$(shell grep AED_API_VERSION include/aed_api.h | head -1 | cut -f2 -d\")

INCLUDES=-I../libaed-water/${incdir}  -I../libaed-water/${moddir}

ifeq ($(AEDWATDIR),)
  AEDWATDIR=../libaed-water
endif
INCLUDES+=-I${AEDWATDIR}/include -I${incdir} -I${AEDWATDIR}/mod

include ../libaed-water/make_defs.inc

EXTFFLAGS=
#EXTFFLAGS=-DNO_RIPARIAN -DNO_LGT -DNO_DEV
ifneq ($(AEDBENDIR),)
  LIBAEDBBEN=aed-benthic
  SOFLAGS+=${AEDBENDIR}/lib/lib${LIBAEDBEN}.a
endif
ifneq ($(AEDDMODIR),)
  LIBAEDBDMO=aed-demo
  SOFLAGS+=${AEDDMODIR}/lib/lib${LIBAEDDMO}.a
endif
ifneq ($(AEDRIPDIR),)
  LIBAEDBRIP=aed-riparian
  SOFLAGS+=${AEDRIPDIR}/lib/lib${LIBAEDRIP}.a
else
  EXTFFLAGS+=-DNO_RIPARIAN
endif
ifneq ($(AEDLGTDIR),)
  LIBAEDBLGT=aed-lighting
  SOFLAGS+=${AEDLGTDIR}/lib/lib${LIBAEDLGT}.a
else
  EXTFFLAGS+=-DNO_LGT
endif
ifneq ($(AEDDEVDIR),)
  LIBAEDBDEV=aed-dev
  SOFLAGS+=${AEDDEVDIR}/lib/lib${LIBAEDDEV}.a
else
  EXTFFLAGS+=-DNO_DEV
endif

OBJS=${objdir}/aed_zones.o \
     ${objdir}/aed_ptm.o   \
     ${objdir}/aed_api.o

#WOBJS=$(shell ar t ../libaed-water/lib/libaed-water.a | sed -e 's!^!../libaed-water/obj/!')
#BOBJS=$(shell ar t ../libaed-benthic/lib/libaed-benthic.a | sed -e 's!^!../libaed-benthic/obj/!')
#DOBJS=$(shell ar t ../libaed-demo/lib/libaed-demo.a | sed -e 's!^!../libaed-demo/obj/!')

#OBJS+=${WOBJS} ${BOBJS} ${DOBJS}

ifeq ($(EXTERNAL_LIBS),shared)
# FFLAGS+=-fPIC
  TARGET = ${libdir}/$(APILIB).${so_ext}
else
#  FFLAGS+=-fPIE
   TARGET = ${libdir}/$(APILIB).a
endif

#LIBS=${SOFLAGS}
#LIBS=../libaed-water/lib/libaed-water.a ../libaed-benthic/lib/libaed-benthic.a ../libaed-demo/lib/libaed-demo.a

include ../libaed-water/make_rules.inc

${objdir}/aed_external.o: ../libaed-water/src/aed_external.F90
	$(F90) $(FFLAGS) $(EXTFFLAGS) $(OMPFLAG) -c $< -o $@

${libdir}/${APILIB}.${so_ext}: ${libdir}/${APILIB}.a
	$(F90) ${SHARED} -o $@.${SOVERS}.${VERS} ${OBJS} ${LDFLAGS} ${SOFLAGS}
	ln -sf ${APILIB}.${so_ext}.${SOVERS}.${VERS} $@
	ln -sf ${APILIB}.${so_ext}.${SOVERS}.${VERS} $@.${SOVERS}

${libdir}/$(APILIB).a: ${objdir} ${moddir} ${libdir} ${OBJS} ${LIBS}
	ar rv $@ ${OBJS} ${LIBS}
	ranlib $@

