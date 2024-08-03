###############################################################################
#                                                                             #
# Makefile to build the aed water quality library API                         #
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

objdir=obj
srcdir=src
libdir=lib
moddir=mod
incdir=include

OSTYPE=$(shell uname -s)

LIBAEDAPI=aed-api
OUTLIB=lib$(LIBAEDAPI)

ifeq ($(F90),)
  F90=gfortran
endif
ifeq ($(MDEBUG),true)
  DEBUG=true
endif

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

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OMPFLAG=-qopenmp
  OPT_FFLAGS=-O3 ${OMPFLAG}
  FFLAGS=-g -fpp -warn all -module ${moddir} -static-intel -fPIE -mp1 -warn nounused $(DEFINES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check all -check noarg_temp_created
  endif
  FFLAGS+=-real-size 64
else ifeq ($(F90),flang)
  ifeq ($(OSTYPE),FreeBSD)
    INCLUDES+=-I../ancillary/freebsd/mod
  endif
  DEBUG_FFLAGS=-g
  OMPFLAG=-fopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-module ${moddir} $(DEFINES) $(INCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
else
  FFLAGS+=-Wall -J ${moddir} -ffree-line-length-none
  FFLAGS+=-std=f2008 -fall-intrinsics -fdefault-real-8 -fdefault-double-8
# FFLAGS+=-Wno-unused-value -Wno-unused-dummy-argument -Wno-unused-function
# FFLAGS+=-Wno-unused-variable -Wno-c-binding-type
endif


LIBWATAED=aed-water
SOFLAGS = ${libdir}/lib${LIBAEDAPI}.a

EXTFLAG=
ifneq ($(AEDBENDIR),)
  LIBBENAED=aed-benthic
  SOFLAGS+=${AEDBENDIR}/lib/lib${LIBBENAED}.a
endif
ifneq ($(AEDRIPDIR),)
  LIBRIPAED=aed-riparian
  SOFLAGS+=${AEDRIPDIR}/lib/lib${LIBRIPAED}.a
else
  EXTFLAG+=-DNO_RIPARIAN
endif
ifneq ($(AEDDMODIR),)
  LIBDMOAED=aed-demo
  SOFLAGS+=${AEDDMODIR}/lib/lib${LIBDMOAED}.a
endif
ifneq ($(AEDLGTDIR),)
  LIBLGTAED=aed-lighting
  SOFLAGS+=${AEDLGTDIR}/lib/lib${LIBLGTAED}.a
else
  EXTFLAG+=-DNO_LGT
endif
ifneq ($(AEDDEVDIR),)
  LIBDEVAED=aed-dev
  SOFLAGS+=${AEDDEVDIR}/lib/lib${LIBDEVAED}.a
else
  EXTFLAG+=-DNO_DEV
endif

FFLAGS+=$(OPT_FFLAGS)

ifeq ($(DEBUG),true)
  FFLAGS+=$(DEBUG_FFLAGS)
endif

ifeq ($(PRECISION),1)
  TFFLAGS += -D_PRECISION=1
else ifeq ($(PRECISION),2)
  TFFLAGS += -D_PRECISION=2
else
  TFFLAGS += -D_PRECISION=1
endif

TFFLAGS += -g -DAED -DEXTERNAL_WQ=2
INCLUDES += -I${moddir}

FFLAGS+=${INCLUDES}

OBJS=${objdir}/aed_zones.o \
     ${objdir}/aed_api.o

ifeq ($(EXTERNAL_LIBS),shared)
  FFLAGS+=-fPIC
  TARGET = ${libdir}/$(OUTLIB).${so_ext}
else
  FFLAGS+=-fPIE
  TARGET = ${libdir}/$(OUTLIB).a
endif

all: ${TARGET}

${libdir}/lib${LIBAEDAPI}.a: ${objdir} ${moddir} ${libdir} ${OBJS}
	ar -rv $@ ${OBJS}
	ranlib $@

${libdir}/${OUTLIB}.${so_ext}: ${libdir}/lib${LIBAEDAPI}.a ${OBJS}
	$(F90) ${SHARED} -o $@.${SOVERS}.${VERS} ${OBJS} ${LDFLAGS} ${SOFLAGS}
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@.${SOVERS}

${objdir}/%.o: ${srcdir}/%.F90 ${incdir}/aed_api.h ${AEDWATDIR}/include/aed.h
	$(F90) -fPIC $(FFLAGS) $(EXTRA_FFLAGS) -D_FORTRAN_SOURCE_ -c $< -o $@

${objdir}/aed_api.o: ${srcdir}/aed_api.F90 ${objdir}/aed_zones.o
	$(F90) $(FFLAGS) $(EXTFLAG) ${INCLUDES} -g -c $< -o $@

${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

${libdir}:
	@mkdir ${libdir}

clean:
	/bin/rm -f *.i90
	/bin/rm -f ${objdir}/*.o
	/bin/rm -f ${moddir}/*.mod
	/bin/rm -f ${libdir}/*.a
	/bin/rm -f ${libdir}/*.${so_ext}*

distclean: clean
	/bin/rm -rf ${libdir} ${moddir} ${objdir}
