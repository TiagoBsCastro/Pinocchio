#  *****************************************************************
#  *                        PINOCCHIO  V4.1                        *
#  *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
#  *****************************************************************
#  
#  This code was written by
#  Pierluigi Monaco
#  Copyright (C) 2016
#  
#  web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


# Directory settings #

RUNDIR = ../example/

# Executable name #

EXEC ?= pinocchio.x

# Compilation options #

DEBUG ?= NO

########## Preprocessor Definitions ########## 

# Displacement LPT order # 

OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# Dynamics of triaxial collapse #

OPTIONS += -DELL_CLASSIC
OPTIONS += -DNORADIATION

####################### CPU flag #####################
ifeq ($(DEBUG), YES)
OPTIONS += -DCPU_DEBUG
endif
########## System-specific Configuration ############
SYSTYPE ?= "my_machine"
###################### MY_MACHINE ####################

ifeq ($(SYSTYPE),"my_machine")
CC          =  mpicc 
CDEBUG      =  -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  =  -g -O3 -Wno-unused-result -march=native $(OMP_FLAG)
FFTW_LIBR   = -L/usr/local/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/usr/local/include 
PFFT_LIBR   = -L/usr/local/lib -lpfft
PFFT_INCL   = 
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L/usr/local/lib -lgsl -lgslcblas
GSL_INCL    = -I/usr/local/include
endif

########## System-specific Configuration end ########## 

# Include paths for libraries #

INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)

# Library flags #

LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR)

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Source files #

OBJECTS = \
	  variables.o \
	  cosmo.o \
	  initialization.o \
	  fmax.o\
	  fmax-pfft.o \
	  GenIC.o \
	  ReadParamfile.o \
	  allocations.o \
	  LPT.o \
	  collapse_times.o 
#write_snapshot.o fragment.o distribute.o build_groups.o write_halos.o collapse_times_GPU.o cubic_spline_interpolation.o

# Main targets and rules

pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)

# Generic rule for compiling source files #

%.o: %.c Makefile pinocchio.h #def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

# Clean target #

clean:
	\rm -f *.o *~ $(EXEC)


