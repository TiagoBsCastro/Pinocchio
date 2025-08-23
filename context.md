# Repository context for LLM

## How this was built
- Files are text only and truncated to a head limit.
- Common binary and build artifacts are excluded.
- Likely secret files are excluded.

## Directory tree (truncated)
```
Pinocchio/
- example/
  - CAMBFiles/
    - hubble.dat
    - pk_cb_000.dat
    - pk_cb_001.dat
    - pk_cb_002.dat
    - pk_cb_003.dat
    - pk_cb_004.dat
    - pk_cb_005.dat
    - pk_cb_006.dat
    - pk_cb_007.dat
    - pk_cb_008.dat
    - pk_cb_009.dat
    - pk_cb_010.dat
    - pk_cb_011.dat
    - pk_cb_012.dat
    - pk_cb_013.dat
    - pk_cb_014.dat
    - pk_cb_015.dat
    - pk_cb_016.dat
    - pk_cb_017.dat
    - pk_cb_018.dat
    - pk_cb_019.dat
    - pk_cb_020.dat
    - pk_cb_021.dat
    - pk_cb_022.dat
    - pk_cb_023.dat
    - pk_cb_024.dat
    - pk_cb_025.dat
    - pk_cb_026.dat
    - pk_cb_027.dat
    - pk_cb_028.dat
    - pk_cb_029.dat
    - pk_cb_030.dat
    - pk_cb_031.dat
    - pk_cb_032.dat
    - pk_cb_033.dat
    - pk_cb_034.dat
    - pk_cb_035.dat
    - pk_cb_036.dat
    - pk_cb_037.dat
    - pk_cb_038.dat
    - pk_cb_039.dat
    - pk_cb_040.dat
    - pk_cb_041.dat
    - pk_cb_042.dat
    - pk_cb_043.dat
    - pk_cb_044.dat
    - pk_cb_045.dat
    - pk_cb_046.dat
    - pk_cb_047.dat
    - pk_cb_048.dat
    - pk_cb_049.dat
    - pk_cb_050.dat
    - pk_cb_051.dat
    - pk_cb_052.dat
    - pk_cb_053.dat
    - pk_cb_054.dat
    - pk_cb_055.dat
    - pk_cb_056.dat
    - pk_cb_057.dat
    - pk_cb_058.dat
    - pk_cb_059.dat
    - pk_cb_060.dat
    - pk_cb_061.dat
    - pk_cb_062.dat
    - pk_cb_063.dat
    - pk_cb_064.dat
    - pk_cb_065.dat
    - pk_cb_066.dat
    - pk_cb_067.dat
    - pk_cb_068.dat
    - pk_cb_069.dat
    - pk_cb_070.dat
    - pk_cb_071.dat
    - pk_cb_072.dat
    - pk_cb_073.dat
    - pk_cb_074.dat
    - pk_cb_075.dat
    - pk_cb_076.dat
    - pk_cb_077.dat
    - pk_cb_078.dat
    - pk_cb_079.dat
    - pk_cb_080.dat
    - pk_cb_081.dat
    - pk_cb_082.dat
    - pk_cb_083.dat
    - pk_cb_084.dat
    - pk_cb_085.dat
    - pk_cb_086.dat
    - pk_cb_087.dat
    - pk_cb_088.dat
    - pk_cb_089.dat
    - pk_cb_090.dat
    - pk_cb_091.dat
    - pk_cb_092.dat
    - pk_cb_093.dat
    - pk_cb_094.dat
    - pk_cb_095.dat
    - pk_cb_096.dat
    - pk_cb_097.dat
    - pk_cb_098.dat
    - pk_cb_099.dat
    - redshifts.dat
  - lss.png
  - mf.png
  - outputs
  - parameter_file
  - pinocchio.0.0000.example.catalog.out
  - pinocchio.0.0000.example.mf.out
  - pinocchio.0.5000.example.catalog.out
  - pinocchio.0.5000.example.mf.out
  - pinocchio.1.0000.example.catalog.out
  - pinocchio.1.0000.example.mf.out
  - pinocchio.2.0000.example.catalog.out
  - pinocchio.2.0000.example.mf.out
  - pinocchio.example-4.cosmology.out
  - pinocchio.example-4.FmaxPDF.out
  - pinocchio.example-4.geometry.out
  - pinocchio.example-4.histories.out
  - pinocchio.example-4.massmaps.positions.lattice.summary.out
  - pinocchio.example-4.massmaps.positions.summary.out
  - pinocchio.example-4.massmaps.positions.task00000.out
  - pinocchio.example-4.massmaps.positions.task00001.out
  - pinocchio.example-4.massmaps.positions.task00002.out
  - pinocchio.example-4.massmaps.positions.task00003.out
  - pinocchio.example-4.nz.out
  - pinocchio.example-4.plc.out
  - pinocchio.example-4.scaledep.out
  - pinocchio.example.cosmology.out
  - pinocchio.example.FmaxPDF.out
  - pinocchio.example.geometry.out
  - pinocchio.example.histories.out
  - pinocchio.example.massmaps.positions.lattice.summary.out
  - pinocchio.example.massmaps.positions.summary.out
  - pinocchio.example.massmaps.positions.task00000.out
  - pinocchio.example.nz.out
  - pinocchio.example.plc.out
  - pinocchio.example.scaledep.out
  - pinocchio.x
  - plc.png
  - PlotExample.py
- HMF_Validation/
  - log_RUN.txt
  - MHF_Validation_with_Watson_fit.png
  - outputs
  - parameter_file
  - pinocchio.0.0000.test.catalog.out
  - pinocchio.0.0000.test.mf.out
  - pinocchio.0.5000.test.catalog.out
  - pinocchio.0.5000.test.mf.out
  - pinocchio.1.0000.test.catalog.out
  - pinocchio.1.0000.test.mf.out
  - pinocchio.2.0000.test.catalog.out
  - pinocchio.2.0000.test.mf.out
  - pinocchio.test.cosmology.out
  - pinocchio.test.FmaxPDF.out
  - pinocchio.test.histories.out
  - pinocchio.x
  - VALIDATION_log.txt
- scripts/
  - benchmark_mass_maps.sh
  - HMF_validation.py
  - outputs
  - parameter_file
  - Pinocchio2fits.py
  - PkCamb.py
  - PlcGeometryplot_3D.py
  - ReadPinocchio5.py
  - ValidateFits.py
- src/
  - allocations.c
  - allocations.o
  - build_groups.c
  - build_groups.o
  - collapse_times.c
  - collapse_times.o
  - cosmo.c
  - cosmo.o
  - def_splines.h
  - distribute.c
  - distribute.o
  - fmax-fftw.c
  - fmax-pfft.c
  - fmax-pfft.o
  - fmax.c
  - fmax.o
  - fragment.c
  - fragment.h
  - fragment.o
  - GenIC.c
  - GenIC.o
  - initialization.c
  - initialization.o
  - LPT.c
  - LPT.o
  - Makefile
  - pinocchio.c
  - pinocchio.h
  - pinocchio.o
  - pinocchio.x
  - Pk_from_CAMB.c
  - ReadParamfile.c
  - ReadParamfile.o
  - ReadWhiteNoise.c
  - variables.c
  - variables.o
  - write_halos.c
  - write_halos.o
  - write_mass_maps.c
  - write_mass_maps.o
  - write_snapshot.c
  - write_snapshot.o
- tests/
  - ICs_piti_vs_pinocchio/
    - Cross_PK_and_Residuals_3LPT_with_vel_correction.png
    - Delta_vel_x_2LPT_with_new_ID_and_vel_correction.png
    - Delta_vel_x_3LPT.png
    - Delta_vel_x_3LPT_with_new_ID_and_vel_correction.png
    - Delta_vel_y_2LPT_with_new_ID_and_vel_correction.png
    - Delta_vel_y_3LPT.png
    - Delta_vel_y_3LPT_with_new_ID_and_vel_correction.png
    - Delta_vel_z_2LPT_with_new_ID_and_vel_correction.png
    - Delta_vel_z_3LPT.png
    - Delta_vel_z_3LPT_with_new_ID_and_vel_correction.png
    - Delta_x_2LPT.png
    - Delta_x_2LPT_with_new_ID_and_vel_correction.png
    - Delta_x_3LPT.png
    - Delta_x_3LPT_with_new_ID_and_vel_correction.png
    - Delta_y_2LPT.png
    - Delta_y_2LPT_with_new_ID_and_vel_correction.png
    - Delta_y_3LPT.png
    - Delta_y_3LPT_with_new_ID_and_vel_correction.png
    - Delta_z_2LPT.png
    - Delta_z_2LPT_with_new_ID_and_vel_correction.png
    - Delta_z_3LPT.png
    - Delta_z_3LPT_with_new_ID_and_vel_correction.png
  - only_HMF_tests/
    - MOD_GRAV_and_SCALE_DEP/
    - MOD_GRAV_and_SCALE_DEP_and_RECOMPUTE/
    - READ_PK_TABLE_and_SCALE_DEP/
    - RECOMPUTE_and_SCALE_DEP/
    - RECOMPUTE_DISPLACEMENTS_LCDM/
    - SCALE_DEP_LCDM/
  - pk_and_HMF_tests/
    - nuLCDM_READ_PK_SCALE_DEP_REC_DISPLACE_vs_old_run/
    - nuLCDM_READ_PK_SCALE_DEP_vs_old_run/
  - Readme_Pinocchio_tests_V5_1.txt
- .gitignore
- context.md
- context.py
- COPYING
- COPYRIGHT
- DOCUMENTATION
- INSTALLATION
- LICENSE
- README.md
```



----- FILE: README.md -----
```text
# Pinocchio V5.1

PINOCCHIO is a fast code to generate catalogues of cosmological dark matter halos with known mass, position, velocity and merger history. As a byproduct, it can also generate the matter density field and, upon post-processing, its lensing potential.

PINOCCHIO starts from the realisation of a linear density field on a regular grid, as in the generation of the initial conditions of a cosmological N-body simulation. It is based on Lagrangian Perturbation Theory (LPT), excursion set theory and ellipsoidal collapse. The code is made of two main parts, the computation of collapse times and LPT displacements for each particle, and the grouping of collapsed particles into haloes (`fragmentation'), with the costruction of halo merger histories and light-cone with continuous time sampling.
Collapse times are computed by Gaussian-smoothing the linear density field on many smoothing radii, then computing the second derivatives of the initial potential with FFTs; these are used to compute the collapse redshift of each particle using ellipsoidal collapse. We define the inverse collapse time as $F=1+z_{\rm c}$, and store its highest value $F_{\rm max}$ for all smoothing radii. At the final smoothing radius $R=0$ (meaning that the variance of the linear density field is only limited by the Lagrangian grid), the LPT displacement fields are computed, amounting to four vectors for each particle (of the three 3LPT displacement fields we only compute the first two, the third rotational term being negligible).
The second part of the code uses the collapse times and displacements to group particles into haloes, with an algorithm that mimics hierarchical clustering. Because collapse here is identified with orbit crossing, collapsed particles are not necessarily contained in haloes, they may be part of the filamentary network that joins haloes. So particles may be classified into uncollapsed (still in single-stream regime), filaments and halo particles. Having recognised the haloes without really running the simulation (we just performed a single 3LPT time-step when needed), we can see PINOCCHIO as a halo finder that works on the Lagrangian space of initial conditions, plus a 3LPT engine to place the haloes at the right position.

PINOCCHIO is distributed under a gnu-gpl v2 license.

[Documentation]https://github.com/pigimonaco/Pinocchio/blob/master/DOCUMENTATION
[Installation]https://github.com/pigimonaco/Pinocchio/blob/master/INSTALLATION

Github: https://github.com/pigimonaco/Pinocchio

Webpage: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

Papers:

1. The PINOCCHIO algorithm: pinpointing orbit-crossing collapsed hierarchical objects in a linear density field, Pierluigi Monaco, Tom Theuns & Giuliano Taffoni, 2002, MNRAS, 331, 587, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.MNRAS.331.587.pdf)
2. Predicting the Number, Spatial Distribution and Merging History of Dark Matter Haloes, Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Fabio Governato, Tom Quinn & Joachim Stadel, 2002, ApJ, 564, 8, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.ApJ.564.8.pdf)
3. PINOCCHIO and the Hierarchical Build-Up of Dark-Matter Haloes, Giuliano Taffoni, Pierluigi Monaco & Tom Theuns, 2002, MNRAS, 333, 623, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/taffoni.2002.MNRAS.333.623.pdf)
4. An accurate tool for the fast generation of dark matter halo catalogs, Pierluigi Monaco, Emiliano Sefusatti, Stefano Borgani, Martin Crocce, Pablo Fosalba, Ravi Sheth & Tom Theuns, 2013, MNRAS, 433, 2389, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2013.MNRAS.433.2389.pdf)
5. Improving the prediction of dark matter halo clustering with higher orders of Lagrangian Perturbation Theory, E. Munari, P. Monaco, E. Sefusatti, E. Castorina, F.G. Mohammad, S. Anselmi & S., 2017, MNRAS, 465,4658 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/munari.2017.MNRAS.465.4658.pdf)
6. Simulating cosmologies beyond LambdaCDM with Pinocchio, L.A. Rizzo, F. Villaescusa-Navarro, P. Monaco, E. Munari, S. Borgani, E. Castorina & E. Sefusatti, 2017, JCAP, 01/2017, 008 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/rizzo.2017.JCAP.0117.008.pdf)
7. Approximated methods for the generation of dark matter halo catalogs in the age of precision cosmology, P.Monaco, 2016, Galaxies, N. 4, 53, special issue: "Dark Matter: Large versus Small Scale Structures", ed. J. Gaite & A. Diaferio. [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2016.Galaxies.4.53.pdf)
8. Fast numerical method to generate halo catalogues in modified gravity (part I): second-order Lagrangian perturbation theory, C. Moretti, S. Mozzon, P. Monaco, E. Munari, M. Baldi, M. 2020, MNRAS, 493, 1153 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/moretti.2020.mnras.493.1153.pdf)
9. Euclid Preparation. Simulating thousands of Euclid spectroscopic skies, Euclid Consortium: P. Monaco, G. Parimbelli, Y. Elkhashab et al., in preparation.
```


----- FILE: src/Makefile -----
```text
# *****************************************************************
# *                        PINOCCHIO  V5.1                        *
# *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
# *****************************************************************
#
# This code was written by
# Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
# Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
# Copyright (C) 2025
#
# github: https://github.com/pigimonaco/Pinocchio
# web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# Directory settings #

RUNDIR = ./

# Executable name #

EXEC = pinocchio.x

# Compilation options #

DEBUG = NO
OMP   = YES

########## Preprocessor Definitions ##########

# Displacement LPT order #

OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# PLC reconstruction

OPTIONS += -DPLC
OPTIONS += -DMASS_MAPS
OPTIONS += -DMASS_MAPS_PROFILE
OPTIONS += -DMASS_MAPS_VALIDATE

# Dynamics of triaxial collapse #

OPTIONS += -DELL_CLASSIC
# OPTIONS += -DELL_SNG
# OPTIONS += -DTABULATED_CT

# Building groups and fragmentation #

# OPTIONS += -DCLASSIC_FRAGMENTATION

# Output #

# OPTIONS += -DSNAPSHOT
# OPTIONS += -DLIGHT_OUTPUT
# OPTIONS += -DLONGIDS
# OPTIONS += -DDOUBLE_PRECISION_PRODUCTS

# Beyond LambdaCDM models

# these are for neutrino cosmology
OPTIONS += -DSCALE_DEPENDENT
OPTIONS += -DREAD_PK_TABLE
OPTIONS += -DRECOMPUTE_DISPLACEMENTS
OPTIONS += -DREAD_HUBBLE_TABLE

# add also these for f(R) gravity
# OPTIONS += -DMOD_GRAV_FR
# OPTIONS += -DFR0=1.e-8

# Other options #

# OPTIONS += -DWHITENOISE
# OPTIONS += -DNORADIATION

########## Preprocessor Definitions end ##########

# OpenMP Configuration #

ifeq ($(OMP), YES)
OPTIONS += -DUSE_FFT_THREADS
OMP_FLAG = -fopenmp
OMP_FFTW = -lfftw3_omp -lfftw3_threads
else
OMP_FLAG =
OMP_FFTW =
endif

############## Mass Maps Library ####################

HEALPIX_ROOT ?= /usr
CFITSIO_ROOT ?= /usr

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio

########## System-specific Configuration ############

SYSTYPE ?= "mylaptop"

################### mylaptop ####################

ifeq ($(SYSTYPE),"mylaptop")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L/usr/lib/ -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/usr/include
PFFT_LIBR   = -L$(HOME)/lib -lpfft
PFFT_INCL   = -I$(HOME)/include

MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L$(HOME)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(HOME)/include/gsl

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio

endif

###################### BASE ######################

ifeq ($(SYSTYPE),"base")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(HOME)/include
PFFT_LIBR   = -L$(HOME)/lib -lpfft
PFFT_INCL   = -I$(HOME)/include

MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L$(HOME)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(HOME)/include/gsl
endif

###################### MARCONI ######################

ifeq ($(SYSTYPE),"Marconi100")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O2 -Wno-unused-result  $(OMP_FLAG)
FFTW_LIBR   = -L$(HOME)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(HOME)/include
PFFT_LIBR   = -lpfft
PFFT_INCL   =
MPI_LIBR    = -L$(SMPI_ROOT)/lib -lmpi_ibm
MPI_INCL    = -I$(SMPI_ROOT)/include
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

###################### PLEIADI ######################

ifeq ($(SYSTYPE),"pleiadi")
CC          =  mpicc
CDEBUG      = -g -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/fftw-3.3.10-bap3ootbcmiypnmo7rouk57jbdhfjuty/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/fftw-3.3.10-bap3ootbcmiypnmo7rouk57jbdhfjuty/include
PFFT_LIBR   = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/pfft-1.0.8-alpha-qy6lqerctfrhjzaoamzjhbpga2ofejqf/lib -lpfft
PFFT_INCL   = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/pfft-1.0.8-alpha-qy6lqerctfrhjzaoamzjhbpga2ofejqf/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/gsl-2.7.1-uehbowahpvubxundtx3cz4a4a5sdwzun/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/gsl-2.7.1-uehbowahpvubxundtx3cz4a4a5sdwzun/include
endif


###################### LEONARDO ######################
ifeq ($(SYSTYPE), "LeonardoBoost")

CC          =  nvc
CDEBUG      = -g -O0 -Minfo=all -v -Mneginfo $(OMP_FLAG)
COPTIMIZED  = -O3 -fast -Minfo=all -v -Mneginfo $(OMP_FLAG)

MPI         = /leonardo/prod/opt/libraries/openmpi/4.1.6/nvhpc--24.3/
MPI_LIBR    = -L$(MPI)/lib -lmpi
MPI_INCL    = -I$(MPI)/include

LIBDIR      = /leonardo_scratch/fast/CNHPC_1498509/lib/nvhpc-23.11
FFTW        = $(LIBDIR)/fftw/fftw-3.3.10
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include

GSL         = $(LIBDIR)/gsl/gsl-2.7.1
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include

PFFT        = $(LIBDIR)/pfft/pfft
PFFT_LIBR   = -L$(PFFT)/lib -lpfft
PFFT_INCL   = -I$(PFFT)/include


PMT         = /leonardo_scratch/fast/CNHPC_1498509/lib/pmt/local
PMT_LIBR    = -L$(PMT)/lib64 -lpmt
PMT_INCL    = -I$(PMT)/include

endif


########## System-specific Configuration end ##########

# Include paths for libraries #

INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)

# Library flags #

###############################################################
# HEALPix / CFITSIO configuration (override from environment)
###############################################################

# Append their include paths
INC += $(HEALPIX_INCL) $(CFITSIO_INCL)

# Append their libraries (placed after math deps)
LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(CFITSIO_LIBR) $(HEALPIX_LIBR)

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif
# Source files #

OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o write_mass_maps.o cosmo.o

# Main targets and rules

pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)

# Generic rule for compiling source files #

%.o: %.c Makefile pinocchio.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

# Clean target #

clean:
	\rm -f *.o *~ $(EXEC) run_planner


```


----- FILE: .gitignore -----
```text
*.o
.vscode
*/pinocchio.x
*json
```


----- FILE: context.py -----
```text
#!/usr/bin/env python3
import argparse, os, sys, subprocess, textwrap, unicodedata
from pathlib import Path

DEFAULT_EXCLUDE_DIRS = {
    ".git","node_modules","dist","build",".venv","venv","__pycache__",
    ".tox",".mypy_cache",".pytest_cache",".idea",".vscode",".DS_Store"
}
DEFAULT_EXCLUDE_EXT = {
    ".png",".jpg",".jpeg",".gif",".webp",".svg",".pdf",".zip",".tar",".gz",".bz2",".7z",
    ".mp3",".wav",".flac",".mp4",".mov",".avi",".mkv",
    ".exe",".dll",".so",".dylib",".bin",".class",".o",".a",".jar",
    ".otf",".ttf",".woff",".woff2",".ico",".lock"
}
DEFAULT_SECRET_BASENAMES = {
    ".env","id_rsa","id_ed25519","credentials.json","service-account.json","secret.json",
    "secrets.json","firebase-service-account.json","google-credentials.json"
}
DEFAULT_PRIORITY_FILES = [
    "README.md","README","README.rst","README.txt",
    "pyproject.toml","requirements.txt","environment.yml","Pipfile","Pipfile.lock",
    "package.json","yarn.lock","pnpm-lock.yaml","go.mod","go.sum",
    "Cargo.toml","Cargo.lock","Gemfile","pom.xml","Makefile","Dockerfile",".env.example"
]

def is_binary(path: Path, sample_bytes=8192) -> bool:
    try:
        with path.open("rb") as f:
            chunk = f.read(sample_bytes)
        if b"\x00" in chunk:
            return True
        # Heuristic: too many non-text characters
        text_ratio = sum(32 <= b <= 126 or b in (9,10,13) for b in chunk) / max(1, len(chunk))
        return text_ratio < 0.80
    except Exception:
        return True

def git_ls_files(repo: Path):
    try:
        out = subprocess.check_output(["git","ls-files"], cwd=repo, text=True)
        return [repo / p for p in out.splitlines()]
    except Exception:
        # Fallback to walk
        return [p for p in repo.rglob("*") if p.is_file()]

def within_size_limits(path: Path, max_bytes_per_file: int) -> bool:
    try:
        return path.stat().st_size <= max_bytes_per_file
    except Exception:
        return False

def looks_texty(path: Path) -> bool:
    # Allow files without extensions if not binary
    ext = path.suffix.lower()
    if ext in DEFAULT_EXCLUDE_EXT:
        return False
    return True

def normalize(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    # Strip trailing spaces to save tokens
    return "\n".join(line.rstrip() for line in s.splitlines())

def read_head(path: Path, max_lines: int) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            lines = []
            for i, line in enumerate(f, 1):
                lines.append(line)
                if i >= max_lines:
                    break
        return normalize("".join(lines))
    except Exception as e:
        return f"[Could not read file: {e}]"

def main():
    ap = argparse.ArgumentParser(description="Build a compact context.md from a repo.")
    ap.add_argument("--repo", default=".", help="Path to repo root")
    ap.add_argument("--out", default="context.md", help="Output file")
    ap.add_argument("--max-lines-per-file", type=int, default=400)
    ap.add_argument("--max-bytes-per-file", type=int, default=200_000)
    ap.add_argument("--max-total-bytes", type=int, default=5_000_000)
    ap.add_argument("--depth", type=int, default=3, help="Tree depth")
    args = ap.parse_args()

    root = Path(args.repo).resolve()
    files = git_ls_files(root)

    # Filter
    candidates = []
    for p in files:
        rel = p.relative_to(root)
        parts = set(rel.parts)
        if any(d in parts for d in DEFAULT_EXCLUDE_DIRS):
            continue
        if p.name in DEFAULT_SECRET_BASENAMES:
            continue
        if not looks_texty(p):
            continue
        if is_binary(p):
            continue
        if not within_size_limits(p, args.max_bytes_per_file):
            continue
        candidates.append(p)

    # Put priority files first
    priority = []
    others = []
    priority_names = set(DEFAULT_PRIORITY_FILES)
    for p in candidates:
        if p.name in priority_names or str(p.relative_to(root)) in DEFAULT_PRIORITY_FILES:
            priority.append(p)
        else:
            others.append(p)
    ordered = priority + sorted(others, key=lambda x: str(x).lower())

    # Build header
    def make_tree(root: Path, depth: int) -> str:
        lines = []
        def walk(base: Path, level: int):
            if level > depth:
                return
            entries = sorted([p for p in base.iterdir() if p.name not in DEFAULT_EXCLUDE_DIRS], key=lambda x: (x.is_file(), x.name.lower()))
            for e in entries:
                indent = "  " * (level - 1)
                lines.append(f"{indent}- {e.name}{'/' if e.is_dir() else ''}")
                if e.is_dir():
                    walk(e, level + 1)
        lines.append(f"{root.name}/")
        walk(root, 1)
        return "\n".join(lines)

    header = [
        "# Repository context for LLM",
        "",
        "## How this was built",
        "- Files are text only and truncated to a head limit.",
        "- Common binary and build artifacts are excluded.",
        "- Likely secret files are excluded.",
        "",
        "## Directory tree (truncated)",
        "```",
        make_tree(root, args.depth),
        "```",
        ""
    ]

    total_bytes = sum(len(s)+1 for s in header)
    chunks = ["\n".join(header)]

    for p in ordered:
        if ("/tests/" in str(p)) or ("/CAMBFiles/" in str(p)) or (".example." in str(p)) or ("HMF_Validation" in str(p)) or ("example/log" in str(p)):
            continue
        rel = str(p.relative_to(root))
        body = read_head(p, args.max_lines_per_file)
        block = f"\n\n----- FILE: {rel} -----\n```text\n{body}\n```"
        b = len(block.encode("utf-8"))
        if total_bytes + b > args.max_total_bytes:
            chunks.append("\n\n[Truncated due to max_total_bytes limit]\n")
            break
        chunks.append(block)
        total_bytes += b

    out = "\n".join(chunks)
    Path(args.out).write_text(out, encoding="utf-8")
    print(f"Wrote {args.out} ({total_bytes} bytes)")

if __name__ == "__main__":
    main()
```


----- FILE: COPYING -----
```text
		    GNU GENERAL PUBLIC LICENSE
		       Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc.
 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

			    Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Lesser General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

		    GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

			    NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

		     END OF TERMS AND CONDITIONS

	    How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
convey the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


Also add information on how to contact you by electronic and paper mail.

If the program is interactive, make it output a short notice like this
when it starts in an interactive mode:

    Gnomovision version 69, Copyright (C) year name of author
    Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, the commands you use may
be called something other than `show w' and `show c'; they could even be
mouse-clicks or menu items--whatever suits your program.

You should also get your employer (if you work as a programmer) or your
school, if any, to sign a "copyright disclaimer" for the program, if
necessary.  Here is a sample; alter the names:

  Yoyodyne, Inc., hereby disclaims all copyright interest in the program
  `Gnomovision' (which makes passes at compilers) written by James Hacker.

  <signature of Ty Coon>, 1 April 1989
  Ty Coon, President of Vice

This General Public License does not permit incorporating your program into
proprietary programs.  If your program is a subroutine library, you may
consider it more useful to permit linking proprietary applications with the
library.  If this is what you want to do, use the GNU Lesser General
Public License instead of this License.
```


----- FILE: COPYRIGHT -----
```text
Copyright (c) 2025
Pierluigi Monaco     Dipartimento di Fisica, Universita` di Trieste

This code was written by
Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
Copyright (C) 2025

github: https://github.com/pigimonaco/Pinocchio
web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


```


----- FILE: DOCUMENTATION -----
```text
*****************************************************************
*                        PINOCCHIO  V5.1                        *
*  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
*****************************************************************

This code was written by
Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
Copyright (C) 2025

github: https://github.com/pigimonaco/Pinocchio
web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


*****************************************************************

   This file contains information on how to run PINOCCHIO v5.1

*****************************************************************


This documentation is meant to give technical support in running the
PINOCCHIO code. It is not meant to explain the scientific and
cosmological background of the code. The user is supposed to be
familiar with the formalism that has been used to develop PINOCCHIO.

Complete information can be found in the papers quoted at the end of
this file, and in the PINOCCHIO web-site

http://adlibitum.oats.inaf.it/monaco/pinocchio.html


*****************************************************************


INDEX

 0. V5.1
 1. Differences with previous versions
 2. The package - V5.1
 3. The parameter file
 4. Pre-compiler directives
 5. Outputs
 6. Designing a run
 7. Scale-dependent growth rate
 8. Special behaviour
 9. PINOCCHIO papers


*****************************************************************


0. V5.1

This version is described in Euclid Collaboration: P. Monaco et al.
2025, in preparation. It has been used to produce the simulations for
the mock catalogs for the Euclid spectroscopic sample. Its main
features are:

 - it creates (1) halo catalogs at fixed time, with their mass
   function; (2) merger histories for all the halos; (3) on request,
   halo catalogs in the past light cone;
 - particle displacements can be computed with Lagrangian Perturbation
   Theory (LPT) up to third order;
 - to solve the memory limitations due to the slab distribution of
   fftw, the code is interfaced with the pfft library;
 - the identification of halos is based on an algorithm that is more
   efficient in memory usage;
 - the code allows the user to have a much better control of memory,
   so as to best exploit a given machine;
 - it can produce a "timeless snapshot", containing the inverse
   collapse times, the redshift at which the particle is accreted on a
   halo and the LPT displacement fields;
 - it can work with massive neutrinos and f(R) gravity.

The companion code Pitiless-slicer uses the output of pinocchio, in
particular the timeless snapshot, to (1) create snapshots of all
particles in a Gadget2 format, at any redshift, (2) to produce the
mass distribution of matter in shells around the observer and use them
to compute lensing convergence and shear.


*****************************************************************


1. DIFFERENCES WITH PREVIOUS VERSIONS

V1 of pinocchio is in fortran 77 and is not parallel. It adopts an
out-of-core strategy to minimize memory requirements, so it needs fast
access to the disc.

V2 has been re-written in fortran 90. The fmax code is fully parallel
but maintains the out-of-core design, so it keeps memory requirements
low but needs fast disc access.  It uses fftw-2.1.5 to compute FFTs,
in place of the code based on Numerical Recipes.  The fragment code is
provided in two versions, a scalar one and a parallel one. This
parallelization has no restrictions in validity but does not have good
scaling properties, the scalar code on a shared-memory machine is
preferable whenever possible.

V3 is fully, parallel, is written in fortran90 + C, and uses Message
Passing Interface (MPI) for communications. It has been designed to
run on hundreds if not thousands of cores of a massively parallel
super-computer. Main features:

- Parameters are passed to the code with a Gadget-like parameter file.

- The two separate codes (fmax and fragment) have been merged and no
  out-of-core strategy is adopted, so i/o is very limited but the
  amount of needed memory rises by a factor of three.  This version of
  the code needs ~110 bytes per particle.

- Fragmentation is performed by dividing the box into sub-volumes and
  distributing one sub-volume to each MPI task. The tasks do not
  communicate during fragmentation, so each sub-volume needs to extend
  the calculation to a boundary layer, where reconstruction is
  affected by boundaries. For small boxes at very high resolution the
  overhead implied by the boundary layer would become dominant; in
  this case the fragmentation code of version 2 would be preferable.

- We have merged PINOCCHIO with the generator of a Gaussian density
  field embedded in the N-GenIC code by V. Springel.

- The code has also been extended to consider a wider range of
  cosmologies including a generic, redshift-dependent equation of
  state of the quintessence.

- The code can generate a full-sky catalog of halos in the
  past-light-cone up to a redshift specified by the user. In this case
  the box is replicated using periodic boundary conditions to fill the
  needed cosmological volume.

A complete description of the code is provided in Monaco et
al. (2013); see Section 7 below.

V4 has been fully rewritten in C. It uses MPI for communications and
can use OpenMP for threading.

- Particle displacements are computed with Zeldovich, 2LPT or 3LPT
  (without the rotational term). The LPT order is decided at
  compilation time. Higher-order displacements are used to compute the
  position of a group in an output catalog, but groups can be
  reconstructed with lower-order LPT if requrired.

- The code is interfaced with fftw3, the obsolete fftw2.1.5 library is
  not adopted any more.

- Parameters are calibrated so as to reproduce the analytic mass
  function of Watson et al. (2013) in a wide range of box size and
  mass resolution.

- Memory usage has been improved by performing the fragmentation of
  collapsed particles (the procedure used to create halo catalogs) in
  slices of the original box. The number of slices is computed by
  requiring that the required memory stays below a limit, provided by
  the user in units of bytes per particle. A separate scalar code
  simulates memory requirements for a given run, so as to best design
  big runs.

- The code can output positions and velocities of all particles as a
  snapshot in Gadget2 format. Particle positions are computed as
  follows: halo particles are distributed around the center as
  Navarro-Frenk-White spheres, while uncollapsed particles and
  particles in filaments are simply displaced from their initial
  position. Alternatively, the code can write an LPT output of all
  particles, again in Gadget2 format, at the order specified at
  compilation. This can be used to generate initial conditions for a
  simulation; this option is presently limited to second-order LPT.

- In the merger history file, trees are given consecutively (before
  all halos that belong to some tree were giving together, and trees
  ought to be extracted from the complete halo list) and the
  quantities provided are slightly different from previous versions.

- The code has been extended to massive neutrinos cosmologies and to
  some modified gravity models.

V4.1.3 is the last stable version of V4, and is available on github.
It has several bug fixes, and improves in the writing of results as
Gadget snapshots. It has been used to produce some of the simulations
for Euclid.

V5 is able to scale up to very large boxes; its charcteristics are
thoroughly described in Euclid Collaboration: P. Monaco et al. 2025,
in preparation.

- To overcome the slab domain decomposition of fftw, it uses pfft that
  can distribute memory in slabs, pencils and subvolume. The
  generation of initial conditions has been completely rewritten and
  adapted to a generic domain decomposition.

- Fragmentation is performed on much larger subvolumes, where only the
  particles that may be relevant are acquired by a task. First each
  task receives only the particles that are expected to collapse by
  the final redshift, then the fragmentation is run once on a minimal
  subvolume; then each task tags all the particles that lie outside
  the domain and are near a halo that touches the border, and requests
  these particles to the other tasks. The fragmentation is then run to
  the end.

- The code can write a "timeless-snapshot" that contains all the
  inverse collapse times Fmax, the LPT displacement fields and the
  time at which the particle enters a halo.

The first two improvements lead to a much better scalability in memory
of the code.


*****************************************************************


2. THE PACKAGE - V5.1

This software has been tested successfully on several machines. You
are most warmly encouraged to use this code and report any problem to
pierluigi.monaco@inaf.it.

The src/ directory of the PINOCCHIO-5.1 package contains source files,
headers, the Makefile, the main code pinocchio.c, and these further
codes:

- run_planner.c to test whether a given configuration will achieve the
  required memory;

The scripts/ directory contains:

- ReadPinocchio5.py to read the catalogs from python;

- PlotPLCGeometry.py to visualize the tiling of boxes used in a run to
  sample the past light cone (when the PLC directive is present).

- Pinocchio2fits.py to translate pinocchio binary output to fits.

The example/ directory contains the results of an example run, see
below and the INSTALLATION file for more details.

Finally, the tests/ contains a set of tests to perform code
validation, see the documentation in that directory.


*****************************************************************


3. THE PARAMETER FILE

Parameters are passed to the code through a parameter file. The file
named "parameter_file" in the example/ directory gives a list of all
allowed parameters. Each line starting with # or % will be interpreted
as a comment. There are two kinds of keywords in this file, those that
require an argument and MUST be present, and logical keywords that
require no argument and, if absent, are assumed as FALSE. As an
exception, the PLCCenter and PLCAxis parameters must be present only
if the logical keyword PLCProvideConeData is specified. The example
parameter file is reported here for clarity. Empty lines and lines
starting with "%" or "#" are ignored.

-------------------------  parameter file  -------------------------

# This is an example parameter file for the Pinocchio 5.1 code

# run properties
RunFlag                example      % name of the run
OutputList             outputs      % name of file with required output redshifts
BoxSize                500          % physical size of the box in Mpc
BoxInH100                           % specify that the box is in Mpc/h
GridSize               128          % number of grid points per side
RandomSeed             486604       % random seed for initial conditions
% FixedIC                           % if present, the modulus in ICs is fixed to the average
% PairedIC                          % if present, the phase in ICs is shifted by PI

# cosmology
Omega0                 0.25         % Omega_0 (total matter)
OmegaLambda            0.75         % Omega_Lambda
OmegaBaryon            0.044        % Omega_b (baryonic matter)
Hubble100              0.70         % little h
Sigma8                 0.8          % sigma8; if 0, it is computed from the provided P(k)
PrimordialIndex        0.96         % n_s
DEw0                   -1.0         % w0 of parametric dark energy equation of state
DEwa                   0.0          % wa of parametric dark energy equation of state
TabulatedEoSfile       no           % equation of state of dark energy tabulated in a file
FileWithInputSpectrum  no           % P(k) tabulated in a file
                                    % "no" means that the fit of Eisenstein & Hu is used

# from N-GenIC
InputSpectrum_UnitLength_in_cm 0    % units of tabulated P(k), or 0 if it is in h/Mpc
WDM_PartMass_in_kev    0.0          % WDM cut following Bode, Ostriker & Turok (2001)

# control of memory requirements
BoundaryLayerFactor    3.0          % width of the boundary layer for fragmentation
MaxMem                 3600         % max available memory to an MPI task in Mbyte
MaxMemPerParticle      150          % max available memory in bytes per particle
PredPeakFactor         0.8          % guess for the number of peaks in the subvolume

# output
CatalogInAscii                      % catalogs are written in ascii and not in binary format
OutputInH100                        % units are in H=100 instead of the true H value
NumFiles               1            % number of files in which each catalog is written
MinHaloMass            10           % smallest halo that is given in output
AnalyticMassFunction   9            % form of analytic mass function given in the .mf.out files

# output options:
% WriteTimelessSnapshot             % writes a Gadget2 snapshot as an output
% DoNotWriteCatalogs                % skips the writing of full catalogs (including PLC)
% DoNotWriteHistories               % skips the writing of merger histories

# past light cone
StartingzForPLC        0.3          % starting (highest) redshift for the past light cone
LastzForPLC            0.0          % final (lowest) redshift for the past light cone
PLCAperture            30           % cone aperture for the past light cone
PLCProvideConeData                  % read vertex and direction of cone from paramter file
PLCCenter 0. 0. 0.                  % cone vertex in the same coordinates as the BoxSize
PLCAxis   1. 1. 0.                  % un-normalized direction of the cone axis

# Table of collapseTime file, needed if the code is compiled with TABULATED_CT
% CTtableFile  none

# CAMB PK tables, needed if the code is compiled with READ_PK_TABLE
% CAMBMatterFileTag     matterpower     % label for matter power spectrum files
% CAMBTransferFileTag   transfer_out    % label for transfer function files
% CAMBRunName           nilcdm_0.3eV    % name of CAMB run
% CAMBRedsfhitsFile     inputredshift   % list of redshifts for which the power spectrum is avail

----------------------------------------------------------------------

Outputs will be given at all the redshifts specified in a file, whose
name is given in the parameter file (OutputList). This file should
contain one redshift per line, in descending (chronological) order.
Empty lines and lines starting with "%" or "#" will be ignored.
The last redshift in this file determines the final redshift for the
run, so at least one valid line must be present.

Please remember that, because pinocchio constructs merger histories
and light cones on the fly, you will rarely need to write out a large
number of outputs.

As a default, Pinocchio uses the analytic fit of Eisenstein & Hu
(1998) to compute the power spectrum; alternatively, the user can
provide a tabulated power spectrum to the code. The filename
containing the power spectrum will be given as an argument of the
FileWithInputSpectrum keyword. The file is supposed to provide two
columns with k and P(k), in units of h/Mpc and (Mpc/h)^3; if the P(k)
is provided in other units, you can specify them using
InputSpectrum_UnitLength_in_cm. It is anyway assumed that units are
for H=100. A P(k) generated by the CAMB code can be directly provided.
The parameter Sigma8 gives the value of the sigma8 parameter, and it
is used to renormalize the power spectrum, but in case the
normalization of the power spectrum provided in a file is already
correct (say, it has been generated by CAMB), you can specify 0.0 for
Sigma8; in this case the code will not compute the normalization
constant of the power spectrum. The value of the sigma8 parameter will
be written in the standard output, so it is easy to check its
correctness.

Other options in the parameter file:

DumpProducts: after fmax is run it writes the products on dump files,
that can be read by another run to bypass the generation of ICs + fmax
and directly go to the fragmentation.

ReadProductsFromDumps: reads the dumped products and jumps to the
fragmentation.

UseTransposedFFT and MimicOldSeed: these two keywords should be added
if one wants to reproduce exactly a simulation whose initial
conditions have been produced with N-GenIC or 2lptIC.

Constrain_dim0, Constrain_dim1, Constrain_dim2: these keywords force
pfft to adopt a given domain decomposition. This can be useful to
force a pencil domain in the case pfft prefers a slab domain despite
it is less memory efficient. The product of these three numbers must
be equal to the number of MPI tasks.

ExitIfExtraParticles: if the memory overhead is not sufficient for
some tasks to store all needed particles during fragmentation, then a
warning is issued but the computation goes on. This flag forces the
computation to stop.


*****************************************************************


4. PRE-COMPILER DIRECTIVES

These are the pre-compiler switches and directives that can be
specified in the Makefile.

```


----- FILE: example/outputs -----
```text
# This file contains the list of output redshifts, in chronological
# (i.e. descending) order. The last value is the final redshift of the
# run.  The past-light cone is NOT generated using these outputs but
# is computed with continuous time sampling.

2.0
1.0
0.5
0.0

```


----- FILE: example/parameter_file -----
```text
# This is an example parameter file for the Pinocchio 5.1 code

# run properties
RunFlag                example      % name of the run
OutputList             outputs      % output list
BoxSize                1800         % physical size of the box in Mpc
BoxInH100                           % specify that the box is in Mpc/h
GridSize               128          % number of grid points per side
RandomSeed             250693       % random seed for initial conditions
% FixedIC                           % if present, the modulus in ICs is fixed to the average
% PairedIC                          % if present, the phase in ICs is shifted by PI

# cosmology
Omega0                 0.3110       % Omega_0 (CDM + Baryons)
OmegaLambda            0.6890       % Omega_Lambda
OmegaBaryon            0.0489       % Omega_b (baryonic matter)
Hubble100              0.6766       % little h
Sigma8                 0.0          % sigma8; if 0, it is computed from the provided P(k)
PrimordialIndex        0            % n_s
DEw0                   -1.0         % w0 of parametric dark energy equation of state
DEwa                   0.0          % wa of parametric dark energy equation of state
TabulatedEoSfile       no           % equation of state of dark energy tabulated in a file
FileWithInputSpectrum  CAMBTable

# from N-GenIC
InputSpectrum_UnitLength_in_cm 0    % units of tabulated P(k), or 0 if it is in h/Mpc
WDM_PartMass_in_kev    0.0          % WDM cut following Bode, Ostriker & Turok (2001)

# control of memory requirements
BoundaryLayerFactor    3.0          % width of the boundary layer for fragmentation
MaxMem                 7200         % max available memory to an MPI task in Mbyte
MaxMemPerParticle      300          % max available memory in bytes per particle
PredPeakFactor         0.8          % guess for the number of peaks in the subvolume

# output
CatalogInAscii                      % catalogs are written in ascii and not in binary format
OutputInH100                        % units are in H=100 instead of the true H value
NumFiles               1            % number of files in which each catalog is written
MinHaloMass            10           % smallest halo that is given in output
AnalyticMassFunction   9            % form of analytic mass function given in the .mf.out files

# output options:
% WriteTimelessSnapshot             % writes a Gadget2 snapshot as an output
% DoNotWriteCatalogs                % skips the writing of full catalogs (including PLC)
% DoNotWriteHistories               % skips the writing of merger histories

# past light cone
StartingzForPLC        0.3          % starting (highest) redshift for the past light cone
LastzForPLC            0.0          % final (lowest) redshift for the past light cone
PLCAperture            180          % cone aperture for the past light cone
PLCProvideConeData                  % read vertex and direction of cone from paramter file
PLCCenter 900 900 900               % cone vertex in the same coordinates as the BoxSize
PLCAxis   0. 0. 1.                  % un-normalized direction of the cone axis
NumMassPlanes          4            % number of mass planes to be constructed
MassMapNSIDE           512          % NSIDE for healpix mass maps
MassMapMasterMaxGB     2.0          % Extra memory for master storing the sheet


# Table of collapseTime file, needed if the code is compiled with TABULATED_CT
% CTtableFile  none

# CAMB PK tables, needed if the code is compiled with READ_PK_TABLE
CAMBMatterFile      CAMBFiles/pk_cb         % label for matter power spectrum files (CDM+Baryons)
CAMBRedshiftsFile   CAMBFiles/redshifts.dat % list of redshifts for the Pk table
HubbleTableFile     CAMBFiles/hubble.dat    % Hubble table file
```


----- FILE: example/PlotExample.py -----
```text
#  *****************************************************************
#  *                        PINOCCHIO  V5.1                        *
#  *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
#  *****************************************************************

#  This code was written by
#  Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
#  Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
#  Copyright (C) 2025

#  github: https://github.com/pigimonaco/Pinocchio
#  web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# This python script gives three plots of the example run provided in the distribution
# and of the same run performed as a test. It requires the run name


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv)<2:
    print("Usage: PlotExample.py [run name]")
    sys.exit(0)

runname=sys.argv[1]

(m,nm,fit)=np.loadtxt(f'pinocchio.0.0000.{runname}.mf.out',unpack=True,usecols=(0,1,5))

plt.figure()

plt.plot(m,m*nm,label='pinocchio MF',ls='-',lw=3,c='green')
plt.plot(m,m*fit,label='Watson fit',c='blue')


plt.xlim([3e13,1.e16])
plt.ylim([1.e-8,1.e-3])
plt.xscale('log')
plt.yscale('log')

plt.legend(frameon=True)
plt.title('Mass function at z=0')
plt.xlabel(r'M (M$_\odot$)',fontsize=16)
plt.ylabel(r'M n(M) (Mpc$^{-3}$)',fontsize=16)

plt.savefig('mf.png')
plt.show()

(x,y,z,m)=np.loadtxt(f'pinocchio.0.0000.{runname}.catalog.out',unpack=True,usecols=(5,6,7,11))

plt.figure()

index=(z<100)
plt.scatter(x[index],y[index],s=m[index],marker='o',c='green',label='pinocchio halos')

plt.xlim([0,500])
plt.ylim([0,500])
plt.xscale('linear')
plt.yscale('linear')

#plt.legend(frameon=True)
plt.title('Large-scale structure at z=0')
plt.xlabel(r'x (Mpc/h)',fontsize=16)
plt.ylabel(r'y (Mpc/h)',fontsize=16)

plt.savefig('lss.png')
plt.show()

(x,y,z,m)=np.loadtxt(f'pinocchio.{runname}.plc.out',unpack=True,usecols=(2,3,4,8))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

index=(m>1.e14)
ax.scatter(x[index],y[index],z[index],marker='o',c='green')

ax.set_xlabel('x (Mpc/h)')
ax.set_ylabel('y (Mpc/h)')
ax.set_zlabel('z (Mpc/h)')

plt.title('Past-light cone in comoving coordinates')

plt.savefig('plc.png')
plt.show()
```


----- FILE: INSTALLATION -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
# This code was written by
# Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
# Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
# Copyright (C) 2025
#
# github: https://github.com/pigimonaco/Pinocchio
# web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



********************************************************************

        INSTALLATION AND USAGE GUIDE FOR PINOCCHIO 5.1

********************************************************************


Pinocchio is a fully parallel hybrid MPI + OpenMP code. This guide
provides step-by-step instructions to install, compile, and run the code.


REQUIREMENTS

Before installing Pinocchio, ensure you have the following
dependencies installed:

- Message Passing Interface (MPI) libraries
- FFTW3 libraries with MPI support (http://www.fftw.org/)
- GNU Scientific Library (GSL) (http://www.gnu.org/software/gsl/)
- PFFT (Parallel FFT Library) (https://github.com/mpip/pfft)

Ensure these libraries are installed and available in your environment

INSTALLATION

- Option 1: Clone from GitHub

1. prompt> git clone https://github.com/pigimonaco/Pinocchio.git

- Option 2: Download and Extract the Tarball

1.  prompt> tar xvzf pinocchio-5.1.tar.gz

    or (if you don't have gnu tar)

    prompt> gunzip pinocchio-5.1.tar.gz ; tar xvf pinocchio-5.1.tar

NEXT STEPS

Once you have the code, follow these steps:

2.  Navigate to the source directory Pinocchio-5.1/src/

    prompt> cd Pinocchio-5.1/src/

3.  Edit the Makefile and set the values of compilation flags and
    variables to those suitable for your system. You find an example
    valid for a generic linux machine, using the gnu C compiler.

5.  Make the executable by typing

    prompt> make

    If the code does not compile, check carefully the Makefile and ask
    support to your local system manager. If everything fails, email
    all the details to monaco@oats.inaf.it.

5.  A utility, called memorytest, has been written to help design a run.
    Compile it with:

    prompt> make memorytest

    Instruction on how to use it are given in the DOCUMENTATION file.

Use the script ReadPinocchio5.py (located in the scripts directory)
for reading the binary output of the code.


RUN THE CODE

The executable is found in the src/ directory after compilation, and
is pinocchio.x (file name specified by the EXEC variable in the
Makefile). In the example/ directory you find a small working example
of a 500 Mpc/h box sampled with 200^3 particles, run with 5 MPI tasks.
It requires ~2Gb of RAM. The code has been compiled with TWO_LPT,
THREE_LPT and PLC options, as in the provided Makefile.

Open a new directory my.example/ and copy there the three files
parameter_file, outputs and PlotExample.py that can be found in the
example/ directory. The parameter file contains all needed parameters
with a quick comment, read the DOCUMENTATION for more details. The
file outputs contains a list of redshifts at which output is required.

chdir to my.directory and run the code:

> mpirun -np 5 ../src/pinocchio.x parameter_file > log

The code will generate four catalogs and four mass function files, a
past-light cone file and a cosmology file that contains basic
cosmological quantities. These files, including the log file, should
be nearly identical to the ones provided in the example/ directory.  A
quick check can be performed with the PlotExample.py python script
contained in the example/ directory.

prompt> ipython

In [1]: %run PlotExample.py

The script will produce three plots, mf.png with the mass function at
z=0, lss.png with a slice of the box and plc.png, with a 3D
visualization of the past-light cone in comoving coordinates; rotation
can be adjusted in the panel produced by python.

```


----- FILE: LICENSE -----
```text
                    GNU GENERAL PUBLIC LICENSE
                       Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc., <http://fsf.org/>
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Lesser General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

                    GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

                            NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

                     END OF TERMS AND CONDITIONS

            How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
convey the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    {description}
    Copyright (C) {year}  {fullname}

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Also add information on how to contact you by electronic and paper mail.

If the program is interactive, make it output a short notice like this
when it starts in an interactive mode:

    Gnomovision version 69, Copyright (C) year name of author
    Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, the commands you use may
be called something other than `show w' and `show c'; they could even be
mouse-clicks or menu items--whatever suits your program.

You should also get your employer (if you work as a programmer) or your
school, if any, to sign a "copyright disclaimer" for the program, if
necessary.  Here is a sample; alter the names:

  Yoyodyne, Inc., hereby disclaims all copyright interest in the program
  `Gnomovision' (which makes passes at compilers) written by James Hacker.

  {signature of Ty Coon}, 1 April 1989
  Ty Coon, President of Vice

This General Public License does not permit incorporating your program into
proprietary programs.  If your program is a subroutine library, you may
consider it more useful to permit linking proprietary applications with the
library.  If this is what you want to do, use the GNU Lesser General
Public License instead of this License.
```


----- FILE: scripts/benchmark_mass_maps.sh -----
```text
#!/usr/bin/env bash
# Automated benchmark for MASS_MAPS solver variants (Brent vs Bisection)
# Usage (from src directory):
#   ../scripts/benchmark_mass_maps.sh <paramfile> "mpirun -np 4" [--decide] [--threshold 0.03]
#
# ENV / Flags:
#   MASS_MAPS_NUM_RUNS   Number of repetitions (default 3)
#   MASS_MAPS_SOLVERS    Space separated list (default "BRENT BISECTION")
#   --decide             Emit recommendation & optional patch suggestion
#   --threshold <frac>   Relative total-time advantage needed to keep slower path (default 0.03 = 3%)
#
# Output:
#   profiling_<SOLVER>.txt  concatenated profiling sections per run
#   profiling_summary.tsv   tab-separated summary of averages
#   Recommendation printed if --decide specified.

set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <paramfile> <mpirun command (quoted)> [--decide] [--threshold 0.03]" >&2
  exit 1
fi

PARAMFILE=$1; shift
MPICMD=$1; shift
DECIDE=0
THRESH=0.03
while [ $# -gt 0 ]; do
  case "$1" in
    --decide) DECIDE=1; shift ;;
    --threshold) THRESH=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

REPS=${MASS_MAPS_NUM_RUNS:-3}
SOLVERS=${MASS_MAPS_SOLVERS:-"BRENT BISECTION"}

if [ ! -f "$PARAMFILE" ]; then
  echo "Parameter file not found: $PARAMFILE" >&2
  exit 3
fi

run_case() { # solver label
  local solver=$1
  echo "==== Solver: ${solver} ===="
  make clean >/dev/null 2>&1 || true
  if [ "$solver" = "BISECTION" ]; then
    make pinocchio MASS_MAPS_SOLVER=BISECTION >/dev/null
  else
    make pinocchio >/dev/null
  fi
  local total=0
  local outfile="profiling_${solver}.txt"
  : > "$outfile"
  for r in $(seq 1 $REPS); do
    echo "-- Run $r/$REPS" >&2
    $MPICMD ./pinocchio.x $PARAMFILE 2>&1 | tee run_${solver}_$r.log | awk '/Mass maps:/{flag=1} flag{print}' >> "$outfile" || true
  done
  echo "Captured profiling in $outfile"
  echo
}

for s in $SOLVERS; do
  run_case "$s"
done

echo "Summary (average of Total and key buckets)" | tee profiling_summary.tsv
awk '/Total \(avg\)/{t=$3} /Scan condition/ {scan=$3} /Solver condition/ {solc=$3} /Solver overhead/ {solo=$3} /Bracket:/ {br=$3} /Counts:/{print FILENAME, t, scan, solc, solo, br}' profiling_*.txt | \
  awk -v out=profiling_summary.tsv 'BEGIN{hdr="solver\tTotal\tScan\tSolCond\tSolOvh\tBracket";print hdr >> out;printf("%-18s %8s %8s %8s %8s %8s\n","solver","Total","Scan","SolCond","SolOvh","Bracket");}
  {cnt[$1]++;T[$1]+=$2;SC[$1]+=$3;SO[$1]+=$4;SH[$1]+=$5;BR[$1]+=$6;}
  END{for(f in T){at=T[f]/cnt[f];asc=SC[f]/cnt[f];aso=SO[f]/cnt[f];ash=SH[f]/cnt[f];abr=BR[f]/cnt[f];printf("%-18s %8.4f %8.4f %8.4f %8.4f %8.4f\n",f,at,asc,aso,ash,abr);
      gsub("profiling_","",f); gsub(".txt","",f);printf("%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",f,at,asc,aso,ash,abr)>>out;}}'

if [ $DECIDE -eq 1 ]; then
  # Extract totals from TSV (expect exactly two rows, else just skip decision)
  mapfile -t rows < <(tail -n +2 profiling_summary.tsv)
  if [ ${#rows[@]} -ge 2 ]; then
    # Pick the smallest total
    bestSolver=""; bestTotal=1e99
    declare -A totals
    for r in "${rows[@]}"; do
      s=$(echo "$r" | awk '{print $1}')
      tot=$(echo "$r" | awk '{print $2}')
      totals[$s]=$tot
      awk -v s=$s -v tot=$tot -v best=$bestTotal 'END{}' >/dev/null
      if awk -v a=$tot -v b=$bestTotal 'BEGIN{exit !(a<b)}'; then
        bestSolver=$s; bestTotal=$tot
      fi
    done
    echo "Decision threshold (relative): $THRESH"
    for s in "${!totals[@]}"; do
      echo "  Total[$s] = ${totals[$s]}"
    done
    # Compare best against others
    keepBoth=0
    for s in "${!totals[@]}"; do
      [ "$s" = "$bestSolver" ] && continue
      rel=$(awk -v a=${totals[$s]} -v b=$bestTotal 'BEGIN{print (a/b)-1.0}')
      echo "  Relative overhead of $s vs $bestSolver: $(printf '%.2f' $(echo $rel))"
      if awk -v r=$rel -v th=$THRESH 'BEGIN{exit !(r<th)}'; then
        : # fine
      else
        keepBoth=1
      fi
    done
    if [ $keepBoth -eq 0 ]; then
      echo "Recommendation: keep ONLY solver '$bestSolver' (others slower by > threshold)." | tee -a profiling_summary.tsv
      echo "To remove the alternative from code: delete its #ifdef path and simplify build (manual step)." | tee -a profiling_summary.tsv
    else
      echo "Recommendation: performance differences < threshold; you may drop either for simplicity (default: keep BRENT)." | tee -a profiling_summary.tsv
    fi
  else
    echo "Decision: insufficient data (need >=2 solvers)." | tee -a profiling_summary.tsv
  fi
fi

if [ $DECIDE -eq 0 ]; then
  cat <<EOF
Next steps:
- Inspect profiling_*.txt for detailed bucket breakdowns.
- (Optional) re-run with --decide to print a recommendation.
EOF
fi
```


----- FILE: scripts/HMF_validation.py -----
```text
"""
Validation Script for PINOCCHIO Simulation
------------------------------------------
This script compiles and runs a test PINOCCHIO simulation, then analyzes the Halo Mass Function (HMF).

It performs the following steps:

- Compiles the PINOCCHIO code using the Makefile with selected compiler options.
  (See below and the PINOCCHIO documentation for available flags)
- Runs the simulation using MPI-OpenMP and an input parameter file.
- Compares the output HMF with either:
    1) Theoretical fit (Watson Fit)
    2) Another mass function from a previous run.
- Saves validation and simulation logs for further analysis.


OUTPUT_DIR Structure:
The script will automatically create an output directory 'OUTPUT_DIR' as specified by the user. Inside this directory, it will:

- Copy the executable from 'MAKEFILE_DIR' to 'OUTPUT_DIR'
- Copy the parameter file from 'PARAM_SOURCE' to 'OUTPUT_DIR'
- Copy the outputs files from 'OUTPUTS_SOURCE' to 'OUTPUT_DIR'
- Modify the parameter file with the necessary settings ('RunFlag', 'BoxSize', 'GridSize').
- Perform HMF validation only at z = 0.

### Usage
Run the script using:
```bash
python HMF_validation.py
```

Dependencies:
- Python: numpy, matplotlib
- Ensure all dependencies to correctly run PINOCCHIO are installed and aligned with `SYSTYPE` in the Makefile.
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import shutil


# =========================== #
#   COMPILER and EXEC name    #
# =========================== #
EXEC_NAME    = "pinocchio.x"
COMPILER_CC  = "mpicc"
COMPILER_CPP = "mpicc++"

# =========================== #
#          USER SETTINGS      #
# =========================== #
RUN_FLAG   = "test"      # Simulation tag
BOX_SIZE   = "128"       # Box size in Mpc/h
GRID_SIZE  = "128"       # Grid size (resolution)

NUM_PROCESSES = 1  # Number of MPI processes
NUM_THREADS   = 1  # Number of OpenMP threads

# =========================== #
#        FILE LOCATIONS       #
# =========================== #
MAKEFILE_DIR   = "../src"                                        # Directory containing the Makefile
OUTPUT_DIR     = "../HMF_Validation"                             # Directory containing the final validation output
LOG_FILE       = os.path.join(OUTPUT_DIR, "VALIDATION_log.txt")  # Name of the validation log

EXEC_SOURCE    = os.path.join(MAKEFILE_DIR, EXEC_NAME)           # Executable file source
EXEC_DEST      = os.path.join(OUTPUT_DIR, EXEC_NAME)             # Executbale file destionation
PARAM_SOURCE   = "../example/parameter_file"                     # Parameter file source
PARAM_DEST     = os.path.join(OUTPUT_DIR, "parameter_file")      # Parameter file destination
OUTPUTS_SOURCE = "../example/outputs"                            # Outputs file source
OUTPUTS_DEST   = os.path.join(OUTPUT_DIR, "outputs")             # Outputs file destination


# =========================== #
#    COMPILATION OPTIONS      #
# =========================== #
COMPILATION_OPTIONS = {
    "DEBUG"        : "NO",
    "OMP"          : "YES"
}

# =========================== #
#     PREPOCESSOR FLAGS       #
# =========================== #
PREPROCESSOR_FLAGS = {

    # Displacement LPT order #
    "-DTWO_LPT"   : "YES",
    "-DTHREE_LPT" : "YES",

    # PLC reconstruction
    "-DPLC": "NO",

    # Dynamics of triaxial collapse #
    "-DELL_CLASSIC"  : "YES",
    "-DELL_SNG"      : "NO",
    "-DTABULATED_CT" : "NO",

    # Building groups and fragmentation #
    "-DCLASSIC_FRAGMENTATION": "NO",

    # Output #
    "-DSNAPSHOT"                  : "NO",
    "-DLIGHT_OUTPUT"              : "NO",
    "-DLONGIDS"                   : "NO",
    "-DDOUBLE_PRECISION_PRODUCTS" : "NO",

    # Beyond LambdaCDM models #

    # These are for neutrino cosmology
    "-DSCALE_DEPENDENT"         : "NO",
    "-DREAD_PK_TABLE"           : "NO",
    "-DONLY_MATTER_POWER"       : "NO",
    "-DRECOMPUTE_DISPLACEMENTS" : "NO",

    # Add also these for f(R) gravity #
    "-DMOD_GRAV_FR": "NO",
    "-DFR0=1.e-8"  : "NO",

    # Other options #
    "-DWHITENOISE" : "NO",
    "-DNORADIATION": "YES",

    # - This option impacts CPU code only
    # - It is used by default using the GPU with OMP
    # "-DCUSTOM_INTERPOLATION": "NO"
}

# =========================== #
#       HELPER FUNCTIONS      #
# =========================== #

def copy_files(files_to_copy):
    """
    Copies multiple files from source to destination.

    Args:
        files_to_copy (dict): A dictionary where keys are source file paths and values are destination paths.
    """
    for src, dest in files_to_copy.items():
        if not os.path.exists(src):
            print(f"Error: File not found at {src}")
            exit(1)

        shutil.copy(src, dest)
        print(f"Copied {src} -> {dest}")

def update_param_file():
    """ Updates RunFlag, BoxSize, and GridSize in the parameter file. """

    if not os.path.isfile(PARAM_DEST):
        print(f"Error: Parameter file not found at {PARAM_DEST}")
        exit(1)

    updated_lines = []
    with open(PARAM_DEST, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) > 1:
                key, value = parts[0], parts[1]
                if key == "RunFlag":
                    updated_lines.append(f"RunFlag                {RUN_FLAG}\n")
                elif key == "BoxSize":
                    updated_lines.append(f"BoxSize                {BOX_SIZE}\n")
                elif key == "GridSize":
                    updated_lines.append(f"GridSize               {GRID_SIZE}\n")
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)  # Keep empty lines unchanged

    with open(PARAM_DEST, "w") as file:
        file.writelines(updated_lines)

    print(f"Updated parameter file at {PARAM_DEST}")

def parse_parameters(param_file):
    """ Reads parameter file to extract RunFlag, BoxSize, and GridSize. """
    run_flag, box_size, grid_size = "example", "128", "128"

    with open(param_file, "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) > 1:
                if "RunFlag" in parts[0]:
                    run_flag = parts[1]
                elif "BoxSize" in parts[0]:
                    box_size = parts[1]
                elif "GridSize" in parts[0]:
                    grid_size = parts[1]

    return run_flag, box_size, grid_size

def compile_code():
    """ Compiles the PINOCCHIO code. """
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # Ensure output dir exists
    print("Compiling the code...")

    # Extract active compilation options
    compilation_str = " ".join(f"{key}={value}" for key, value in COMPILATION_OPTIONS.items())

    # Extract active preprocessor flags
    preprocessor_str = " ".join(flag for flag, enabled in PREPROCESSOR_FLAGS.items() if enabled == "YES")

    # Run `make clean`
    clean_result = subprocess.run(["make", "clean"], cwd=MAKEFILE_DIR, capture_output=True, text=True)
    if clean_result.returncode != 0:
        print("Error during 'make clean':", clean_result.stderr)
        exit(1)

    # Run `make` with dynamically selected options
    build_command = [
        "make",
        f"EXEC={EXEC_NAME}",
        f"CC={COMPILER_CC}",
        compilation_str,
        f"OPTIONS={preprocessor_str}",
        "-j16"  # Parallel compilation
    ]
    result = subprocess.run(build_command, cwd=MAKEFILE_DIR, capture_output=True, text=True)

    # Save compilation log
    with open(LOG_FILE, "w") as log:
        log.write("=== Compilation Log ===\n")
        log.write(" ".join(build_command) + "\n")
        log.write(result.stdout)
        log.write("\n")

    if result.returncode == 0:
        print("Compilation successful.")
    else:
        print("Error during compilation:", result.stderr)
        exit(1)

def run_simulation():
    """ Runs the PINOCCHIO simulation with MPI support and parameter file. """

    # Creating source-destination dictionary
    files_to_copy = {
        EXEC_SOURCE:    EXEC_DEST,
        PARAM_SOURCE:   PARAM_DEST,
        OUTPUTS_SOURCE: OUTPUTS_DEST
    }

    # Copy needed file in the OUTPUT_DIR (including the executable)
    copy_files(files_to_copy)

    # Update parameter file
    update_param_file()

    # Extract run tags from parameter file
    run_flag, box_size, grid_size = parse_parameters(PARAM_DEST)

    print("Running PINOCCHIO simulation...")

    # Set environment variable OMP_NUM_THREADS
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(NUM_THREADS)

    # Run the simulation
    run_log_file = os.path.join(OUTPUT_DIR, "log_RUN.txt")

    command = f"mpirun -np {NUM_PROCESSES} {EXEC_DEST} {PARAM_DEST} | tee {run_log_file}"
    result = subprocess.run(command, cwd=OUTPUT_DIR, shell=True, capture_output=True, env=env)

    # Save execution logs
    with open(LOG_FILE, "a") as log:
        log.write("\n=== Simulation INFO ===\n")
        log.write(f"Using RunFlag: {run_flag}, BoxSize: {box_size}, GridSize: {grid_size}")
        log.write(f"\nSimulation log saved to: {run_log_file}\n")

    if result.returncode == 0:
        print("Simulation completed successfully.")
    else:
        print("Error in simulation:", result.stderr)
        exit(1)

def plot_hmf(compare_with_fit=True, previous_run_flag=None):
    """
    Plots the Halo Mass Function (HMF) and allows the user to compare it either with:
    - The Watson Fit (default)
    - Another mass function from a previous run

    Args:
        compare_with_fit (bool): If True, compares with Watson Fit; otherwise, compares with another MF.
        previous_run_flag (str, optional): RunFlag of a previous simulation to use for comparison.
    """
    run_flag, _, _ = parse_parameters(PARAM_DEST)
    hmf_file = os.path.join(OUTPUT_DIR, f"pinocchio.0.0000.{run_flag}.mf.out")

    if not os.path.exists(hmf_file):
        print(f"Error: HMF output file {hmf_file} not found!")
        return None

    # Load data for the current run
    m, nm, fit = np.loadtxt(hmf_file, unpack=True, usecols=(0, 1, 5))

    if compare_with_fit:
        comparison_label = "Watson Fit"
        comparison_data  = fit
        HMF_plot_name    = "MHF_Validation_with_Watson_fit.png"
    else:
        while True:
            prev_hmf_file = os.path.join(OUTPUT_DIR, f"pinocchio.0.0000.{previous_run_flag}.mf.out")

            if os.path.exists(prev_hmf_file):
                break  # File exists, proceed
            else:
                print(f"Error: The file 'pinocchio.0.0000.{previous_run_flag}.mf.out' is not present. Please check and enter a correct RunFlag.")
                previous_run_flag = input("Enter the previous RunFlag for comparison: ").strip()
        _ , prev_nm, _ = np.loadtxt(prev_hmf_file, unpack=True, usecols=(0, 1, 5))
        comparison_label = f"Reference Run ({previous_run_flag})"
        comparison_data = prev_nm
        HMF_plot_name    = "MHF_Validation_with_Reference_Run.png"

    # Compute residuals
    residuals = (nm - comparison_data) / comparison_data

    # Mask out residuals that are NaN or exactly -1
    valid_residuals = residuals[~np.isnan(residuals) & (residuals != -1)]

    # Compute the average of the valid residuals
    avg_residual = np.mean(np.abs(valid_residuals))

    # Create figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    # First subplot: Mass function
    axs[0].plot(m, m * nm, label=f'PINOCCHIO ({run_flag})', ls='-', lw=2, c='crimson')
    axs[0].plot(m, m * comparison_data, label=comparison_label, ls='--', lw=2, c='dodgerblue')
    axs[0].set_ylabel(r'M n(M) (Mpc$^{-3}$)', fontsize=14)
    axs[0].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[0].set_title('Mass function at z=0')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].legend(frameon=True)

    # Second subplot: Residuals
    axs[1].plot(m, residuals, linestyle='-', c='coral', label='Residuals')
    axs[1].axhline(0, color='black', linestyle='-', lw=1)  # Reference line at 0
    axs[1].set_ylabel(r'(PINOCCHIO - Comparison)/Comparison', fontsize=14)
    axs[1].set_xlabel(r'M (M$_\odot$)', fontsize=14)
    axs[1].set_title('Residuals at z=0')
    axs[1].set_xscale('log')

    # Ensure the output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save plot in the output directory
    plot_path = os.path.join(OUTPUT_DIR, f"{HMF_plot_name}")
    plt.savefig(plot_path)
    print(f"Plot saved to {plot_path}")

    # Append results to log
    with open(LOG_FILE, "a") as log:
        log.write("\n=== Analysis Results ===\n")
        log.write(f"HMF Average Residual: {avg_residual:.5e}\n")
        log.write(f"Plot saved to: {plot_path}\n")

        log.write("\nCompilation Options:\n")
        for key, value in COMPILATION_OPTIONS.items():
            log.write(f"  {key}: {value}\n")

        log.write("\nPreprocessor Flags:\n")
        for flag, enabled in PREPROCESSOR_FLAGS.items():
            log.write(f"  {flag}: {enabled}\n")


    return avg_residual

# Main execution
def main():

    compile_code()
    run_simulation()

    # Ask the user whether to compare with the Watson Fit or a previous RunFlag
    compare_with_fit  = input("Compare the new HMF with Watson Fit? (yes/no): ").strip().lower() == 'yes'
    previous_run_flag = None

    if not compare_with_fit:
        while True:  # Keep asking until they enter something
            previous_run_flag = input("Enter the previous RunFlag for comparison: ").strip()
            if previous_run_flag:
                break
            print("Error: You must enter a valid RunFlag.")

    avg_residual = plot_hmf(compare_with_fit, previous_run_flag)

    # Print final summary
    print("\n=== Summary ===")
    print("Compilation Options:")
    for key, value in COMPILATION_OPTIONS.items():
        print(f"  {key}: {value}")

    print("\nPreprocessor Flags:")
    for flag, enabled in PREPROCESSOR_FLAGS.items():
        print(f"  {flag}: {enabled}")
```


----- FILE: scripts/outputs -----
```text
# This file contains the list of output redshifts, in chronological
# (i.e. descending) order. The last value is the final redshift of the
# run.  The past-light cone is NOT generated using these outputs but
# is computed with continuous time sampling.

2.0
1.0
0.5
0.0

```


----- FILE: scripts/parameter_file -----
```text
# This is an example parameter file for the Pinocchio 5.1 code

# run properties
RunFlag                example      % name of the run
OutputList             outputs      % output list
BoxSize                500          % physical size of the box in Mpc
BoxInH100                           % specify that the box is in Mpc/h
GridSize               64           % number of grid points per side
RandomSeed             486604       % random seed for initial conditions
% FixedIC                           % if present, the modulus in ICs is fixed to the average
% PairedIC                          % if present, the phase in ICs is shifted by PI

# cosmology
Omega0                 0.3110       % Omega_0 (CDM + Baryons)
OmegaLambda            0.6890       % Omega_Lambda
OmegaBaryon            0.0489       % Omega_b (baryonic matter)
Hubble100              0.6766       % little h
Sigma8                 0.0          % sigma8; if 0, it is computed from the provided P(k)
PrimordialIndex        0            % n_s
DEw0                   -1.0         % w0 of parametric dark energy equation of state
DEwa                   0.0          % wa of parametric dark energy equation of state
TabulatedEoSfile       no           % equation of state of dark energy tabulated in a file
FileWithInputSpectrum  CAMBTable

# from N-GenIC
InputSpectrum_UnitLength_in_cm 0    % units of tabulated P(k), or 0 if it is in h/Mpc
WDM_PartMass_in_kev    0.0          % WDM cut following Bode, Ostriker & Turok (2001)

# control of memory requirements
BoundaryLayerFactor    3.0          % width of the boundary layer for fragmentation
MaxMem                 3600         % max available memory to an MPI task in Mbyte
MaxMemPerParticle      300          % max available memory in bytes per particle
PredPeakFactor         0.8          % guess for the number of peaks in the subvolume

# output
CatalogInAscii                      % catalogs are written in ascii and not in binary format
OutputInH100                        % units are in H=100 instead of the true H value
NumFiles               1            % number of files in which each catalog is written
MinHaloMass            10           % smallest halo that is given in output
AnalyticMassFunction   9            % form of analytic mass function given in the .mf.out files

# output options:
% WriteTimelessSnapshot             % writes a Gadget2 snapshot as an output
% DoNotWriteCatalogs                % skips the writing of full catalogs (including PLC)
% DoNotWriteHistories               % skips the writing of merger histories

# past light cone
StartingzForPLC        0.3          % starting (highest) redshift for the past light cone
LastzForPLC            0.0          % final (lowest) redshift for the past light cone
PLCAperture            30           % cone aperture for the past light cone
PLCProvideConeData                  % read vertex and direction of cone from paramter file
PLCCenter 250. 250. 250.            % cone vertex in the same coordinates as the BoxSize
PLCAxis   0. 0. 1.                  % un-normalized direction of the cone axis
NumMassPlanes          4            % number of mass planes to be constructed

# Table of collapseTime file, needed if the code is compiled with TABULATED_CT
% CTtableFile  none

# CAMB PK tables, needed if the code is compiled with READ_PK_TABLE
CAMBMatterFile      CAMBFiles/pk_cb         % label for matter power spectrum files (CDM+Baryons)
CAMBRedshiftsFile   CAMBFiles/redshifts.dat % list of redshifts for the Pk table
HubbleTableFile     CAMBFiles/hubble.dat    % Hubble table file
```


----- FILE: scripts/Pinocchio2fits.py -----
```text
import numpy as np
from astropy.io import fits
import ReadPinocchio5 as rp
import sys

'''
This script converts Pinocchio binary output to fits

Usage: python3 Pinocchio2fits.py [pinocchio parameter file]
'''

VERSION='v5.0'
CATALOGS=True
PLC=True
HISTORIES=True

params={
    'RunFlag'                        :['     ','string  ','name of the run',''],
    'OutputList'                     :['     ','string  ','filename with required output redshifts',''],
    'BoxSize'                        :['     ','int     ','physical size of the box [Mpc or Mpc/h]',''],
    'BoxInH100'                      :['False','bool    ','specify that the box is in Mpc/h',''],
    'GridSize'                       :['     ','int     ','number of grid points per side',''],
    'RandomSeed'                     :['     ','int     ','random seed for initial conditions',''],
    'FixedIC'                        :['False','bool    ','modulus in ICs is fixed to the average',''],
    'PairedIC'                       :['False','bool    ','phase in ICs is shifted by PI',''],
    'Omega0'                         :['     ','float   ','Omega_0 (total matter)',''],
    'OmegaLambda'                    :['     ','float   ','Omega_Lambda',''],
    'OmegaBaryon'                    :['     ','float   ','Omega_b (baryonic matter)',''],
    'Hubble100'                      :['     ','float   ','little h',''],
    'Sigma8'                         :['     ','float   ','sigma8; if 0, it is computed from P(k)',''],
    'PrimordialIndex'                :['     ','float   ','spectral index n_s',''],
    'DEw0'                           :['     ','float   ','w0 of dark energy EoS',''],
    'DEwa'                           :['     ','float   ','wa of dark energy EoS',''],
    'TabulatedEoSfile'               :['     ','string  ','dark energy EoS tabulated in file',''],
    'FileWithInputSpectrum'          :['     ','string  ','P(k) tabulated in file',''],
    'InputSpectrum_UnitLength_in_cm' :['     ','float   ','units of tabulated P(k) [cm]',''],
    'WDM_PartMass_in_kev'            :['     ','float   ','WDM cut [keV]',''],
    'BoundaryLayerFactor'            :['     ','float   ','width of boundary layer for fragmentation',''],
    'MaxMem'                         :['     ','int     ','max available memory to an MPI task [Mbyte]',''],
    'MaxMemPerParticle'              :['     ','int     ','max available memory [bytes per particle]',''],
    'PredPeakFactor'                 :['     ','float   ','guess for the number of peaks in subvolume',''],
    'CatalogInAscii'                 :['False','bool    ','catalogs are written in ascii, not binary',''],
    'OutputInH100'                   :['False','bool    ','units are in H=100 instead of true H value',''],
    'NumFiles'                       :['     ','int     ','number of files for each catalog',''],
    'MinHaloMass'                    :['     ','int     ','smallest halo in output [particles]',''],
    'AnalyticMassFunction'           :['     ','int     ','analytic mass function given in *.mf.out',''],
    'WriteTimelessSnapshot'          :['False','bool    ','write timeless snapshot',''],
    'DoNotWriteCatalogs'             :['False','bool    ','skip writing catalogs (including PLC)',''],
    'DoNotWriteHistories'            :['False','bool    ','skip writing merger histories',''],
    'StartingzForPLC'                :['     ','float   ','starting (highest) redshift for PLC',''],
    'LastzForPLC'                    :['     ','float   ','final (lowest) redshift for PLC',''],
    'PLCAperture'                    :['     ','float   ','cone aperture for PLC [deg]',''],
    'PLCProvideConeData'             :['False','bool    ','read vertex and direction of cone',''],
    'PLCCenter'                      :['     ','3 floats','cone vertex [same coordinates as BoxSize]',''],
    'PLCAxis'                        :['     ','3 floats','un-normalized direction of cone axis',''],
    'CTtableFile'                    :['     ','string  ','filename with collapse time table to read in',''],
    'CAMBMatterFileTag'              :['     ','string  ','label for CAMB matter power spectrum files',''],
    'CAMBTransferFileTag'            :['     ','string  ','label for CAMB transfer function files',''],
    'CAMBRunName'                    :['     ','string  ','name of CAMB run',''],
    'CAMBRedshiftsFile'              :['     ','string  ','list of redshifts of CAMB power spectra','']}


def read_param_file(param_fname):


    try:
        f=open(param_fname)
    except:
        print(f'Error: parameter file {param_fname} not found')
        return 1

    print(f'Parameters in file: {param_fname}')

    for line in f:

        if line[0]=='#' or line[0]=='%':
            continue

        keys=line.split()

        if len(keys)==0 or keys[0]=='%':
            continue

        if keys[0] not in params:
            print(f"WARNING: unrecognized parameter, {keys[0]}")

        else:

            if params[keys[0]][0]=='False':
                params[keys[0]][0]='True'
            elif keys[0]=='PLCCenter' or keys[0]=='PLCAxis':
                params[keys[0]][0]=keys[1]+' '+keys[2]+' '+keys[3]
            else:
                params[keys[0]][0]=keys[1]

    f.close()

    return 0


def create_fits_file(output_file, column_data, column_names):

    # Create a FITS binary table extension
    columns = [fits.Column(name=col, format='D', array=np.array(data)) for col, data in zip(column_names, column_data)]
    table = fits.BinTableHDU.from_columns(columns)

    # Save the FITS file
    table.writeto(output_file, overwrite=True)


if len(sys.argv)<2:

    print('Usage: python3 Pinocchio2fits.py [pinocchio parameter file]')
    sys.exit(0)

param_fname=sys.argv[1]

if read_param_file(param_fname):
    sys.exit(0)

runflag=params['RunFlag'][0]
print(f'RunFlag: {runflag}')

output_fname=params['OutputList'][0]

print(f'Output list in file: {output_fname}')

try:
    outputs=np.loadtxt(output_fname,unpack=True)
except:
    print(f'Error: outputs file {output_fname} not found')
    sys.exit(0)

print(f'Output list: {outputs}')

def write_fits(fname,cat_type,cat):

    print(f'Writing fits catalog {fname}')

    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.append(('CODE',f'Pinocchio {VERSION}','https://github.com/pigimonaco/Pinocchio'))
    hdul = fits.HDUList([primary_hdu])
    hdu = fits.BinTableHDU(cat.data)
    hdu.name=cat_type
    counter=1
    hdu.header.append(('NHALOS',cat.Nhalos,'Number of halos in catalog'))
    for key in params:
        hdu.header.append((f'PAR{counter}',key,params[key][2]))
        hdu.header.append((f'VAL{counter}',params[key][0],params[key][3]))
        counter+=1
    hdul.append(hdu)
    hdul.writeto(fname, overwrite=True)

    print(f'Writing done')


def write_fits_histories(fname,cat):

    print(f'Writing fits catalog {fname}')

    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.append(('CODE',f'Pinocchio {VERSION}','https://github.com/pigimonaco/Pinocchio'))
    hdul = fits.HDUList([primary_hdu])
    data = fits.BinTableHDU(cat.data)
    myarr=np.empty((cat.Ntrees,),dtype=[('Nbranches', np.uint32), ('pointers', np.uint32)])
    myarr['Nbranches']=cat.Nbranches
    myarr['pointers']=cat.pointers
    pointers = fits.BinTableHDU(myarr)
    data.name='histories'
    pointers.name='pointers'
    counter=1
    data.header.append(('NTREES',cat.Ntrees,'number of trees'))
    data.header.append(('NBRANCH',cat.Nbranches_tot,'number of branches'))
    for key in params:
        data.header.append((f'PAR{counter}',key,params[key][2]))
        data.header.append((f'VAL{counter}',params[key][0],params[key][3]))
        counter+=1
    pointers.header.append(('NTREES',cat.Ntrees,'number of trees'))
    pointers.header.append(('NBRANCH',cat.Nbranches_tot,'number of branches'))
    hdul.append(data)
    hdul.append(pointers)
    hdul.writeto(fname, overwrite=True)

    print(f'Writing done')

#
# catalogs
#

if CATALOGS:
    for z in outputs:

        pin_fname=f'pinocchio.{z:6.4f}.{runflag}.catalog.out'
        print(f'Reading catalog {pin_fname}')

        cat=rp.catalog(pin_fname)

        print(f'The catalog contains the following fields: {cat.data.dtype.names}')

        fits_fname=pin_fname[:-3]+'fits'

        write_fits(fits_fname,'catalog',cat)


if PLC:
    pin_fname=f'pinocchio.{runflag}.plc.out'
    print(f'Reading catalog {pin_fname}')

    cat=rp.plc(pin_fname)

    print(f'The catalog contains the following fields: {cat.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    write_fits(fits_fname,'plc',cat)


if HISTORIES:
    pin_fname=f'pinocchio.{runflag}.histories.out'
    print(f'Reading catalog {pin_fname}')

    cat=rp.histories(pin_fname)

    print(f'The catalog contains the following fields: {cat.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    write_fits_histories(fits_fname,cat)




```


----- FILE: scripts/PkCamb.py -----
```text
#!/usr/bin/env python3
"""
Compute P_cb(k,z) with CAMB and (optionally) also total matter P(k) and transfer functions.

IMPORTANT CHANGE: Background tables (H(z) and optional Omega_cb(z)) are now written on a FIXED scale-factor grid
independent of the user-chosen P(k) redshift sampling. This avoids coupling background resolution to the (often
coarser) spectrum output list.

Grids:
  P(k) redshifts: controlled by --zmin/--zmax/--nz (must be >= 0 for spectra; negative z ignored for P(k)).
  Background grid: 200 linearly spaced scale factors a in [0.0001, 1.51356]. For each a we output z = 1/a - 1.
                   This spans very high redshift (z ~ 9999) down to z ~ -0.339 (a > 1 future times).

Outputs per P(k) redshift index i (always for z>=0):
    pk_*_cb_{i:03d}.dat      -> columns: k  P_cb  (prefix via --pk-prefix)

Additional outputs when --write-optional-files is set (for z>=0):
    pk_*_tot_{i:03d}.dat     -> total matter (cb + massive nu) spectrum
    tf_*_{i:03d}.dat         -> columns: k  T_cdm  T_baryon  T_photon  T_neutrino  T_nu  T_tot

Background outputs (independent grid):
    redshifts.dat            -> index  z_pk (only the P(k) list; unchanged)
    <hubble-file>            -> columns: z_bg  H(z_bg)   (fixed background grid)
    omega_cb.dat (optional)  -> columns: z_bg  Omega_cb(z_bg) (same background grid)

Unit control:
    --units h         : k in h/Mpc, P in (Mpc/h)^3 (default)
    --units physical  : k in 1/Mpc,  P in Mpc^3

Notes:
  - Power spectrum uses CAMB's interpolator var1=var2='delta_nonu' (CDM+baryons) and optionally 'delta_tot'.
  - Transfer functions: results.get_matter_transfer_data(); native k is k/h. Converted if physical units requested.
  - Transfer is log-interpolated onto the chosen k_grid.
  - CAMB is asked to compute matter power at the union of requested P(k) redshifts and the maximum background redshift
    so that H(z) is reliable up to z~9999.
"""

import argparse
import os
import numpy as np

import camb
from camb import model

def parse_redshifts(args) -> np.ndarray:
    """Return P(k) redshift list from user controls (may include negatives, but only z>=0 used for P(k))."""
    amin = 1.0 / (1 + args.zmax)
    amax = 1.0 / (1 + args.zmin)
    atab = np.linspace(amin, amax, args.nz, dtype=float)
    return 1.0 / atab - 1.0

def background_scale_factor_grid() -> tuple[np.ndarray, np.ndarray]:
    """Return (a_bg, z_bg) for the fixed background grid.

    a: 200 linearly spaced points in [1e-5, 1.5].
    z: 1/a - 1 for each of those scale factors.
    """
    a_bg = np.geomspace(1e-5, 1.5, 500, dtype=float)
    z_bg = 1.0 / a_bg - 1.0
    return a_bg, z_bg

def main():
    p = argparse.ArgumentParser(description="Compute CDM+baryons P(k) and transfer T(k) with CAMB.")
    # Redshift controls
    p.add_argument("--zmin", type=float, default=0.00, help="Minimum redshift (used with --zmax/--nz).")
    p.add_argument("--zmax", type=float, default=99.0, help="Maximum redshift (used with --zmin/--nz).")
    p.add_argument("--nz", type=int, default=100, help="Number of redshifts (used with --zmin/--zmax).")
    # Units choice
    p.add_argument("--units", choices=["h", "physical"], default="h",
                   help="Output units: 'h' => k[h/Mpc], P[(Mpc/h)^3]; 'physical' => k[1/Mpc], P[Mpc^3]. Default: h")
    # k-grid limits (interpreted according to --units)
    p.add_argument("--kmin", type=float, default=1e-4,
                   help="Minimum k in h/Mpc if --units h, or 1/Mpc if --units physical (default: 1e-4).")
    p.add_argument("--kmax", type=float, default=10.0,
                   help="Maximum k in h/Mpc if --units h, or 1/Mpc if --units physical (default: 10).")
    p.add_argument("--nk", type=int, default=512, help="Number of k points (log-spaced; default: 512).")
    # Cosmology (Planck-like defaults)
    p.add_argument("--H0", type=float, default=67.66, help="H0 in km/s/Mpc (default: 67.66).")
    p.add_argument("--Ob0", type=float, default=0.0489, help="Omega_b (default: 0.0489).")
    p.add_argument("--Ocdm0", type=float, default=0.2621, help="Omega_cdm (default: 0.2621).")
    p.add_argument("--ns", type=float, default=0.9649, help="Scalar spectral index n_s (default: 0.9649).")
    p.add_argument("--As", type=float, default=2.1e-9, help="Primordial amplitude A_s (default: 2.1e-9).")
    p.add_argument("--mnu", type=float, default=0.06, help="Sum of neutrino masses in eV (default: 0.06).")
    p.add_argument("--Neff", type=float, default=3.046, help="Effective N_eff for massless species (default: 3.046).")
    p.add_argument("--tau", type=float, default=0.054, help="Optical depth tau (default: 0.054).")
    # Linear / non-linear choice for P(k)
    p.add_argument("--nonlinear", action="store_true",
                   help="If set, include non-linear corrections for matter P(k). Linear by default.")
    # Output controls
    p.add_argument("--outdir", type=str, default=".", help="Output directory (default: current).")
    p.add_argument("--pk-prefix", type=str, default="pk", help="P(k) filename prefix (default: 'pk').")
    p.add_argument("--tf-prefix", type=str, default="tf", help="T(k) filename prefix (default: 'tf').")
    p.add_argument("--write-optional-files", action="store_true",
                   help="If set, also write total matter P(k) and transfer function files. Default: disabled (only P_cb and redshifts.dat).")
    p.add_argument("--hubble-file", type=str, default="hubble.dat",
                   help="Filename for the Hubble table (two columns: z  H(z) [km/s/Mpc]). Default: hubble.dat")

    args = p.parse_args()

    # P(k) redshift grid (user controlled)
    zs = parse_redshifts(args)
    Nz = len(zs)
    # Background grid (fixed)
    a_bg, z_bg = background_scale_factor_grid()
    os.makedirs(args.outdir, exist_ok=True)

    # Choose k grid & CAMB settings based on units
    h = args.H0 / 100.0
    if args.units == "h":
        k_grid = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)  # h/Mpc
        kmax_mpc = args.kmax * h  # convert to 1/Mpc for CAMB backend
        hubble_units = True
        k_hunit = True
    else:  # physical units
        k_grid = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)  # 1/Mpc
        kmax_mpc = args.kmax  # already 1/Mpc
        hubble_units = False
        k_hunit = False

    # CAMB parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=args.H0, ombh2=args.Ob0 * h**2, omch2=args.Ocdm0 * h**2,
                       mnu=args.mnu, nnu=args.Neff, tau=args.tau)
    pars.InitPower.set_params(As=args.As, ns=args.ns)
    pars.NonLinear = model.NonLinear_pk if args.nonlinear else model.NonLinear_none

    # CAMB redshift list: include all non-negative P(k) redshifts plus the maximum background redshift for H(z)
    pk_zs = [float(z) for z in zs if z >= 0.0]
    if len(pk_zs) == 0:
        pk_zs = [0.0]
    z_bg_max = float(np.max(z_bg))  # very high (~9999)
    if z_bg_max not in pk_zs:
        camb_redshifts = sorted(set(pk_zs + [z_bg_max]), reverse=False)
    else:
        camb_redshifts = pk_zs
    pars.set_matter_power(redshifts=camb_redshifts, kmax=kmax_mpc)  # enables transfer + P(k)

    # Run CAMB once
    results = camb.get_results(pars)

    # Interpolator for P_cb,cb with hubble_units=True and k_hunit=True
    PK_cb = results.get_matter_power_interpolator(
        nonlinear=args.nonlinear,
        var1="delta_nonu",  # CDM+baryons
        var2="delta_nonu",
        hubble_units=hubble_units,
        k_hunit=k_hunit
    )
    PK_tot = results.get_matter_power_interpolator(
        nonlinear=args.nonlinear,
        var1="delta_tot",   # CDM+baryons+massive nu
        var2="delta_tot",
        hubble_units=hubble_units,
        k_hunit=k_hunit
    )

    # Matter transfer data (linear). Provides arrays per redshift and k/h grid.
    # CAMB always stores the wavenumbers for transfer data as k/h (h/Mpc).
    mt = results.get_matter_transfer_data()
    for i, z in enumerate(zs):
        # Always write P_cb
        if z >= 0.0:
            Pk_cb_vals = PK_cb.P(z, k_grid)
            pk_cb_fname = os.path.join(args.outdir, f"{args.pk_prefix}_cb_{i:03d}.dat")
            with open(pk_cb_fname, "w") as f:
                for kk, pp in zip(k_grid, Pk_cb_vals):
                    f.write(f"{kk:.8e} {pp:.8e}\n")

        if args.write_optional_files and z >= 0.0:
            # Total matter P(k)
            Pk_tot_vals = PK_tot.P(z, k_grid)
            pk_tot_fname = os.path.join(args.outdir, f"{args.pk_prefix}_tot_{i:03d}.dat")
            with open(pk_tot_fname, "w") as f:
                for kk, pp in zip(k_grid, Pk_tot_vals):
                    f.write(f"{kk:.8e} {pp:.8e}\n")

            # Transfer functions
            zgrid = np.asarray(results.transfer_redshifts)
            iz = int(np.abs(zgrid - z).argmin())
            if abs(zgrid[iz] - z) > 1e-6:
                print(f"[warn] using nearest transfer z={zgrid[iz]:.6f} for requested z={z:.6f}")
            kh_native = mt.transfer_z("k/h", z_index=iz)
            k_native = kh_native * h if args.units == "physical" else kh_native
            T_native = [mt.transfer_z(label, z_index=iz) for label in ["delta_cdm", "delta_baryon", "delta_photon", "delta_neutrino", "delta_nu", "delta_tot"]]
            mask = k_native > 0
            Tk = [np.interp(np.log(k_grid), np.log(k_native[mask]), T[mask]) for T in T_native]
            tf_fname = os.path.join(args.outdir, f"{args.tf_prefix}_{i:03d}.dat")
            with open(tf_fname, "w") as f:
                for kk, tt0, tt1, tt2, tt3, tt4, tt5 in zip(k_grid, Tk[0], Tk[1], Tk[2], Tk[3], Tk[4], Tk[5]):
                    f.write(f"{kk:.8e} {tt0:.8e} {tt1:.8e} {tt2:.8e} {tt3:.8e} {tt4:.8e} {tt5:.8e}\n")
    # P(k) redshift index file (unchanged semantics)
    z_fname = os.path.join(args.outdir, "redshifts.dat")
    with open(z_fname, "w") as f:
        for i, zz in enumerate(zs):
            f.write(f"{i:03d} {zz:.8e}\n")

    # Write Hubble table: two columns z  H(z) [km/s/Mpc]
    hubble_path = os.path.join(args.outdir, args.hubble_file)
    with open(hubble_path, "w") as f:
        for zz in z_bg:
            Hz = results.hubble_parameter(zz)/args.H0  # km/s/Mpc
            f.write(f"{zz:.8e} {Hz:.8e}\n")

    # Optionally write Omega_cb(z) table: two columns z  Omega_cb(z)
    if args.write_optional_files:
        omega_path = os.path.join(args.outdir, "omega_cb.dat")
        with open(omega_path, "w") as f:
            for zz in z_bg:
                Ez = results.hubble_parameter(zz) / args.H0
                omega_cb = (args.Ob0 + args.Ocdm0) * (1.0 + zz) ** 3 / (Ez * Ez)
                f.write(f"{zz:.8e} {omega_cb:.8e}\n")

    extra = []
    if args.write_optional_files:
        extra.append(f"{Nz} T(k) files with prefix '{args.tf_prefix}'")
        extra.append("omega_cb.dat")
    extras = ((" and " + " and ".join(extra)) if extra else "")
    print(
        "Wrote "
        f"{Nz} P(k) files (prefix '{args.pk_prefix}'){extras}, redshifts.dat, "
        f"background H(z) grid ({args.hubble_file}) and {'omega_cb.dat, ' if args.write_optional_files else ''}"
        f"using fixed a-grid (200 points 0.0001->1.51356). Output dir: {os.path.abspath(args.outdir)} (units: {args.units})."
    )

if __name__ == "__main__":
    main()
```


----- FILE: scripts/PlcGeometryplot_3D.py -----
```text

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import argparse


########################################## Functions definition ##########################################

# Plot 3D cubes

def plot_cube( shifted_origin , side_lenght , fate , ax ):

    # Dictionary with a different colors for a different fate of the cube face

    if ( fate > 0 ):

        color={ 0: "gray", # Perché stiamo usando anche zero se c'è un if che dice > 0 ?
                1: "blue", # Starting cube
                2: "cyan", #
                3: "palegreen", #
                4: "turquoise" } #

    # Control Transparency

        alpha={ 0: 0.1,
                1: 0.3,
                2: 0.3,
                3: 0.3,
                4: 0.3 }

    # Generating vertices of a unitary cube

        base_vector = [ 0.0 , 1.0 ]
        vertices = np.array(list(product( base_vector, base_vector , base_vector ))) # Cartesian product of the base vector

    # Translation of generic vertices on the vertices of box replications that are used to sample the light cone

        for i in range(3):

            vertices[ : , i ] *= side_lenght[i]     # L is the side length of the cube
            vertices[ : , i ] += shifted_origin[i]  # C is the starting vertex of the cube (starting and replications)

    # Scatter plot of all vertices generated before

        ax.scatter3D( vertices[: , 0], vertices[: , 1] , vertices[: , 2] , s=5 )

    # Creating X , Y , Z base data point for plotting the cubes using plot_surface
    # X and Y are a 2D array of points x and y while Z is used to indicate the 2D array of heights for x and y points

        X , Y = np.meshgrid( base_vector , base_vector )
        Z = np.ones(( 2 , 2))

    # Plotting the faces of the cubes
    # We have to plot each face separately because it could potentially have a different color from the others according to the fate value

        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] , Y*side_lenght[1] + shifted_origin[1] , Z*shifted_origin[2]                    ,  alpha = alpha[fate], color = color[fate])   # Bottom Z face
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] , Y*side_lenght[1] + shifted_origin[1] , Z*(side_lenght[2] + shifted_origin[2]) ,  alpha = alpha[fate], color = color[fate])   # Top Z face
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] ,                    shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Right Y face
        ax.plot_surface(X*side_lenght[0] + shifted_origin[0] ,   side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Left Y face
        ax.plot_surface(                   shifted_origin[0] , X*side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Left X face
        ax.plot_surface(  side_lenght[0] + shifted_origin[0] , X*side_lenght[1] + shifted_origin[1] , Y*side_lenght[2] +  shifted_origin[2]  ,  alpha = alpha[fate], color = color[fate])   # Right X face

# Plot distance

def plot_dist( C , L , V , r , ax ):

    vers  = ( C + L/2 - V ) /np.linalg.norm( C + L/2 - V )

    ax.plot( np.array([ V[0] , V[0] + vers[0] * r ]),
             np.array([ V[1] , V[1] + vers[1] * r ]),
             np.array([ V[2] , V[2] + vers[2] * r ]))

# Plot the cone lateral surface and the dome

def plot_cone( starting_vertex , cone_direction , cone_aperture , ax, rstart = 0.0 , rstop = 1.0 , N_r = 5.0 , N_phi = 10.0 ):

    # Step between radius of concentric circumferences

    deltar = ( rstop - rstart) / N_r

    # Normalizing the cone axis of symmetry

    norm_direction = cone_direction/np.linalg.norm( cone_direction )

    # Find the perpendicular direction to the axis of symmetry

    x_versor = np.array([ 1. , 0. , 0. ])
    p = np.cross( norm_direction , x_versor ) # Cross product of norm_direction and x_versor in R^3 is a vector perpendicular to both norm_direction and x_versor

    # Check perpendicularity

    if np.linalg.norm(p) == 0: # i.e. p vector is parallel to norm_direction

        z_versor = np.array([ 0. , 0. , 1.])
        p = np.cross( norm_direction , z_versor )

    # Normalizing p

    p /= np.linalg.norm(p)

    # Find the perpendicular direction to p

    q = np.cross( norm_direction , p )

    # p and q form a basis for the plane perpendicular to the cone axis

    # Plot lateral surface of the cone

    AA = min( cone_aperture , np.pi ) # The aperture angle must be less than or equal to pi. If the aperture angle is less than pi, the cone is truncated at a certain height

    if cone_aperture < np.pi :

        # Loop over polar angle

        for i_phi in range( N_phi ):  # N_phi is the number of polar angles at which the generatrixes of the cone will be sampled

            phi_1 = i_phi * 2.0 * np.pi / (N_phi)   # Start value
            phi_2 = ( i_phi + 1 ) * 2.0 * np.pi / (N_phi)  # Increased value

            # Loop over the radius

            for i_r in range( N_r ): # N_r is the number of radii

                r_1 = i_r * deltar + rstart # Start value
                r_2= ( i_r + 1 ) * deltar + rstart # Increased value

                # For each combination of angles and radii, it calculates the position of two points on the surface of the cone using trigonometry and the basis vectors calculated earlier
                # It then plots a line segment connecting these two points on the given matplotlib axis object. The result is a series of line segments that together form the lateral
                # surface of the cone. The aperture angle is taken into account when determining the positions of the points, and the number of samples along the surface of the cone is
                # determined by the input parameters N_r and N_phi

                X = np.array([[ starting_vertex[0] + r_1 * np.cos(AA) * norm_direction[0] + r_1 * np.sin(AA) * ( p[0] * np.cos(phi_1) + q[0] * np.sin(phi_1) ),
                                starting_vertex[0] + r_2 * np.cos(AA) * norm_direction[0] + r_2 * np.sin(AA) * ( p[0] * np.cos(phi_1) + q[0] * np.sin(phi_1) )],

                              [ starting_vertex[0] + r_1 * np.cos(AA) * norm_direction[0] + r_1 * np.sin(AA) * ( p[0] * np.cos(phi_2) + q[0] * np.sin(phi_2) ),
                                starting_vertex[0] + r_2 * np.cos(AA) * norm_direction[0] + r_2 * np.sin(AA) * ( p[0] * np.cos(phi_2) + q[0] * np.sin(phi_2) )]])

                Y = np.array([[ starting_vertex[1] + r_1 * np.cos(AA) * norm_direction[1] + r_1 * np.sin(AA) * ( p[1] * np.cos(phi_1) + q[1] * np.sin(phi_1) ),
                                starting_vertex[1] + r_2 * np.cos(AA) * norm_direction[1] + r_2 * np.sin(AA) * ( p[1] * np.cos(phi_1) + q[1] * np.sin(phi_1) )],

                              [ starting_vertex[1] + r_1 * np.cos(AA) * norm_direction[1] + r_1 * np.sin(AA) * ( p[1] * np.cos(phi_2) + q[1] * np.sin(phi_2) ),
                                starting_vertex[1] + r_2 * np.cos(AA) * norm_direction[1] + r_2 * np.sin(AA) * ( p[1] * np.cos(phi_2) + q[1] * np.sin(phi_2) )]])

                Z = np.array([[ starting_vertex[2] + r_1 * np.cos(AA) * norm_direction[2] + r_1 * np.sin(AA) * ( p[2] * np.cos(phi_1) + q[2] * np.sin(phi_1) ),
                                starting_vertex[2] + r_2 * np.cos(AA) * norm_direction[2] + r_2 * np.sin(AA) * ( p[2] * np.cos(phi_1) + q[2] * np.sin(phi_1) )],

                              [ starting_vertex[2] + r_1 * np.cos(AA) * norm_direction[2] + r_1 * np.sin(AA) * ( p[2] * np.cos(phi_2) + q[2] * np.sin(phi_2) ),
                                starting_vertex[2] + r_2 * np.cos(AA) * norm_direction[2] + r_2 * np.sin(AA) * ( p[2] * np.cos(phi_2) + q[2] * np.sin(phi_2) )]])

                # Plot of each line segment of the lateral surface

                ax.plot_wireframe(X,Y,Z, rcount=1, ccount=1, color='yellow' , alpha=0.3)

    # Plot the dome of the cone as a 3D wireframe surface with axis along the vector norm_direction, centered at the point starting_vertex
    # The wireframe is generated by sampling points on the surface using a polar coordinate system and then converting them to Cartesian coordinates

    N_cap = max( int( 1.5 * N_phi * AA / np.pi) + 1, 5 ) # Determines the number of points along the circular cross-section of the wireframe that is being plotted

    # Specify the step sizes of the meshgrid in a single variable

    a = complex( 0 , N_phi+1 ) # real part 0 and imaginary part N_phi+1, which specifies the step size along the phi dimension
    b = complex( 0 , N_cap ) # real part 0 and imaginary part N_cap, which specifies the step size along the r dimension

    # Creates a dense grid of points in the polar coordinate system

    u , v = np.mgrid[ 0 : 2*np.pi : a, 0 : AA : b ] # u and v represent the polar coordinates of the sampled points. Meshgrid_shape = ( N_phi+1 , Ncap )
                                                    # u is the angle [0:2*np.pi]
                                                    # v is the radius (in rad)

    # Cartesian coordinates

    x = starting_vertex[0] + rstop * ( np.cos(v) * norm_direction[0] + np.sin(v) * ( p[0] * np.cos(u) + q[0] * np.sin(u) ) )
    y = starting_vertex[1] + rstop * ( np.cos(v) * norm_direction[1] + np.sin(v) * ( p[1] * np.cos(u) + q[1] * np.sin(u) ) )
    z = starting_vertex[2] + rstop * ( np.cos(v) * norm_direction[2] + np.sin(v) * ( p[2] * np.cos(u) + q[2] * np.sin(u) ) )

    # Plot the dome surface

    ax.plot_wireframe(x,y,z, rstride=5, cstride=5, color='yellow', alpha=0.5)


# Input #

parser = argparse.ArgumentParser()
parser.add_argument('--Geometry_cat', type=str, help='Input plc geometry catalog')
args = parser.parse_args()


################################# Read Input ###################################################

catalog = open( args.Geometry_cat , "r" )  # Open catalog
lines = catalog.readlines() # read lines

# Gathering the cube and PLC geometric infromation

# Cube

N_cubes = int( lines[0].split()[3] )                                 # Number of replicated cubes
Cube_side_lenght = np.array( lines[4].split()[3:] ).astype( float )  # Side lenght of the cube
IPD = float( lines[6].split()[3] )                                   # Intra particle distance

# PLC

rstart = float( lines[1].split()[3] )                                # Starting distance i.e. cone vertex
rstop = float( lines[1].split()[4] )                                 # Max distance from cone vertex

Plc_vertex = np.array( lines[2].split()[3:] ).astype( float )        # Plc vertex
Plc_direction = np.array( lines[3].split()[3:] ).astype( float )     # Plc simmetry axis

Aperture_angle = float( lines[5].split()[3] )                        # Aperture angle

# Gathering simulation output

replication_origin, rmin, rmax, fate = [], [], [], []

for i in range(8 , 8 + N_cubes) :

    # Origin of the replicated cube

    replication_origin.append( np.array(lines[i].split()[1:4]).astype(int) )

    # Min and max distance from the origin

    rmin.append( float( lines[i].split()[4] ) )
    rmax.append( float( lines[i].split()[5] ) )

    # Fate of the cube. See plot_cube for more information

    fate.append( float( lines[i].split()[6] ) )

# Closing geometry catalog

catalog.close()

# Formatting simulation output as array for plotting aims

replication_origin = np.asarray( replication_origin )
rmin = np.asarray( rmin )
rmax = np.asarray( rmax )
fate = np.asarray( fate )

# Convert the aperture angle in rad

Aperture_angle *= np.pi/180.

# Rescaling all the geometric parameters using the IPD

Plc_vertex *= IPD
Cube_side_lenght *= IPD
rstart *= IPD
rstop *= IPD
rmin *= IPD
rmax *= IPD

################################# Make Plot ###################################################

# General Plot style

plt.style.use('dark_background')
fig = plt.figure(figsize=(12, 9), dpi=800)
fig.tight_layout()

##################### First Plot ############################

# 3D view

# General settings

ax = fig.add_subplot(221, projection='3d')
ax.set_aspect("equal")
ax.grid(True)

# X-axis

ax.xaxis.pane.fill = False
ax.set_xlabel('X')
ax.xaxis.set_tick_params(labelsize=6)

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zlabel('Z')
ax.zaxis.set_tick_params(labelsize=6)

# Point of view

ax.view_init(25, -155) # First value elevation, second azimuth

# Plot cubes

for i in range( N_cubes ) :

    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Second Plot ############################

# X-Y plane projection

# General settings

ax = fig.add_subplot(222, projection='3d')
ax.set_aspect("equal")
ax.grid(True)

# X-axis

ax.xaxis.pane.fill = False
ax.set_xlabel('X')
ax.xaxis.set_tick_params(labelsize=6)

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zticks([])  # For graphical reason we set zticks to void

# Point of view == X-Y projection

ax.view_init(90, -90)

# Plot cubes

for i in range( N_cubes ) :

    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Third Plot ############################

# Y-Z plane projection

# General settings

ax = fig.add_subplot(223, projection='3d')
ax.set_aspect("equal")
ax.grid(True)

# X-axis

ax.xaxis.pane.fill = False
ax.set_xticks([]) # For graphical reason we set xticks to void

# Y-axis

ax.yaxis.pane.fill = False
ax.set_ylabel('Y')
ax.yaxis.set_tick_params(labelsize=6)

# Z-axis

ax.zaxis.pane.fill = False
ax.set_zlabel('Z')
ax.zaxis.set_tick_params(labelsize=6)

# Point of view == Y-Z projection

ax.view_init(0, 180)

# Plot cubes

for i in range( N_cubes ) :

    plot_cube( replication_origin[i] * Cube_side_lenght , Cube_side_lenght , fate[i] , ax )

# Plot cone surface

plot_cone( Plc_vertex , Plc_direction , Aperture_angle , ax , rstart=rstart , rstop = rstop , N_r = 10 , N_phi = 50 )

# Optimized scale on the axis

plt.autoscale()

##################### Fourth Plot ############################

```


----- FILE: scripts/ReadPinocchio5.py -----
```text
# HEADER

"""
Routine for reading Pinocchio's binary catalogs and PLCs.
Usage:

# CATALOGS
import ReadPinocchio5 as rp
mycat = rp.catalog("pinocchio.0.0000.example.catalog.out")
print(mycat.Mass)

# PLC
import ReadPinocchio5 as rp
myplc = rp.plc("pinocchio.example.plc.out")
print(myplc.redshift)

# HISTORIES
import ReadPinocchio5 as rp
myhist = rp.histories("pinocchio.example.histories.out")
print(myhist.name)

Written by Pierluigi Monaco, Matteo Biagetti and Emiliano Munari

LAST MODIFIED:
Pierlugi Monaco 13/3/2025

"""


import numpy as np
import os
import sys
import copy
import struct


VERBOSE=False

class catalog:

    '''
    Reads a pinocchio catalog at fixed redshift
    Usage: mycat=rp.catalog(filename, silent=False, first_file=None, last_file=None)

    filename: name of catalog file, ending with .out. The code will recognise if the catalog is writen in several files
    silent: if True, the script does not give messages in the standard output
    first_file ,last_file: if the catalog if written in several files, it allows to select a limited range of output files

    Returns:
    mycat.data: the catalog, a structured numpy array
    mycat.cat_dtype: the data type of the catalog:
        [('name', numpy.int64),
        ('Mass', numpy.float32),
        ('pos', numpy.float32, 3),
        ('vel', numpy.float32, 3),
        ('posin', numpy.float32, 3),
        ('npart', numpy.int32)]       (absent in light format)
    mycat.Nhalos: number of halos in the catalog
    mycat.Nfiles:  number of files in which the catalog is written

    mycat.Mass: Mass field (for faster access)
    mycat.pos: pos field
    mycat.vel: vel field
    mycat.posin: posin field
    mycat.npart: npart field            (absent in light format)

    '''

    def __init__(self,filename,silent=False,first_file=None,last_file=None):

        if VERBOSE:
            silent=False

        # checks that the filename contains 'catalog'
        if not 'catalog' in filename:
            print("Are you sure you are providing the right file name?")
            if 'plc' in filename:
                print("...this looks like a plc file, please read it with rp.plc(filename)")
            elif 'histories' in filename:
                print("...this looks like a histories file, please read it with rp.histories(filename)")
            return None

        # checks that the input file ends by ".out"
        last_ext=filename.rfind('.')
        if filename[last_ext:]!='.out':
            print("The catalog file should end with .out, the file number extension will be checked by the code")
            return None

        # checks that the file exists, or that there are multiple files, and in case count them
        if not os.path.exists(filename):

            if not os.path.exists(filename+'.0'):

                print("file {} or {} not found:".format(filename,filename+'.0'))
                return None

            else:

                Nfiles=1
                while os.path.exists(filename+'.{}'.format(Nfiles)):
                    Nfiles+=1
                if not silent:
                    print("The catalog is written in {} files".format(Nfiles))

        else:

            Nfiles=1
            if not silent:
                print("The catalog is written in 1 file")

        # opens the (first) file and reads the record length
        if Nfiles==1:
            if not silent:
                print('reading header of file '+filename)
            reading = np.fromfile(filename,dtype=np.int32,count=10)
        else:
            if not silent:
                print('reading header of file '+filename+'.0')
            reading = np.fromfile(filename+'.0',dtype=np.int32,count=10)

        # number of tasks that write into a single file
        NTasksPerFile = reading[1]
        if not silent:
            print('This file has been written by {} tasks'.format(NTasksPerFile))

        # the header gives either the number of slices (always < 10) or the record length
        # in the first case the record length can be read from the 8th integer
        if reading[2]>10:
            newRun=True
            record_length=reading[2]
            Nslices=1
            if not silent:
                print('This is new output format, record length: {}'.format(record_length))
        else:
            newRun=False
            record_length=reading[7]
            Nslices=reading[2]
            if not silent:
                print('This is classic output format, record length: {}'.format(record_length))
                if Nslices==1:
                    print('The box has been fragmented in 1 slice')
                else:
                    print('The box has been fragmented in {} slices'.format(Nslices))


        # sets the record
        if record_length==96:          # this is the classic format in double precision

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float64),
                             ('posin', np.float64, 3),
                             ('pos',   np.float64, 3),
                             ('vel',   np.float64, 3),
                             ('npart', np.int32) ]

            stored_dtype  =[ ('fort',  np.int32),
                             ('name',  np.int64),
                             ('Mass',  np.float64),
                             ('posin', np.float64, 3),
                             ('pos',   np.float64, 3),
                             ('vel',   np.float64, 3),
                             ('npart', np.int32),
                             ('pad'  , np.int32),
                             ('trof',  np.int32) ]


        elif record_length==56:

            if newRun:                     # the new format has posin in a different position

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('posin', np.float32, 3),
                                 ('npart', np.int32) ]

                stored_dtype  =[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('posin', np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32) ]

            else:                          # this is the classic format in single precision

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('posin', np.float32, 3),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('npart', np.int32) ]

                stored_dtype  =[ ('fort',  np.int32),
                                 ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('posin', np.float32, 3),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32),
                                 ('trof',  np.int32) ]

        elif record_length==48:        # this is the new light format

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('posin', np.float32, 3) ]

            stored_dtype = self.cat_dtype


        elif record_length==40:        # this was used in NewClusterMocks

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('npart', np.int32) ]

            stored_dtype  =[ ('fort',  np.int32),
                             ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('npart', np.int32),
                             ('trof',  np.int32) ]

        else:
            print("sorry, I do not recognize this record length")
            return None

        # decides what files to read
        if Nfiles>1:
            if first_file is None:
                first_file=0
            elif first_file<0:
                first_file=0
            elif first_file>Nfiles:
                first_file=Nfiles
            if last_file is None:
                last_file=Nfiles
            else:
                last_file += 1
                if last_file<first_file:
                    last_file=first_file+1
                elif last_file>Nfiles:
                    last_file=Nfiles
            # this is to be used to define a pythonic range
            if not silent:
                print("I will read files in the python range from {} to {}".format(first_file,last_file))
        else:
            first_file=0
            last_file=Nfiles
        self.Nfiles=Nfiles

        # prepares to read the file(s)
        self.data=None
        NhalosPerFile=np.zeros(Nfiles,dtype=np.int64)

        # loops on the files to be read
        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = filename+'.{}'.format(myfile)

            if VERBOSE:
                print("reading file {}".format(myfname))

            # reads the file as a binary object
            with open(myfname,'rb') as f:
                bindata=f.read()
                f.close()
            filesize = len(bindata)

            # reads the number of halos contained in each header
            Nblocks=NTasksPerFile * Nslices     # this is the number of blocks to read
            Nwritten=0                          # some tasks may have no halos to write
            pos=16
            vechalo=[]
            while pos < filesize:
                # reads the number of halos to read
                vec=struct.Struct('iii').unpack(bindata[pos:pos+12])[1]
                pos+=12
                vechalo.append(vec)
                if vec>0:
                    if newRun:
                        pos += 8+vec*record_length
                    else:
                        pos += vec*(record_length+8)

            vechalo = np.asarray(vechalo)
            #print(vechalo.size, Nblocks, vechalo)
            NhalosPerFile[myfile]=vechalo.sum()
            Nwritten=(vechalo>0).sum()

            if VERBOSE:
                print(f"Found {Nwritten} non-void blocks in this file over {Nblocks}")

            # checks that the lenght of data is as expected
            if newRun:
                file_length = 16 + Nblocks*12 + Nwritten*8 + NhalosPerFile[myfile]*record_length
            else:
                file_length = 16 + Nblocks*12 + NhalosPerFile[myfile]*(record_length+8)

            if file_length != filesize:
                print(f'ERROR: inconsistency in the file size, should be {file_length} but I found {filesize} bytes')
                return None
            elif VERBOSE:
                print(f'predicted file length {file_length} matches with file size {filesize}')

            # format to select information from the byte object
            if newRun:
                cleanForm='16x '
                for block in range(Nblocks):
                    if vechalo[block]>0:
                        cleanForm+='16x {}s 4x '.format(vechalo[block]*record_length)
                    else:
                        cleanForm+='12x '

            else:
                cleanForm='16x '
                for block in range(Nblocks):
                    if vechalo[block]>0:
                        cleanForm+='12x {}s '.format(vechalo[block]*(record_length+8))
                    else:
                        cleanForm+='12x '

            # removes all unwanted information from the binary structure
            try:
                cleaned = b''.join(struct.unpack(cleanForm, bindata))
            except:
                print("ERROR: I do not recognise the data structure!")
                return None
            del bindata

            # reads the catalog from the cleaned bynary structure
            thiscat = np.frombuffer(cleaned, dtype=stored_dtype)
            del cleaned

            # removes unwanted columns from the catalog
            if self.data is None:
                self.data = np.zeros(NhalosPerFile[myfile], dtype=self.cat_dtype)
            else:
                self.data.resize(self.data.shape[0]+NhalosPerFile[myfile])

            for name in self.data.dtype.names:
                self.data[name][-NhalosPerFile[myfile]:]=thiscat[name]
            del thiscat

            if not silent:
                print("done with file {}".format(myfname))

        if not silent:
            print("Reading catalog done, {} groups found".format(NhalosPerFile.sum()))

        self.Nhalos = len(self.data)

        # Create few pointers to make it compatible with previous versions
        self.Mass  = self.data['Mass']
        self.pos   = self.data['pos']
        if record_length>48:
            self.Npart = self.data['npart']
        self.vel   = self.data['vel']

class plc:

    '''
    Reads a pinocchio catalog on the past light cone
    Usage: myplc=rp.plc(filename, silent=False, first_file=None, last_file=None, onlyNfiles=False)

    filename: name of catalog file, ending with .out. The code will recognise if the catalog is writen in several files
    silent: if True, the script does not give messages in the standard output
    first_file ,last_file: if the catalog if written in several files, it allows to select a limited range of output files
    onlyNfiles: it computes Nfiles and exits returning its value, without reading the catalog

    Returns:
    myplc.data: the catalog, a structured numpy array
    myplc.cat_dtype: the data type of the catalog:
        [('name', numpy.uint64),
        ('truez', numpy.float32),
        ('pos', numpy.float32, 3),       (absent in light format)
        ('vel', numpy.float32, 3),       (absent in light format)
        ('Mass', numpy.float32),
        ('theta', numpy.float32),
        ('phi', numpy.float32),
        ('vlos', numpy.float32),         (absent in light format)
        ('obsz', numpy.float32)]
    myplc.Nhalos: number of halos in the catalog
    myplc.Nfiles: number of files in which the catalog is written

    '''

    def __init__(self,filename,silent=False,first_file=None,last_file=None,onlyNfiles=False):

```


----- FILE: scripts/ValidateFits.py -----
```text
import numpy as np
from astropy.io import fits
import ReadPinocchio5 as rp
import sys

if len(sys.argv)<2:

    print('Usage: python3 validate_fits.py [pinocchio parameter file]')
    sys.exit(0)

VERSION='v5.0'
CATALOGS=True
PLC=True
HISTORIES=True

params={
    'RunFlag'                        :['     ','string  ','name of the run','{}'],
    'OutputList'                     :['     ','string  ','filename with required output redshifts','{}'],
    'BoxSize'                        :['     ','int     ','physical size of the box [Mpc or Mpc/h]','{}'],
    'BoxInH100'                      :['False','bool    ','specify that the box is in Mpc/h','{}'],
    'GridSize'                       :['     ','int     ','number of grid points per side','{}'],
    'RandomSeed'                     :['     ','int     ','random seed for initial conditions','{}'],
    'FixedIC'                        :['False','bool    ','modulus in ICs is fixed to the average','{}'],
    'PairedIC'                       :['False','bool    ','phase in ICs is shifted by PI','{}'],
    'Omega0'                         :['     ','float   ','Omega_0 (total matter)','{}'],
    'OmegaLambda'                    :['     ','float   ','Omega_Lambda','{}'],
    'OmegaBaryon'                    :['     ','float   ','Omega_b (baryonic matter)','{}'],
    'Hubble100'                      :['     ','float   ','little h','{}'],
    'Sigma8'                         :['     ','float   ','sigma8; if 0, it is computed from P(k)','{}'],
    'PrimordialIndex'                :['     ','float   ','spectral index n_s','{}'],
    'DEw0'                           :['     ','float   ','w0 of dark energy EoS','{}'],
    'DEwa'                           :['     ','float   ','wa of dark energy EoS','{}'],
    'TabulatedEoSfile'               :['     ','string  ','dark energy EoS tabulated in file','{}'],
    'FileWithInputSpectrum'          :['     ','string  ','P(k) tabulated in file','{}'],
    'InputSpectrum_UnitLength_in_cm' :['     ','float   ','units of tabulated P(k) [cm]','{}'],
    'WDM_PartMass_in_kev'            :['     ','float   ','WDM cut [keV]','{}'],
    'BoundaryLayerFactor'            :['     ','float   ','width of boundary layer for fragmentation','{}'],
    'MaxMem'                         :['     ','int     ','max available memory to an MPI task [Mbyte]','{}'],
    'MaxMemPerParticle'              :['     ','int     ','max available memory [bytes per particle]','{}'],
    'PredPeakFactor'                 :['     ','float   ','guess for the number of peaks in subvolume','{}'],
    'CatalogInAscii'                 :['False','bool    ','catalogs are written in ascii, not binary','{}'],
    'OutputInH100'                   :['False','bool    ','units are in H=100 instead of true H value','{}'],
    'NumFiles'                       :['     ','int     ','number of files for each catalog','{}'],
    'MinHaloMass'                    :['     ','int     ','smallest halo in output [particles]','{}'],
    'AnalyticMassFunction'           :['     ','int     ','analytic mass function given in *.mf.out','{}'],
    'WriteTimelessSnapshot'          :['False','bool    ','write timeless snapshot','{}'],
    'DoNotWriteCatalogs'             :['False','bool    ','skip writing catalogs (including PLC)','{}'],
    'DoNotWriteHistories'            :['False','bool    ','skip writing merger histories','{}'],
    'StartingzForPLC'                :['     ','float   ','starting (highest) redshift for PLC','{}'],
    'LastzForPLC'                    :['     ','float   ','final (lowest) redshift for PLC','{}'],
    'PLCAperture'                    :['     ','float   ','cone aperture for PLC [deg]','{}'],
    'PLCProvideConeData'             :['False','bool    ','read vertex and direction of cone','{}'],
    'PLCCenter'                      :['     ','3 floats','cone vertex [same coordinates as BoxSize]','{}'],
    'PLCAxis'                        :['     ','3 floats','un-normalized direction of cone axis','{}'],
    'CTtableFile'                    :['     ','string  ','filename with collapse time table to read in','{}'],
    'CAMBMatterFileTag'              :['     ','string  ','label for CAMB matter power spectrum files','{}'],
    'CAMBTransferFileTag'            :['     ','string  ','label for CAMB transfer function files','{}'],
    'CAMBRunName'                    :['     ','string  ','name of CAMB run','{}'],
    'CAMBRedshiftsFile'              :['     ','string  ','list of redshifts of CAMB power spectra','{}']}


def read_param_file(param_fname):


    try:
        f=open(param_fname)
    except:
        print(f'Error: parameter file {param_fname} not found')
        return 1

    #print(f'Parameters in file: {param_fname}')

    for line in f:

        if line[0]=='#' or line[0]=='%':
            continue

        keys=line.split()

        if len(keys)==0 or keys[0]=='%':
            continue

        if keys[0] not in params:
            print(f"WARNING: unrecognized parameter, {keys[0]}")

        else:

            if params[keys[0]][0]=='False':
                params[keys[0]][0]='True'
            elif keys[0]=='PLCCenter' or keys[0]=='PLCAxis':
                params[keys[0]][0]=keys[1]+' '+keys[2]+' '+keys[3]
            else:
                params[keys[0]][0]=keys[1]

    f.close()

    return 0

if len(sys.argv)<2:

    print('Usage: python3 validate_fits.py [pinocchio parameter file]')
    sys.exit(0)

param_fname=sys.argv[1]

if read_param_file(param_fname):
    sys.exit(0)

runflag=params['RunFlag'][0]
print(f'RunFlag: {runflag}')

output_fname=params['OutputList'][0]

print(f'Output list in file: {output_fname}')

try:
    outputs=np.loadtxt(output_fname,unpack=True)
except:
    print(f'Error: outputs file {output_fname} not found')
    sys.exit(0)

print(f'Output list: {outputs}')

Nerrors=0

if CATALOGS:
    for z in outputs:

        pin_fname=f'pinocchio.{z:6.4f}.{runflag}.catalog.out'
        print(f'Reading catalog {pin_fname}')

        cat1=rp.catalog(pin_fname,silent=True)

        print(f'The pinocchio catalog contains the following fields: {cat1.data.dtype.names}')

        fits_fname=pin_fname[:-3]+'fits'

        f=fits.open(fits_fname)

        names=[]
        for i in range(1,10):
            n=f'TTYPE{i}'
            if n in f[1].header:
                names.append(f[1].header[n])
            else:
                break
        print(f'The fits catalog contains the following fields: {names}')

        print('number of halos:',cat1.Nhalos,int(f[1].header['NHALOS']))
        if cat1.Nhalos != int(f[1].header['NHALOS']):
            Nerrors+=1
            print ("ERROR: number of halos does not match!")

        for field in names:
            if ~(cat1.data[field]==f[1].data[field]).sum() >0:
                Nerrors+=1
                print(f'error in {field}')

if PLC:
    pin_fname=f'pinocchio.{runflag}.plc.out'
    print(f'Reading catalog {pin_fname}')

    cat1=rp.plc(pin_fname)

    print(f'The catalog contains the following fields: {cat1.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    f=fits.open(fits_fname)

    names=[]
    for i in range(1,10):
        n=f'TTYPE{i}'
        if n in f[1].header:
            names.append(f[1].header[n])
        else:
            break
    print(f'The fits catalog contains the following fields: {names}')

    print('number of halos:',cat1.Nhalos,int(f[1].header['NHALOS']))
    if cat1.Nhalos != int(f[1].header['NHALOS']):
        Nerrors+=1
        print ("ERROR: number of halos does not match!")

    for field in names:
        if ~(cat1.data[field]==f[1].data[field]).sum() >0:
            Nerrors+=1
            print(f'error in {field}')


if HISTORIES:
    pin_fname=f'pinocchio.{runflag}.histories.out'
    print(f'Reading catalog {pin_fname}')

    cat1=rp.histories(pin_fname)

    print(f'The catalog contains the following fields: {cat1.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    f=fits.open(fits_fname)

    names=[]
    for i in range(1,10):
        n=f'TTYPE{i}'
        if n in f[1].header:
            names.append(f[1].header[n])
        else:
            break
    print(f'The fits catalog contains the following fields: {names}')

    print('number of trees:',cat1.Ntrees,int(f[1].header['NTREES']))
    if cat1.Ntrees != int(f[1].header['NTREES']):
        Nerrors+=1
        print ("ERROR: number of trees does not match!")

    for field in names:
        if ~(cat1.data[field]==f[1].data[field]).sum() >0:
            Nerrors+=1
            print(f'error in {field}')

    if ~(cat1.Nbranches==f[2].data['Nbranches']).sum() >0:
        Nerrors+=1
        print(f'error in Nbranches')

    if ~(cat1.pointers==f[2].data['pointers']).sum() >0:
        Nerrors+=1
        print(f'error in pointers')



if Nerrors==0:
    print('*********')
    print('no errors')
    print('*********')


```


----- FILE: src/allocations.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

// #define VERBOSE

#define ALIGN_MEMORY_BLOCK(M) \
  while ((M) % ALIGN)         \
  M++

int organize_main_memory()
{

  /* here it computes the size of required memory
     and returns the number of allocatable particles for fragmentation

     ZELDOVICH DISPLACEMENTS:

     name                sizeof                    N

     products            4 float + 1 int = 20      MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size

     first_derivatives and density will share the same memory as second_derivatives

     these will be allocated by pfft routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     This is needed by GenIC:
     seedtable           1 int = 4                 MyGrids[*].GSglobal[_x_] * MyGrids[*].GSglobal[_y_]

     Total for ZEL: 76 + 16 = 92 + overhead for fragmentation + smaller things (seedtable etc)

     2LPT/3LPT DISPLACEMENTS

     name                sizeof                    N

     products (2LPT):    7 float + 1 int = 32      MyGrids[0].total_local_size
     products (3LPT)    13 float + 1 int = 56     MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size
     kvector_2LPT        1 double = 8              MyGrids[0].total_local_size_fft
     kvector_3LPT_1      1 double = 8              MyGrids[0].total_local_size_fft
     kvector_3LPT_2      1 double = 8              MyGrids[0].total_local_size_fft

     first_derivatives and density will share the same memory as second_derivatives
     source_?LPT* will share the same memory as kvector_?LPT*

     This is needed by GenIC:
     seedtable           1 int = 4                 MyGrids[*].GSglobal[_x_] * MyGrids[*].GSglobal[_y_]

     These will be allocated by pfft routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     Total for 2LPT: 32 + 80 + 16 = 128 + overhead for fragmentation + smaller things
     Total for 3LPT: 56 + 80 + 16 = 152 + overhead for fragmentation + smaller things

     Timeless snapshot adds 4 bytes to the count.

     In case or recomputation of displacements:

     products (ZEL):     7 float + 1 int = 32      MyGrids[0].total_local_size
     products (2LPT):    13 float + 1 int = 56     MyGrids[0].total_local_size
     products (3LPT):    25 float + 1 int = 104    MyGrids[0].total_local_size
     Total for ZEL:   32 + 80 + 16 = 128 + overhead for fragmentation
     Total for 2LPT:  56 + 80 + 16 = 152 + overhead for fragmentation
     Total for 3LPT: 104 + 80 + 16 = 200 + overhead for fragmentation

     Products for fragmentation will require the same memory as those for fmax
     plus the overhead for boundary layers.

     Group catalogs will require less memory than the fmax products. This will be
     checked with an estimate of the number of halos.

     memory.prods:            fmax products in fft space
     memory.fields_to_keep    LPT sources in kspace
     memory.fields            second_derivatives, // seedtable
     memory.first_allocated:  the three above
     memory.fft:              needed by pfft vectors
     memory.fmax_total:       needed by fmax, the sum of fft and first_allocated
     memory.frag_prods:       fmax products in subbox space
     memory.frag_arrays:      linking list and group ID
     memory.groups:           groups catalog, including histories, map and PLC
     memory.frag_allocated:   needed by fragment, allocated
     memory.frag_total:       needed by fragment, including FFTs
     memory.all_allocated:    total allocated memory (not including FFTs)
     memory.all:              total memory needed (including FFTs)

   */

  int igrid;
  unsigned int myNalloc, Nalloc;

  /***********************************************/
  /* COLLAPSE TIMES PART                         */
  /* this is the memory needed to store products */
  memory.prods = MyGrids[0].total_local_size * sizeof(product_data); /* products */
  ALIGN_MEMORY_BLOCK(memory.prods);
  memory.fields = memory.fields_to_keep = 0;

  /* these fields are needed after products if displacements are recomputed */
  for (igrid = 0; igrid < Ngrids; igrid++)
  {
    memory.fields_to_keep += MyGrids[igrid].total_local_size_fft * sizeof(double); /* kdensity */
    ALIGN_MEMORY_BLOCK(memory.fields_to_keep);
  }
#ifdef TWO_LPT
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double); /* kvector_2LPT */
  ALIGN_MEMORY_BLOCK(memory.fields_to_keep);
#ifdef THREE_LPT
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double); /* kvector_3LPT_1 */
  ALIGN_MEMORY_BLOCK(memory.fields_to_keep);
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double); /* kvector_3LPT_2 */
  ALIGN_MEMORY_BLOCK(memory.fields_to_keep);
#endif
#endif

  /* these fields are needed only by Fmax */
  for (igrid = 0; igrid < Ngrids; igrid++)
    for (int dim = 0; dim < 6; dim++)
    {
      memory.fields += MyGrids[igrid].total_local_size * sizeof(double); /* second_derivatives */
      ALIGN_MEMORY_BLOCK(memory.fields);
    }

  /* for (igrid=0; igrid<Ngrids; igrid++) */
  /*   memory.fields += MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);  /\* seedtable *\/ */

  /* memory allocated by PFFT */
  for (igrid = 0, memory.fft = 0; igrid < Ngrids; igrid++)
  {
    memory.fft += MyGrids[igrid].total_local_size_fft * sizeof(double);
    ALIGN_MEMORY_BLOCK(memory.fft);
    memory.fft += MyGrids[igrid].total_local_size_fft * sizeof(double);
    ALIGN_MEMORY_BLOCK(memory.fft);
  }

  memory.first_allocated = memory.prods + memory.fields_to_keep + memory.fields;
  memory.fmax_total = memory.first_allocated + memory.fft;

  /********************************************************************************/
  /* FRAGMENTATION PART                                                           */
  /* this is the memory needed to fragment the collapsed medium, including PLC data
     and the map for needed particles on the boundary */
  memory.groups = subbox.PredNpeaks * sizeof(group_data);
  ALIGN_MEMORY_BLOCK(memory.groups);
  memory.groups += subbox.PredNpeaks * +sizeof(histories_data);
  ALIGN_MEMORY_BLOCK(memory.groups);
#ifdef PLC
  memory.groups += plc.Nmax * sizeof(plcgroup_data);
  ALIGN_MEMORY_BLOCK(memory.groups);
#endif
  memory.groups += subbox.maplength * sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK(memory.groups);
  memory.groups += subbox.maplength * sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK(memory.groups);

  /* Here it computes the number of particles it can allocate for the subbox */
  size_t other_mem = memory.prods + memory.groups
#ifdef RECOMPUTE_DISPLACEMENTS
                     + memory.fields_to_keep + memory.fft
#endif
      ;

  if (other_mem > MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle)
    myNalloc = 0;
  else
    /* -10 is to compensate for remainder in integer division */
    myNalloc = (MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle - other_mem) / (sizeof(product_data) + FRAGFIELDS * sizeof(int)) - 10;

  /* Nalloc will be the smallest among all tasks */
  MPI_Reduce(&myNalloc, &Nalloc, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nalloc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  memory.frag_prods = Nalloc * sizeof(product_data);
  memory.frag_arrays = Nalloc * FRAGFIELDS * sizeof(int);
  /* this is the memory occupied by fragmentation */
  memory.frag_allocated = memory.prods + memory.frag_prods + memory.groups + memory.frag_arrays
#ifdef RECOMPUTE_DISPLACEMENTS
                          + memory.fields_to_keep
#endif
      ;
  memory.frag_total = memory.frag_allocated
#ifdef RECOMPUTE_DISPLACEMENTS
                      + memory.fft
#endif
      ;

#ifdef CLASSIC_FRAGMENTATION
  if (myNalloc < subbox.Npart)
  {
    printf("ERROR: Task %d can allocate only %d subbox particles, while %d are needed;\n",
           ThisTask, myNalloc, subbox.Npart);
    printf("       a large overhead is probably needed, please increase MaxMemPerParticle\n");
    return (size_t)0;
  }
#endif

  /* this is the largest amount of memory needed by the code */
  memory.all_allocated = (memory.first_allocated > memory.frag_allocated ? memory.first_allocated : memory.frag_allocated);
  memory.all = (memory.fmax_total > memory.frag_total ? memory.fmax_total : memory.frag_total);

  double myfraction = (double)Nalloc / (double)subbox.Ngood;
  double minfraction;
  MPI_Reduce(&myfraction, &minfraction, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if (!ThisTask)
  {
    if (minfraction < 0.5)
    {
      printf("ERROR: for some tasks the number of allocatable particles is only a factor %8f\n", minfraction);
      printf("       times the number of good particles in the subbox; this is far too small!\n");
      printf("       Please increase MaxMemPerParticle\n");
      fflush(stdout);
      return (size_t)0;
    }
  }

  return Nalloc;
}

int allocate_main_memory()
{
  /*
     This routine allocates all the memory needed to contain information
     on particles and groups.
  */

  int igrid, i;
  size_t count_memory;
  struct map
  {
    int lx, ly, lz;
    double oh, am, m1, m2, m3, m4, m5, m6, m7, m8;
  } mymap;
  MPI_Status status;
#ifdef VERBOSE
  size_t *bcast;
  int bsize, n;
  void *last;
#endif

  /* map of memory usage for all tasks */
  mymap.lx = MyGrids[0].GSlocal[_x_];
  mymap.ly = MyGrids[0].GSlocal[_y_];
  mymap.lz = MyGrids[0].GSlocal[_z_];
  mymap.am = (double)memory.all / MBYTE;
  mymap.oh = (double)subbox.Nalloc / (double)MyGrids[0].ParticlesPerTask;
  mymap.m1 = (double)memory.prods / (double)MyGrids[0].ParticlesPerTask;
  mymap.m2 = (double)(memory.fields + memory.fields_to_keep) / (double)MyGrids[0].ParticlesPerTask;
  mymap.m3 = (double)memory.fft / (double)MyGrids[0].ParticlesPerTask;
  mymap.m4 = (double)memory.fmax_total / (double)MyGrids[0].ParticlesPerTask;
  mymap.m5 = (double)memory.frag_prods / (double)MyGrids[0].ParticlesPerTask;
  mymap.m6 = (double)(memory.groups + memory.frag_arrays) / (double)MyGrids[0].ParticlesPerTask;
  mymap.m7 = (double)memory.frag_total / (double)MyGrids[0].ParticlesPerTask;
  mymap.m8 = (double)memory.all / (double)MyGrids[0].ParticlesPerTask;

  if (!ThisTask)
  {
    printf("\n");
    printf("Map of memory usage for all MPI tasks\n");
    printf("Task N.   FFT domain      mem(MB) overhead   products   fields     ffts     fmax  frag pr.  groups fragment  total bytes per particle\n");
  }
  for (i = 0; i < NTasks; i++)
  {
    if (!ThisTask)
    {
      if (i)
        MPI_Recv((void *)&mymap, sizeof(struct map), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
      printf("%6d  %4d-%4d-%4d  %8.0f  %6.1f       %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f\n",
             i, mymap.lx, mymap.ly, mymap.lz, mymap.am, mymap.oh, mymap.m1, mymap.m2, mymap.m3, mymap.m4, mymap.m5, mymap.m6, mymap.m7, mymap.m8);
    }
    else if (ThisTask == i)
      MPI_Send((void *)&mymap, sizeof(struct map), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  }
  if (!ThisTask)
  {
    printf("\n");
    fflush(stdout);
  }

  /* this is just to have a little margin */
  memory.all += 10;

  /* it tests than there is enough space to allocate all needed memory */
  main_memory = (char *)calloc(memory.all, sizeof(char));
  if (main_memory == 0x0)
  {
    printf("ERROR on task %d: I cannot allocate memory\n", ThisTask);
    fflush(stdout);
    return 1;
  }
  free(main_memory);

  /* allocates main memory */
  main_memory = (char *)calloc(memory.first_allocated, sizeof(char));
  if (main_memory == 0x0)
  {
    printf("ERROR on taks %d: I cannot allocate memory after first successful attempt\n", ThisTask);
    fflush(stdout);
    return 1;
  }

  /* sets pointers to the allocated memory */
  products = (product_data *)main_memory;
  count_memory = MyGrids[0].total_local_size * sizeof(product_data);
  ALIGN_MEMORY_BLOCK(count_memory);
  for (igrid = 0; igrid < Ngrids; igrid++)
  {
    kdensity[igrid] = (double *)(main_memory + count_memory);
    count_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
    ALIGN_MEMORY_BLOCK(count_memory);
  }
#ifdef TWO_LPT
  kvector_2LPT = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK(count_memory);
  source_2LPT = kvector_2LPT;
#ifdef THREE_LPT
  kvector_3LPT_1 = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK(count_memory);
  kvector_3LPT_2 = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK(count_memory);
  source_3LPT_1 = kvector_3LPT_1;
  source_3LPT_2 = kvector_3LPT_2;
#endif
#endif
  for (igrid = 0; igrid < Ngrids; igrid++)
  {
    for (i = 0; i < 6; i++)
    {
      second_derivatives[igrid][i] = (double *)(main_memory + count_memory);
      count_memory += MyGrids[igrid].total_local_size * sizeof(double);
      ALIGN_MEMORY_BLOCK(count_memory);
    }
    for (i = 0; i < 3; i++)
      first_derivatives[igrid][i] = second_derivatives[igrid][i];
    density[igrid] = second_derivatives[igrid][0];
  }
  //  for (igrid=0; igrid<Ngrids; igrid++)
  //    {
  //      seedtable[igrid] = (unsigned int*)(main_memory + count_memory);
  //      count_memory += MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);
  //    }
  //
  /* allocates fft vectors */
  for (igrid = 0; igrid < Ngrids; igrid++)
  {
    rvector_fft[igrid] = pfft_alloc_real(MyGrids[igrid].total_local_size_fft);
    cvector_fft[igrid] = pfft_alloc_complex(MyGrids[igrid].total_local_size_fft / 2);
    // printf("Task %d has got alignment for grid %d [ %llu %llu %llu  -  %llu %llu %llu ]\n",
    // ThisTask, igrid, (unsigned long long int)rvector_fft[igrid] % 256, (unsigned long long int)rvector_fft[igrid] % 128, (unsigned long long int)rvector_fft[igrid] % 64,
    //(unsigned long long int)cvector_fft[igrid] % 256, (unsigned long long int)cvector_fft[igrid] % 128, (unsigned long long int)cvector_fft[igrid] % 64);

    if (rvector_fft[igrid] == 0x0 || cvector_fft[igrid] == 0x0)
    {
      printf("ERROR on taks %d: I cannot allocate memory for fft vectors\n", ThisTask);
      fflush(stdout);
      return 1;
    }
  }

  if (!ThisTask)
    printf("[%s] Memory has been successfully allocated\n", fdate());

#ifdef VERBOSE

```


----- FILE: src/build_groups.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
   This file contains code to construct groups (dark matter halos)
   from the list of collapsed particles.
   The main function is build_groups, that performs the group
   construction down to redshift zstop.
   quick_build_groups implements a quick version of build_groups, used
   to determine which part of the boundary region should be
   communicated to perform group construction.
*/

#include "pinocchio.h"
#include <gsl/gsl_sf.h>

#define NCOUNTERS 15

unsigned long long int particle_name;
int good_particle;
pos_data obj1,obj2;

void set_weight(pos_data *);

#ifdef PLC
const gsl_root_fsolver_type *brent;
gsl_root_fsolver *solver;
gsl_function cPLC;
int thisgroup;
int replicate[3];
double brent_err;
#define SAFEPLC 3

double condition_PLC(PRODFLOAT);
double condition_F(double, void *);
int store_PLC(PRODFLOAT);
double find_brent(double, double);
#endif

int build_groups(int Npeaks, double zstop, int first_call)
{

  /* The algorithm for group construction runs as follows:
     loop on all collapsed particle
     + for each particle check its neighbours
     + if it is a peak of Fmax, make it a one-particle halo
     + if it touches only one halo:
       -> check if it should be accreted
       -> if it should, accrete it
       -> otherwise, tag it as filament
     + if it touches more than one halo:
       -> check if it should be accreted to one halo
       -> if it should, choose the one it gets nearest to
       -> then check if all the halo pairs should merge
       -> if they should, merge them
       -> if the particle was not accreted before, re-check it
       -> otherwise, tag it as filament
     + if it touches only filaments, tag it as filament
   */


  int merge[NV][NV], neigh[NV], fil_list[NV][4];
  static int iout,nstep,nstep_p;
  int nn,ifil,this_z,neigrp,nf,pos;
  int iz,i1,j1,k1,skip;
  int ig3,small,large,to_group,accgrp,ig1,ig2;
  int accrflag,nmerge,peak_cond;
  double ratio,best_ratio,d2,r2,cputmp;
  int merge_flag;
  int ibox,jbox,kbox;

  static int last_z=0;

  /*
    counters of various events:
    0 number of peaks
    i number of particles with i neighbours (i<=6)
    7 number of accretion events
    8 number of accretion events before checking a merger
    9 number of accretion events after checking a merger
    10 number of merging events
    11 number of major mergers
    12 number of filament particles
    13 number of accreted filament particles
    14 number of good halos
  */

  static unsigned long long counters[NCOUNTERS],all_counters[NCOUNTERS];

#ifdef PLC
  static int plc_started=0, last_check_done=0;
  int irep, save, mysave, storex, storey, storez;
  static double NextF_PLC, DeltaF_PLC;
  double aa, bb, Fplc;
#endif

  if (first_call)
    {
      /* this part contains initializations that should be done at
         first call */
      ngroups=FILAMENT;           /* group list starts from FILAMENTS */
      /* filaments are not grouped */
      for (i1=0; i1<=FILAMENT; i1++)
        {
          groups[i1].point=-1;
          groups[i1].bottom=-1;
          groups[i1].good=0;
        }
      iout=0;

      /* sets the counter to zero */
      for (i1=0; i1<NCOUNTERS; i1++)
        {
          counters[i1]=0;
          all_counters[i1]=0;
        }


      if (!ThisTask)
        printf("[%s] Starting the fragmentation process to redshift %7.4f\n",fdate(),zstop);

      /* nstep is the number of collapsed particles that will be
         checked by the code. In classic fragmentation this is the
         number of particles that collapse by the end of the run. In
         default fragmentation this selection is performed at
         distribution time, so this is the number of stored particles.
      */
#ifdef CLASSIC_FRAGMENTATION
      nstep=0;
      while (frag[indices[nstep]].Fmax >= outputs.Flast)
        nstep++;
#else
      nstep=subbox.Nstored;
#endif
      nstep_p=nstep/20;

#ifdef PLC
      /* initialization of the root finding routine used by the PLC code */
      cPLC.function = &condition_F;
      brent  = gsl_root_fsolver_brent;
      solver = gsl_root_fsolver_alloc (brent);
      DeltaF_PLC  = 0.9; // CAPIRE COME FISSARLO
      NextF_PLC   = plc.Fstart * DeltaF_PLC;
      brent_err   = 1.e-2 * params.InterPartDist;
      plc.Nstored = plc.Nstored_last = 0;
#endif

      first_call=0;
    }
  else
    {
      if (!ThisTask)
        printf("[%s] Restarting the fragmentation process to redshift %7.4f\n",fdate(),zstop);
    }

  /************************************************************************
                    START OF THE CYCLE ON COLLAPSED PARTICLES
   ************************************************************************/
  for (this_z=last_z; this_z<nstep; this_z++)
    {
      /* In classic fragmentation the particles must be addressed in
         order of collapse time, while in default fragmentation they
         are already in that order */
#ifdef CLASSIC_FRAGMENTATION
      iz=indices[this_z];
#else
      iz=this_z;
#endif


#ifdef PLC
      /* here it checks if it is time to write the stored PLC halos.
         If it is the case, it waits for all the tasks to get to the
         same point and writes the catalog up to this point
       */
      if (plc_started && frag[iz].Fmax < NextF_PLC && frag[iz].Fmax >= plc.Fstop)
        {
          if (!ThisTask)
            printf("[%s] Syncing tasks for PLC...  Nmax=%d, Nstored=%d, Nstored_last=%d, Fmax=%f\n",fdate(),plc.Nmax, plc.Nstored, plc.Nstored_last,frag[iz].Fmax);

          /* each task checks if the plc buffer is full and it is
             necessary to write the plc catalog. The criterion is that
             you need at least space for SAFEPLC times the number of
             halos that have been updated after last check */
          mysave = (plc.Nmax - plc.Nstored < SAFEPLC * (plc.Nstored - plc.Nstored_last));

          MPI_Reduce(&mysave, &save, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Bcast(&save, 1, MPI_INT, 0, MPI_COMM_WORLD);

          if (!ThisTask)
            printf("[%s] %d tasks require to store PLC halos at F=%6.3F...\n",fdate(),save,frag[iz].Fmax);

          if (save)
            {
              /* Save PLC halos if at least one task asks to */
              if (write_PLC(0))
                return 1;
              plc.Nstored=0;
            }

          /* next update will be at this time */
          NextF_PLC *= DeltaF_PLC;
          plc.Nstored_last = plc.Nstored;

        }
#endif

      /* More initializations */

      neigrp=0;               /* number of neighbouring groups */
      nf=0;                   /* number of neighbouring filament points */
      accrflag=0;             /* if =1 all the neighbouring filaments are accreted */
      for (i1=0; i1<NV; i1++)
        neigh[i1]=0;          /* number of neighbours */

      /*******************************************/
      /* PARTICLE COORDINATES AND PEAK CONDITION */
      /*******************************************/
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);
#endif

      /* skips the peak condition if the point is at the border (and PBCs are not active) */
      skip=0;
      if ( !subbox.pbc[_x_] && (ibox==0 || ibox==subbox.Lgwbl[_x_]-1) ) ++skip;
      if ( !subbox.pbc[_y_] && (jbox==0 || jbox==subbox.Lgwbl[_y_]-1) ) ++skip;
      if ( !subbox.pbc[_z_] && (kbox==0 || kbox==subbox.Lgwbl[_z_]-1) ) ++skip;

      particle_name =
        COORD_TO_INDEX((long long)((ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]),
                       (long long)((jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]),
                       (long long)((kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]),
                       MyGrids[0].GSglobal);

      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] &&
                        jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] &&
                        kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      if (!skip)
        {
          peak_cond=1;
          /* checks whether the neighbouring particles collapse later */

          for (nn=0; nn<NV; nn++)
            {
              /* coordinates of the neighbouring particle */
              switch (nn)
                {
                case 0:
                  i1=( subbox.pbc[_x_] && ibox==0 ? subbox.Lgwbl[_x_]-1 : ibox-1 );
                  j1=jbox;
                  k1=kbox;
                  break;
                case 1:
                  i1=( subbox.pbc[_x_] && ibox==subbox.Lgwbl[_x_]-1 ? 0 : ibox+1 );
                  j1=jbox;
                  k1=kbox;
                  break;
                case 2:
                  i1=ibox;
                  j1=( subbox.pbc[_y_] && jbox==0 ? subbox.Lgwbl[_y_]-1 : jbox-1 );
                  k1=kbox;
                  break;
                case 3:
                  i1=ibox;
                  j1=( subbox.pbc[_y_] && jbox==subbox.Lgwbl[_y_]-1 ? 0 : jbox+1 );
                  k1=kbox;
                  break;
                case 4:
                  i1=ibox;
                  j1=jbox;
                  k1=( subbox.pbc[_z_] && kbox==0 ? subbox.Lgwbl[_z_]-1 : kbox-1 );
                  break;
                case 5:
                  i1=ibox;
                  j1=jbox;
                  k1=( subbox.pbc[_z_] && kbox==subbox.Lgwbl[_z_]-1 ? 0 : kbox+1 );
                  break;
                }

              /* accessing the neighbouring particle differs in
                 classic and default fragmentation: in classic
                 fragmentation the information on the neighbour is
                 immediately obtained, in standard fragmentation the
                 neighbour must be seeked in the particle list */
#ifdef CLASSIC_FRAGMENTATION
              pos = COORD_TO_INDEX(i1,j1,k1,subbox.Lgwbl);
              neigh[nn] = group_ID[pos];
              peak_cond &= (frag[iz].Fmax > frag[pos].Fmax);
#else
              pos = find_location(i1,j1,k1);
              if (pos>=0)
                {
                  neigh[nn] = group_ID[indices[pos]];
                  peak_cond &= (frag[iz].Fmax > frag[indices[pos]].Fmax);
                }
              else
                neigh[nn] = 0;
#endif

              /* neighbouring filaments are stored separately */
              if (neigh[nn]==FILAMENT)
                {
                  neigh[nn]=0;
                  fil_list[nf][0]=i1;
                  fil_list[nf][1]=j1;
                  fil_list[nf][2]=k1;
#ifdef CLASSIC_FRAGMENTATION
                  fil_list[nf][3]=pos;
#else
                  fil_list[nf][3]=indices[pos];
#endif
                  nf++;
                }

            }

          /* Cleans the list of neighbouring groups removing duplicates */
          clean_list(neigh);

          /* Number of neighbouring groups */
          for (nn=neigrp=0; nn<NV; nn++)
            if (neigh[nn]>FILAMENT)
              neigrp++;

          if (neigrp>0 && good_particle)
            counters[neigrp]++;

#ifdef PLC
          /* Past light cone on-the-fly reconstruction: */

          /* is it time to reconstruct the PLC? */
          if (frag[iz].Fmax<plc.Fstart && frag[iz].Fmax>=plc.Fstop)
            {

              /* check if this is the first call */
              if (!plc_started)
                {
                  plc_started=1;
                  if (!ThisTask)
                    printf("[%s] Starting PLC reconstruction, Task 0 will store at most %d halos\n",fdate(),plc.Nmax);
                  cputmp=MPI_Wtime();
                }

              /* PBCs are switched off for this check */
              storex=subbox.pbc[_x_];
              storey=subbox.pbc[_y_];
              storez=subbox.pbc[_z_];
              subbox.pbc[_x_]=subbox.pbc[_y_]=subbox.pbc[_z_]=0;

              /* the check is performed on all neighbouring groups */
              for (ig1=0; ig1<neigrp; ig1++)
                {
                  /* is the group good and massive enough? */
                  if (neigh[ig1] > FILAMENT &&
                      groups[neigh[ig1]].good &&
                      groups[neigh[ig1]].Mass >= params.MinHaloMass)
                    {
                      thisgroup=neigh[ig1];

                      /* loop on replications */
                      /* this loop may be threaded or ported to GPUs,
                         though it's relatively fast */
                      for (irep=0; irep<plc.Nreplications; irep++)
                        /* checks that the redshift falls in the
                           range of the replication */
                        if (!(frag[iz].Fmax > plc.repls[irep].F1 ||
                              groups[thisgroup].Flast < plc.repls[irep].F2))
                          {
                            replicate[0]=plc.repls[irep].i;
                            replicate[1]=plc.repls[irep].j;
                            replicate[2]=plc.repls[irep].k;
                            /* this computes the difference between
```


----- FILE: src/collapse_times.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "pinocchio.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <immintrin.h>

/*------------------------------------------------------- Macros declaration --------------------------------------------------------*/


#define SMALL ((double)1.e-20)
//#define TRILINEAR
#define BILINEAR_SPLINE
//#define ALL_SPLINE
//#define ASCII
//#define HISTO


/*------------------------------------------------------- Functions definition -------------------------------------------------------*/


/* Ellipsoidal solution at 3rd perturbative order in two different ways */

__attribute__((always_inline)) double ell_classic               ( int, double, double, double);
__attribute__((always_inline)) double ell_sng                   ( int, double, double, double);

/* Collapase time calculation */

__attribute__((always_inline)) double ell                       ( int, double, double, double);

/* Calculation of inverse collpase time */

__attribute__((always_inline)) double inverse_collapse_time     ( int, double *, double *, double *, double *, int *);

/* Interpolation of collpase time from a table containing collpase times */

#ifdef TABULATED_CT

__attribute__((always_inline)) double interpolate_collapse_time (int, double, double, double);

#endif

#ifdef MOD_GRAV_FR

__attribute__((always_inline)) double ForceModification (double, double ,double);

#endif
/* Order inverse collpase time */

__attribute__((always_inline)) void ord                         (double *,double *,double *);


/*------------------------------------------------------- Global Variables declaration ----------------------------------------------------*/


/* Declaring the replicated variables that will be use by each single thread */
/* NOTE: Threadprivate directive specifies that variables are replicated, with each thread having its own copy. It's a declarative directive */

#if defined( _OPENMP )

/* Cpu_time for elliptical collapse */

double cputime_ell;
#pragma omp threadprivate(cputime_ell) // threadprivate directive specifies that variables are replicated, with each thread having its own copy


/* It will be use as a fail "flag" indicating that for some reason the calculation of collapse time failed */

int fails;
#pragma omp threadprivate(fails)

/* If there's no _OPENMP macro then define cputime_ell */

#else
#define cputime_ell            cputime.ell
#endif

/* Declaring smoothing radius */

int ismooth;


/*------------------------------------------------------- Functions implementation --------------------------------------------------------*/

/* Classical ellipsoidal collapse solution */

inline double ell_classic(int ismooth, double l1, double l2, double l3) {

    /*
    This routine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
    */

    /* Local variables declaration */
    double ell;
    double del = l1 + l2 + l3;
    double det = l1 * l2 * l3;

    /* Vanishing lambda1 eigenvalue case */
    if (fabs(l1) < SMALL) {
        ell = -0.1;
    }
    /* Not vanishing lambda1 eigenvalue case */
    else {
        double den = det / 126. + 5. * l1 * del * (del - l1) / 84.;
        /* Check 1st perturbative order conditions */
        if (fabs(den) < SMALL) {
            if (fabs(del - l1) < SMALL) {
                if (l1 > 0.0) {
                    /* Zel'dovich approximation */
                    ell = 1. / l1;
                } else {
                    ell = -.1;
                }
            } else {
                /* Check 2nd perturbative order conditions */
                double dis = 7. * l1 * (l1 + 6. * del);
                if (dis < 0.0) {
                    ell = -.1;
                } else {
                    /* 2nd order solution */
                    ell = (7. * l1 - sqrt(dis)) / (3. * l1 * (l1 - del));
                    if (ell < 0.) {
                        ell = -.1;
                    }
                }
            }
        } else {
            /* 3rd order perturbative solution. For more details about the equations implemented, see Monaco 1996a */

            /* Intermediate values */
            double rden = 1.0 / den;
            double a1 = 3. * l1 * (del - l1) / 14. * rden;
            double a1_2 = a1 * a1;
            double a2 = l1 * rden;
            double a3 = -1.0 * rden;

            /* The collapse time b_c will be a combination of R, Q, and D == den */
            double q = (a1_2 - 3. * a2) / 9.;
            double r = (2. * a1_2 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
            double r_2_q_3 = r * r - q * q * q;

            /* Check 3rd perturbative order conditions */

            /* ---------------- Case 1 --------------- */
            /* If R^2 - Q^2 > 0, which is valid for spherical and quasi-spherical perturbations */

            /* 3rd order solution */
            if (r_2_q_3 > 0) {
                double fabs_r = fabs(r);
                double sq = pow(sqrt(r_2_q_3) + fabs_r, 0.333333333333333);
                ell = -fabs_r / r * (sq + q / sq) - a1 / 3.;
                if (ell < 0.) {
                    ell = -.1;
                }
            }

            /* ---------------- Case 2 --------------- */
            /* The solution has to be chosen as the smallest non-negative one between s1, s2, and s3 */

            /* 3rd order solution */
            else {
                double sq = 2 * sqrt(q);
                double inv_3 = 1.0 / 3;
                double t = acos(2 * r / q / sq);
                double s1 = -sq * cos(t * inv_3) - a1 * inv_3;
                double s2 = -sq * cos((t + 2. * PI) * inv_3) - a1 * inv_3;
                double s3 = -sq * cos((t + 4. * PI) * inv_3) - a1 * inv_3;
                if (s1 < 0.) {
                    s1 = 1.e10;
                }
                if (s2 < 0.) {
                    s2 = 1.e10;
                }
                if (s3 < 0.) {
                    s3 = 1.e10;
                }
                ell = (s1 < s2 ? s1 : s2);
                ell = (s3 < ell ? s3 : ell);
                if (ell == 1.e10) {
                    ell = -.1;
                }
            }
        }
    }

    if (del > 0. && ell > 0.) {
        double inv_del = 1.0 / del;
        ell += -.364 * inv_del * exp(-6.5 * (l1 - l2) * inv_del - 2.8 * (l2 - l3) * inv_del);
    }

    return ell;
}

/* Ellipsoidal collapse following Nadkarni-Ghosh & Singhal (2016) */

/* We follow the dynamics of triaxial collapse in terms of eigenvalues of the deformation tensor (lambda_a), the velocity derivative tensor (lambda_v) and the gravity Hessian
(lambda_d). The idea is: starting from BM96, where the dynamic is characterized by the evolution of the three principal axes (the physical coordinates is r=a_i(t)q(i = 1,2,3)).
In this case the evolution of the ellipse is completely determined once six parameters are known: the three axes lengths and their velocities at some initial epoch a_init.
An alternate description of the ellipse can be given by a set of nine dimensionless parameters lambda_a, lambda_v, lambda_d. In this description when an axis is collpasing
lamda_a ----> 1, whereas for an expanding axes lambda_a -----> -inf.
delta = lambda_d_1 + lambda_d_2 + lambda_d_3
The nine (dimensionless) eigenvalues completely characterize the density, velocity and shape perturbations in this model of ellipsoidal collapse.
For the net system for the nine eigenvalues see Nadkarni-Ghosh & Singhal (2016) pag. 5 */


/* We have then to solve a system of differential equations          */
/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */

int sng_system(double t, const double y[], double f[], void *sng_par) {

    /* Local variables declaration */
	int i,j;
	double sum;

	/* Needed cosmological parameters calculation at the given z */
	double omegam = OmegaMatter(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */
	double omegal = OmegaLambda(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */

	/* y array contains the eigenvalues lambda_a, lambda_v and lambda_d */
	/* In this case y[6] = lambda_d_1, y[7] = lambda_d_2, y[8] = lambda_d_3 */
	double delta = y[6] + y[7] + y[8];

    /* Calculating r.h.s sums of lambda_d_i evolution */
	/* Loop from 0 -----> 3 beacause lambda_i has 3 component (it's the same for lambda_a and lambda_v) */

	/* ---------------- Sum  i != j --------------- */
	 for (i = 0; i < 3; i++) {
        sum = 0.;
        for (j = 0; j < 3; j++) {
            if (i == j || y[i] == y[j]) {
                continue;
            } else {
                sum += (y[j + 6] - y[i + 6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i + 3]) -
                       (1. - y[j]) * (1. - y[j]) * (1. + y[j + 3])) /
                       ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
            }
        }

		/* ---------------- Gathering r.h.s of lambda_a_i equation --------------- */
        f[i] = (y[i + 3] * (y[i] - 1.0)) / t;

		  /* ---------------- Gathering r.h.s of lambda_v_i equation --------------- */
        /* NOTE: here is the part where models of modified gravity can be introduced */
        f[i + 3] = (0.5 * (y[i + 3] * (omegam - 2.0 * omegal - 2.0)

#ifdef MOD_GRAV_FR
                       - 3.0 * omegam * y[i + 6] * (1. + ForceModification(*(double *)sng_par, t, delta))
#else
                       - 3.0 * omegam * y[i + 6]
#endif
                       - 2.0 * y[i + 3] * y[i + 3])) / t;

		/* ---------------- Gathering r.h.s of lambda_d_i equation --------------- */
        f[i + 6] = ((5. / 6. + y[i + 6]) *
                    ((3. + y[3] + y[4] + y[5]) - (1. + delta) / (2.5 + delta) * (y[3] + y[4] + y[5])) -
                    (2.5 + delta) * (1. + y[i + 3]) + sum) / t;
    }

    return GSL_SUCCESS;
}

/* --------------------------------------------------- Modified gravity model --------------------------------------------------------*/

#ifdef MOD_GRAV_FR

inline double ForceModification(double size, double a, double delta) {

    double ff        = 4. * params.OmegaLambda / params.Omega0;
    double thickness = FR0 / params.Omega0 / pow(H_over_c * size, 2.0) *
                       pow(a, 7.) * pow((1. + delta), -1. / 3.) *
                       (pow((1.0 + ff) / (1.0 + ff * pow(a, 3.)), 2.0) -
                       pow((1.0 + ff) / (1.0 + delta + ff * pow(a, 3.)), 2.0));

    double F3 = (thickness * (3. + thickness * (-3. + thickness)));
    if (F3 < 0.) {
        F3 = 0.;
    }
    return (F3 < 1. ? F3 / 3. : 1. / 3);
}

#endif

/* Solving the ODE system for the ellipsoidal collapse following SNG */

inline double  ell_sng(int ismooth, double l1, double l2, double l3) {

	/* Needed variables for the integrations step */
	double ode_param, hh = 1.e-6;  // hh = initial step-size, ode_param = arbitrary parameters of the system == Smoothing_radius in our case

    /* Setting integration time */
	double amin = 1.e-5, amax = 5.0;
	double mya = amin;

	/* Step type: Explicit embedded Runge-Kutta-Fehlberg (4, 5) method */
	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

	/* Newly allocated instance of a stepping function of type T for a system of dim = 9 */
	gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9);

	/* The control function examines the proposed change to the solution produced by a stepping function and attempts to determine the optimal step-size for a user-specified level of error.
	The standard control object is a four parameter heuristic based on absolute and relative errors eps_abs and eps_rel,and scaling factors a_y and a_dydt
	for the system state y(t) and derivatives y'(t) respectively. */
	gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);

	/* The evolution function combines the results of a stepping function and control function to reliably advance the solution forward one step using an acceptable step-size.
	This function returns a pointer to a newly allocated instance of an evolution function for a system of dim = 9 */
	gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);

	/* A system of equations is defined using the gsl_odeiv2_system datatype
	This data type defines a general ODE system with arbitrary parameters. In this case we need the r.h.s of the system that will be solved: sng_system
	The vector of derivatives elements: jac
	The dimension of the system of equations: 9
	A pointer to the arbitrary parameters of the system: (void*)&ode_param */
	gsl_odeiv2_system   ode_sys   = {sng_system, jac, 9, (void*)&ode_param};

/* GrowingMode is linear growing mode, interpolated on the grid. See cosmo.c for details (double GrowingMode(double z, double k))*/
#ifdef SCALE_DEPENDENT

	double D_in = GrowingMode(1./amin-1.,params.k_for_GM/Smoothing.Radius[ismooth]*params.InterPartDist);

#else

	double D_in =  GrowingMode(1./amin-1.,1./Smoothing.Radius[ismooth]);

#endif

	double y[9] = {l1*D_in, l2*D_in, l3*D_in,
		       l1*D_in/(l1*D_in - 1.), l2*D_in/(l2*D_in - 1.), l3*D_in/(l3*D_in - 1.),
		       l1*D_in, l2*D_in, l3*D_in};

/* Assignment of ode_param in two different cases */
#ifdef MOD_GRAV_FR

	if (ismooth<Smoothing.Nsmooth-1)
	  {
	    ode_param = Smoothing.Radius[ismooth];
	  }
	else
	  {
	    ode_param = Smoothing.Radius[ismooth-1];
	  }
#endif

    /*---------------- Integration step --------------- */
	double olda = mya;
	double oldlam = y[0];

	while (mya < amax) {
        /* gsl_odeiv2_evolve_apply advances the system from time t and position y using the stepping function selected before i.e. rkf45.
		The new time and position are stored in t and y on output*/
        int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &mya, amax, &hh, y);

        /* Correct integration check */
        if (status != GSL_SUCCESS) {
            printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
            fflush(stdout);
            return -1;
        }

        /* ---------------- Update ellipsoid axis --------------- */
        if (y[0] >= 0.99999) {

            return olda + (1. - oldlam) * (mya - olda) / (y[0] - oldlam);

        }
    }

    /* In this case the ellipsoid does not collapse */
    return 0;
}
```


----- FILE: src/cosmo.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"
#include "def_splines.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define NVAR (1 + NkBINS * 8)
#define NBB 10
#ifdef NORADIATION
#define OMEGARAD_H2 ((double)0.0)
#else
#define OMEGARAD_H2 ((double)4.2e-5)
#endif
#define UnitLength_in_cm ((double)3.085678e24)
#define HUBBLETIME_GYR ((double)3.085678e24 / (double)1.e7 / (double)3.1558150e16)
#define DELTA_C ((double)1.686)
#define SHAPE_EFST ((double)0.21)

// #define FOMEGA_GAMMA 0.554
// TUTTE LE CONDIZIONI SULLE DIRETTIVE DEVONO ESSERE MESSE INSIEME
#if defined(FOMEGA_GAMMA) && defined(SCALE_DEPENDENT)
#error Do not use FOMEGA_GAMMA with SCALE_DEPENDENT
#endif

static int Today;
static int WhichSpectrum, NPowerTable = 0, NtabEoS = 0;
static double PkNorm, MatterDensity, OmegaK, OmegaRad;

/* declaration of gsl quantities */

#ifdef SCALE_DEPENDENT
static double kmin, kmax;
#endif

int system_of_ODEs(double, const double[], double *, void *);
int system_of_ODEs_small(double, const double[], double *, void *);
int read_TabulatedEoS(void);
#ifdef READ_HUBBLE_TABLE
int read_TabulatedHubble(void);
#endif
int initialize_PowerSpectrum(void);
int normalize_PowerSpectrum(void);
int read_Pk_from_file(void);
double IntegrandForEoS(double, void *);
double DE_EquationOfState(double);
double IntegrandComovingDistance(double, void *);
double ComputeMassVariance(double);
double ComputeDisplVariance(double);
#ifdef READ_PK_TABLE
int read_Pk_table_from_CAMB(double *, double *, double *, double *, double *, double *, double *, double *, double *);
#endif

/**************************/
/* INITIALIZATION SECTION */
/**************************/

int initialize_cosmology()
{

  /*
    Computes the following functions:
    Scale factor, growth first, second and third-order LPT, cosmic time,
    comoving and diameter distance on a grid of values to be interpolated
  */

  double ode_param;
  double y[NVAR], x1, x2, hh, norm, result, error, SqrtOmegaK, R0, k;
  int status = GSL_SUCCESS, i, j;
#ifdef SCALE_DEPENDENT
  int ik;
#endif
  char filename[LBLENGTH];
  FILE *fd;
  double log_amin = -4., dloga = -log_amin / (double)(NBINS - NBB);

  double *scalef, *cosmtime, *grow1, *grow2, *IntEoS, *comvdist, *diamdist,
      *fomega1, *fomega2, *grow31, *grow32, *fomega31, *fomega32;

  gsl_function Function_cosmo;

#ifdef MOD_GRAV_FR
  H_over_c = 100. / SPEEDOFLIGHT;
#endif
  OmegaRad = OMEGARAD_H2 / params.Hubble100 / params.Hubble100;
  OmegaK = 1.0 - params.Omega0 - params.OmegaLambda - OmegaRad;
  SqrtOmegaK = sqrt(fabs(OmegaK));
  MatterDensity = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0;
  if (params.DEw0 == -1 && params.DEwa == 0 && !strcmp(params.TabulatedEoSfile, "no"))
    params.simpleLambda = 1;
  else
    params.simpleLambda = 0;

  /* allocation of (most) splines */
  SPLINE = (gsl_spline **)calloc(NSPLINES, sizeof(gsl_spline *));
  for (i = 0; i < NSPLINES - 3; i++)
    if (i != SP_COMVDIST && i != SP_DIAMDIST)
      SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS);
    else
      SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS - NBB);
  ACCEL = (gsl_interp_accel **)calloc(NSPLINES, sizeof(gsl_interp_accel *));
  for (i = 0; i < NSPLINES; i++)
    ACCEL[i] = gsl_interp_accel_alloc();

  /* if needed, read the tabulated Equation of State of the dark energy
     and initialize its spline, then compute the integrand for the DE EoS */
  if (!params.simpleLambda)
  {
    if (strcmp(params.TabulatedEoSfile, "no"))
    {
      if (read_TabulatedEoS()) /* this allocates and initializes the spline */
        return 1;
    }
    else
      NtabEoS = 0;

    SPLINE[SP_INTEOS] = gsl_spline_alloc(gsl_interp_cspline, NBINS);

    scalef = (double *)malloc(NBINS * sizeof(double));
    IntEoS = (double *)malloc(NBINS * sizeof(double));
    Function_cosmo.function = &IntegrandForEoS;
    for (i = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + (i + 1) * dloga);
      gsl_integration_qags(&Function_cosmo, x2, 1.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      scalef[i] = log_amin + (i + 1) * dloga;
      IntEoS[i] = result;
    }

    gsl_spline_init(SPLINE[SP_INTEOS], scalef, IntEoS, NBINS);

    free(IntEoS);
    free(scalef);
  }

#ifdef READ_HUBBLE_TABLE
  /* If requested, read tabulated H(z) and create its spline over log10(a). */
  if (read_TabulatedHubble())
    return 1;
#endif

#ifdef SCALE_DEPENDENT
  kmin = pow(10., LOGKMIN);
  kmax = pow(10., LOGKMIN + (NkBINS - 1) * DELTALOGK);
#endif

  /* The power spectrum is initialized; in case, it is read from file(s) */
  if (initialize_PowerSpectrum())
    return 1;

  /* allocation of vectors for interpolation */
  scalef = (double *)malloc(NBINS * sizeof(double));
  cosmtime = (double *)malloc(NBINS * sizeof(double));
  comvdist = (double *)malloc((NBINS - NBB) * sizeof(double));
  diamdist = (double *)malloc((NBINS - NBB) * sizeof(double));
  grow1 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow2 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow31 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  grow32 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega1 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega2 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega31 = (double *)malloc(NBINS * NkBINS * sizeof(double));
  fomega32 = (double *)malloc(NBINS * NkBINS * sizeof(double));

#ifdef READ_PK_TABLE
  if (WhichSpectrum != 5)
#endif
  {
    /* Runge-Kutta integration of cosmic time and growth rate */
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *ode_s = gsl_odeiv2_step_alloc(T, NVAR);
    gsl_odeiv2_control *ode_c = gsl_odeiv2_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
    gsl_odeiv2_evolve *ode_e = gsl_odeiv2_evolve_alloc(NVAR);
    gsl_odeiv2_system ode_sys = {system_of_ODEs, jac, NVAR, (void *)&ode_param};

    /* ICs for the runge-kutta integration */
    x1 = pow(10., log_amin - 2.);  /*  initial value of scale factor a */
    y[0] = 2. / 3. * pow(x1, 1.5); /*  initial value of t(a)*Hubble0   */

    /* this is valid for scale-independent and scale-dependent functions */
    for (j = 0; j < NkBINS; j++)
    {
      y[1 + j * 8] = 1.0;                      /*  initial value of dD1/da   */
      y[2 + j * 8] = x1;                       /*  initial value of D1(a)    */
      y[3 + j * 8] = -6. / 7. * x1;            /*  initial value of dD2/da   */
      y[4 + j * 8] = -3. / 7. * x1 * x1;       /*  initial value of D2(a)    */
      y[5 + j * 8] = -x1 * x1;                 /*  initial value of dD^3a/da */
      y[6 + j * 8] = -1. / 3. * x1 * x1 * x1;  /*  initial value of D^3a     */
      y[7 + j * 8] = 10. / 7. * x1 * x1;       /*  initial value of dD^3b/da */
      y[8 + j * 8] = 10. / 21. * x1 * x1 * x1; /*  initial value of D^3b     */
    }

    hh = x1 / 10.; /*  initial guess of time-step */

    /* this function will be integrated within the loop */
    Function_cosmo.function = &IntegrandComovingDistance;

    /***********************************************/
    /* ODE integration of time-dependent functions */
    /***********************************************/
    for (i = 0, Today = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + i * dloga);
      if (fabs(log_amin + i * dloga) < dloga / 10.)
        x2 = 1.0;

      /* integration of ODE system */
      while (x1 < x2 && status == GSL_SUCCESS)
      {
        status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
        if (status != GSL_SUCCESS)
        {
          printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
          fflush(stdout);
          return 1;
        }
      }

      scalef[i] = x2;
      cosmtime[i] = log10(y[0] * HUBBLETIME_GYR / params.Hubble100);
      if (!Today && x2 >= 1.)
        Today = i;

      for (j = 0; j < NkBINS; j++) /* First-order growth rate */
        grow1[i + j * NBINS] = y[2 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Second-order growth rate */
        grow2[i + j * NBINS] = -y[4 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order first growth rate */
        grow31[i + j * NBINS] = -y[6 + j * 8] / 3.;
      for (j = 0; j < NkBINS; j++) /* Third-order second growth rate */
        grow32[i + j * NBINS] = y[8 + j * 8] / 4.;

      for (j = 0; j < NkBINS; j++) /* First-order f(Omega) */
        fomega1[i + j * NBINS] = x2 * y[1 + j * 8] / y[2 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Second-order f(Omega) */
        fomega2[i + j * NBINS] = x2 * y[3 + j * 8] / y[4 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order first f(Omega) */
        fomega31[i + j * NBINS] = x2 * y[5 + j * 8] / y[6 + j * 8];
      for (j = 0; j < NkBINS; j++) /* Third-order second f(Omega) */
        fomega32[i + j * NBINS] = x2 * y[7 + j * 8] / y[8 + j * 8];

      /* calculation of comoving distance (Mpc) for a generic cosmology (i.e. flat, open or closed)
         and generic equation of state for the DE component  */
      if (i < NBINS - NBB)
      {
        gsl_integration_qags(&Function_cosmo, 0.0, 1. / x2 - 1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
        comvdist[i] = SPEEDOFLIGHT * result;
        if (fabs(OmegaK) < 1.e-4)
          diamdist[i] = x2 * comvdist[i];
        else if (OmegaK < 0)
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sin(comvdist[i] / R0);
        }
        else
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sinh(comvdist[i] / R0);
        }
      }

      /* closing the loop on integrations */
      x1 = x2;
    }

    /* normalization of the first- and second- order growth rate */
    /* this is valid for LambdaCDM; for scale-dependent growth due to modified gravity,
       the power spectrum is given as the LambdaCDM P(k) extrapolated at z=0, but this
       is valid only at high redshift; this normalization is still correct
       when the k=0 growth rate (identical to LambdaCDM) is used at all scales */
    norm = grow1[Today];
    for (i = 0; i < NBINS * NkBINS; i++)
    {
      grow1[i] /= norm;
      grow2[i] /= norm * norm;
      grow31[i] /= norm * norm * norm;
      grow32[i] /= norm * norm * norm;
    }
  }
#ifdef READ_PK_TABLE
  else
  {
    if (!ThisTask)
      printf("Only the cosmic time is integrated, the growth rate is read from CAMB files\n");

    /* in this case the integration is limited only to the cosmic time */
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *ode_s = gsl_odeiv2_step_alloc(T, 1);
    gsl_odeiv2_control *ode_c = gsl_odeiv2_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
    gsl_odeiv2_evolve *ode_e = gsl_odeiv2_evolve_alloc(1);
    gsl_odeiv2_system ode_sys = {system_of_ODEs_small, jac, 1, (void *)&ode_param};

    /* ICs for the runge-kutta integration of the cosmic time only */
    x1 = pow(10., log_amin - 2.);  /*  initial value of scale factor a */
    y[0] = 2. / 3. * pow(x1, 1.5); /*  initial value of t(a)*Hubble0   */
    hh = x1 / 10.;                 /*  initial guess of time-step */

    /* this function will be integrated within the loop */
    Function_cosmo.function = &IntegrandComovingDistance;

    /***********************************************/
    /* ODE integration of time-dependent functions */
    /***********************************************/
    for (i = 0, Today = 0; i < NBINS; i++)
    {
      x2 = pow(10., log_amin + i * dloga);
      if (fabs(log_amin + i * dloga) < dloga / 10.)
        x2 = 1.0;

      /* integration of ODE system */
      while (x1 < x2 && status == GSL_SUCCESS)
      {
        status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
        if (status != GSL_SUCCESS)
        {
          printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
          fflush(stdout);
          return 1;
        }
      }

      scalef[i] = x2;
      cosmtime[i] = log10(y[0] * HUBBLETIME_GYR / params.Hubble100);
      if (!Today && x2 >= 1.)
        Today = i;

      /* calculation of comoving distance (Mpc) for a generic cosmology (i.e. flat, open or closed)
         and generic equation of state for the DE component  */
      if (i < NBINS - NBB)
      {
        gsl_integration_qags(&Function_cosmo, 0.0, 1. / x2 - 1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
        comvdist[i] = SPEEDOFLIGHT * result;
        if (fabs(OmegaK) < 1.e-4)
          diamdist[i] = x2 * comvdist[i];
        else if (OmegaK < 0)
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sin(comvdist[i] / R0);
        }
        else
        {
          R0 = SPEEDOFLIGHT / params.Hubble100 / 100. / SqrtOmegaK;
          diamdist[i] = x2 * R0 * sinh(comvdist[i] / R0);
        }
      }

      /* closing the loop on integrations */
      x1 = x2;
    }

    /* the growth rates are set by reading the CAMB power spectra */
    if (read_Pk_table_from_CAMB(scalef, grow1, grow2, grow31, grow32, fomega1, fomega2, fomega31, fomega32))
      return 1;
  }
#endif

  /* these quantities will be interpolated logarithmically */
  for (i = 0; i < NBINS; i++)
    scalef[i] = log10(scalef[i]);
  for (j = 0; j < NkBINS; j++)
    for (i = 0; i < NBINS; i++)
    {
      grow1[i + j * NBINS] = log10(grow1[i + j * NBINS]);
      grow2[i + j * NBINS] = log10(grow2[i + j * NBINS]);
      grow31[i + j * NBINS] = log10(grow31[i + j * NBINS]);
      grow32[i + j * NBINS] = log10(grow32[i + j * NBINS]);
    }

  /* initialization of spline interpolations of time-dependent quantities */
  gsl_spline_init(SPLINE[SP_TIME], scalef, cosmtime, NBINS);
  gsl_spline_init(SPLINE[SP_INVTIME], cosmtime, scalef, NBINS);
  gsl_spline_init(SPLINE[SP_COMVDIST], scalef, comvdist, NBINS - NBB);
  gsl_spline_init(SPLINE[SP_DIAMDIST], scalef, diamdist, NBINS - NBB);
  /* inverse grow is always defined on the first growth rate */
```


----- FILE: src/def_splines.h -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* #if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC) */
/* #error Trying to compile with ELL_CLASSIC and SCALE_DEPENDENT together */
/* #endif */

#if defined(MOD_GRAV_FR) && !defined(SCALE_DEPENDENT)
#define SCALE_DEPENDENT
#endif

#ifndef SCALE_DEPENDENT
#define NkBINS 1
#else
#define NkBINS 10
#define LOGKMIN ((double)-3.0)
#define DELTALOGK ((double)0.5)
#endif

#define SP_TIME 0
#define SP_INVTIME 1
#define SP_COMVDIST 2
#define SP_DIAMDIST 3
#define SP_INVGROW 4
#define SP_MASSVAR 5
#define SP_DISPVAR 6
#define SP_RADIUS 7
#define SP_DVARDR 8

#define SP_GROW1 (9)
#define SP_GROW2 (9 + NkBINS)
#define SP_GROW31 (9 + 2 * NkBINS)
#define SP_GROW32 (9 + 3 * NkBINS)
#define SP_FOMEGA1 (9 + 4 * NkBINS)
#define SP_FOMEGA2 (9 + 5 * NkBINS)
#define SP_FOMEGA31 (9 + 6 * NkBINS)
#define SP_FOMEGA32 (9 + 7 * NkBINS)

#define SP_PK (9 + 8 * NkBINS)
#define SP_EOS (9 + 8 * NkBINS + 1)
#define SP_INTEOS (9 + 8 * NkBINS + 2)

#define NSPLINES_BASE (9 + 8 * NkBINS + 3)

#ifdef READ_HUBBLE_TABLE
#define SP_EXT_HUBBLE (NSPLINES_BASE)
#define NSPLINES (NSPLINES_BASE + 1)
#else
#define NSPLINES (NSPLINES_BASE)
#endif
```


----- FILE: src/distribute.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "pinocchio.h"

#define BUFLEN 50000

//#define DEBUG

static product_data *comm_buffer;
#ifndef CLASSIC_FRAGMENTATION
static unsigned long long frag_offset;
#endif

int intersection(int *, int *, int *);
int send_data(int *, int);
int recv_data(int *, int);
int keep_data(int *, int *);
unsigned int fft_space_index(unsigned int, int *);
unsigned int subbox_space_index(unsigned int, int *);
#ifndef CLASSIC_FRAGMENTATION
int get_distmap_bit(unsigned int *, unsigned int);
void set_distmap_bit(unsigned int *, unsigned int, int);
void build_distmap(unsigned int *, int *);
void update_distmap(unsigned int *, int *);
#endif

#ifdef DEBUG
FILE *DBGFD;
#endif

int distribute(void)
{
  /* Distributes products from fft-space to sub-volumes */

  int my_fft_box[6],my_subbox[6];
  int log_ntask, bit, receiver, sender;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* this defines the box that the task possesses in the subbox space
     (the starting coordinate may be negative) */
  my_subbox[0] = subbox.stabl[_x_];
  my_subbox[1] = subbox.stabl[_y_];
  my_subbox[2] = subbox.stabl[_z_];
  my_subbox[3] = subbox.Lgwbl[_x_];
  my_subbox[4] = subbox.Lgwbl[_y_];
  my_subbox[5] = subbox.Lgwbl[_z_];

#ifdef DEBUG
  char fname[SBLENGTH];
  sprintf(fname,"Task%d.dbg",ThisTask);
  DBGFD = fopen(fname,"a");
#endif

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* these are the already stored particles */
#ifndef CLASSIC_FRAGMENTATION
  frag_offset=subbox.Nneeded;
#endif

  /* stores the data relative to the intersection of the fft box and subbox */
  if (keep_data(my_fft_box, my_subbox))
    return 1;

  comm_buffer=(product_data*)calloc(BUFLEN , sizeof(product_data));
  if (comm_buffer==0x0)
    {
      printf("ERROR on task %d: could not allocate comm_buffer in distribute\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  /* hypercubic communication scheme */
  for (log_ntask=0; log_ntask<1000; log_ntask++)
    if (1<<log_ntask >= NTasks)
      break;

  /* loop on hypercube dimension */
  for (bit=1; bit<1<<log_ntask; bit++)
    {
    /* loop on tasks */
      for (sender=0; sender<NTasks; sender++)
	{
	  /* receiver task is computed with a bitwise xor */
	  receiver = sender ^ bit;

	  /* condition on sender and receiver */
	  if (receiver < NTasks && sender < receiver)
	    {

	      /* the communication will be done
		 first sender -> receiver and then receiver -> sender */

	      if (ThisTask==sender)
		send_data(my_fft_box, receiver);
	      else if (ThisTask==receiver)
		{
		  if (recv_data(my_subbox, sender))
		    return 1;
		}

	      if (ThisTask==receiver)
		send_data(my_fft_box, sender);
	      else if (ThisTask==sender)
		{
		  if (recv_data(my_subbox, receiver))
		    return 1;
		}

	    }
	}
    }


  /* updates the number of stored particles */
#ifdef CLASSIC_FRAGMENTATION
  subbox.Nstored=subbox.Npart;
#else

  if (frag_offset > subbox.Nalloc)
    subbox.Nstored=subbox.Nalloc;
  else
    subbox.Nstored=frag_offset;

  subbox.Nneeded=frag_offset;

#endif

  free(comm_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
  fclose(DBGFD);
#endif

  return 0;
}


int intersection(int *fbox, int *sbox, int *ibox)
{

  unsigned int istart[3], istop[3], istart2[3], istop2[3], dim, Nint,
    stop1, stop2, this, off;
  unsigned int ISTARTx, ISTOPx, ISTARTy, ISTOPy, ISTARTz, ISTOPz, ax, ay, az;

  /* intersection of the two boxes is considered dimension by dimension */
  Nint=1;
  for (dim=0; dim<3; dim++)
    {
      /* in case, fix negative starting point */
      if (sbox[dim]<0)
	sbox[dim]+=MyGrids[0].GSglobal[dim];

      /* first stopping point for subbox is at most the global box edge */
      stop1=fbox[dim]+fbox[dim+3];
      if (sbox[dim]+sbox[dim+3] > MyGrids[0].GSglobal[dim])
	stop2=MyGrids[0].GSglobal[dim];
      else
	stop2=sbox[dim]+sbox[dim+3];

      /* intersection up to the global box edge */
      istart[dim] = ( fbox[dim] > sbox[dim] ? fbox[dim] : sbox[dim] );
      istop[dim] = ( stop1 < stop2 ? stop1 : stop2 );

      /* if the subbox goes beyond the global box edge,
	 apply PBCs to the other segment and check the intersection */
      if ( (stop2=sbox[dim]+sbox[dim+3]) > MyGrids[0].GSglobal[dim] )
	{
	  stop2 = stop2%MyGrids[0].GSglobal[dim];
	  istart2[dim] = ( fbox[dim] > 0 ? fbox[dim] : 0);
	  istop2[dim] = ( stop1 < stop2 ? stop1 : stop2 );
	}
      else
	{
	  istart2[dim]=1;
	  istop2[dim]=0;
	}

      /* this dimension contributes 0, 1 or 2 */
      Nint *= (istart[dim] < istop[dim]) + (istart2[dim] < istop2[dim]);
    }

  /* store all intersections, looping on the two options for each dimension */
  if (Nint)
    {
      this=0;
      for (ax=0; ax<2; ax++)
	{
	  if (ax)
	    {
	      ISTARTx=istart[0];
	      ISTOPx=istop[0];
	    }
	  else
	    {
	      ISTARTx=istart2[0];
	      ISTOPx=istop2[0];
	    }
	  for (ay=0; ay<2; ay++)
	    {
	      if (ay)
		{
		  ISTARTy=istart[1];
		  ISTOPy=istop[1];
		}
	      else
		{
		  ISTARTy=istart2[1];
		  ISTOPy=istop2[1];
		}
	      for (az=0; az<2; az++)
		{
		  if (az)
		    {
		      ISTARTz=istart[2];
		      ISTOPz=istop[2];
		    }
		  else
		    {
		      ISTARTz=istart2[2];
		      ISTOPz=istop2[2];
		    }

		  if ( (ISTARTx<ISTOPx) & (ISTARTy<ISTOPy) & (ISTARTz<ISTOPz) )
		    {
		      off=this*6;
		      ibox[  off]=ISTARTx;
		      ibox[1+off]=ISTARTy;
		      ibox[2+off]=ISTARTz;
		      ibox[3+off]=ISTOPx-ISTARTx;
		      ibox[4+off]=ISTOPy-ISTARTy;
		      ibox[5+off]=ISTOPz-ISTARTz;
		      this++;
		    }
		}
	    }
	}
    }


#ifdef DEBUG
  fprintf(DBGFD,"INTERSECTION\n");
  fprintf(DBGFD,"Task %d, fft box:        %d %d %d   %d %d %d\n",
	 ThisTask,fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);
  fprintf(DBGFD,"         subvolume:      %d %d %d   %d %d %d\n",
	 sbox[0],sbox[1],sbox[2],sbox[3],sbox[4],sbox[5]);
  for (this=0; this<Nint; this++)
    {
      off=this*6;
      fprintf(DBGFD,"         intersection %d: %d %d %d   %d %d %d\n",
	     this,ibox[off],ibox[1+off],ibox[2+off],ibox[3+off],ibox[4+off],ibox[5+off]);
    }
  if (!Nint)
    fprintf(DBGFD,"        no intersections\n");
#endif

  return Nint;
}


int send_data(int *mybox, int target)
{
  /* This routine sends the content of a box to a target task.
     Communication is divided in these stages:
     1) the sender communicates the start and length of its box,
     2) the sender receives the start and the length of the needed box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the sender receives a map of needed particles
     5) the sender loops on particles and sends them in chunks of size BUFLEN
  */

#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
#endif
  unsigned int size, bufcount;
  int targetbox[6], interbox[48], i, Nint, box, off;
  MPI_Status status;

  MPI_Send(mybox    , 6, MPI_INT, target, 0, MPI_COMM_WORLD);
  MPI_Recv(targetbox, 6, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

  /* up to 8 intersections of the two boxes */
  Nint=intersection(mybox, targetbox, interbox);
  for (box=0; box<Nint; box++)
    {
      off=box*6;
      size = interbox[3+off] * interbox[4+off] * interbox[5+off];

#ifdef DEBUG
      fprintf(DBGFD,"Task %d will send this box: %d %d %d -- %d %d %d\n",
	     ThisTask, interbox[0+off],interbox[1+off],interbox[2+off],
	     interbox[3+off],interbox[4+off],interbox[5+off]);
#endif

#ifndef CLASSIC_FRAGMENTATION
      /* receive the map of needed particles */
      mapl = size/8 + 1;
      map = (unsigned int*)calloc(mapl, sizeof(unsigned int));

      /* the receiver has built its map and is sending it */
      MPI_Recv(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

      /* updates the map and sends it back to the target */
      update_distmap(map, interbox+off);
      MPI_Send(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD);
#endif

      /* load particles on the buffer and send them to the target */
      bufcount=0;
      for (i=0; i<size; i++)
	{
#ifndef CLASSIC_FRAGMENTATION
	  if (get_distmap_bit(map, i))
#endif
	    {
	      memcpy(&comm_buffer[bufcount],&products[fft_space_index(i, interbox+off)],sizeof(product_data));

	      bufcount++;
	    }
	  if (bufcount==BUFLEN)
	    {
#ifdef DEBUG
	      fprintf(DBGFD,"...sending %d products to Task %d... %d\n",bufcount,target,(int)sizeof(product_data));
	      if (bufcount<10)
		{
		  for (int u=0; u<bufcount; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
	      else
		{
		  for (int u=0; u<5; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD," ... ");
		  for (int u=bufcount-5; u<bufcount; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
#endif
	      MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
	      MPI_Send(comm_buffer, bufcount * sizeof(product_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
	      bufcount=0;
	    }
	}

      if (bufcount)
	{
#ifdef DEBUG
	  fprintf(DBGFD,"...sending %d products to Task %d...\n",bufcount,target);
	  if (bufcount<10)
	    {
	      for (int u=0; u<bufcount; u++)
		fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
	      fprintf(DBGFD,"\n");
	    }
	  else
	    {
	      for (int u=0; u<5; u++)
		fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
	      fprintf(DBGFD," ... ");
	      for (int u=bufcount-5; u<bufcount; u++)
```


----- FILE: src/fmax-fftw.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

//#define DEBUG

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif


double greens_function(double *, double, int, int);

int set_one_grid(int ThisGrid)
{
  ptrdiff_t L, M, N;
  ptrdiff_t alloc_local, local_n0, local_0_start;
  ptrdiff_t local_n1, local_1_start;

  fftw_mpi_init();

  MyGrids[ThisGrid].norm = 1.0
    /(double)MyGrids[ThisGrid].GSglobal_x
    /(double)MyGrids[ThisGrid].GSglobal_y
    /(double)MyGrids[ThisGrid].GSglobal_z;
  MyGrids[ThisGrid].CellSize = MyGrids[ThisGrid].BoxSize/(double)MyGrids[ThisGrid].GSglobal_x;

  L=MyGrids[ThisGrid].GSglobal_x;
  M=MyGrids[ThisGrid].GSglobal_y;
  N=MyGrids[ThisGrid].GSglobal_z;

  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_3d_transposed(L, M, N/2+1, MPI_COMM_WORLD,
						  &local_n0, &local_0_start,
						  &local_n1, &local_1_start);

  MyGrids[ThisGrid].GSlocal_x = MyGrids[ThisGrid].GSglobal_x;
  MyGrids[ThisGrid].GSlocal_y = MyGrids[ThisGrid].GSglobal_y;
  MyGrids[ThisGrid].GSlocal_z = local_n0;

  MyGrids[ThisGrid].GSstart_x = 0;
  MyGrids[ThisGrid].GSstart_y = 0;
  MyGrids[ThisGrid].GSstart_z = local_0_start;

  MyGrids[ThisGrid].GSlocal_k_x = MyGrids[ThisGrid].GSglobal_x;
  MyGrids[ThisGrid].GSlocal_k_y = local_n1;
  MyGrids[ThisGrid].GSlocal_k_z = MyGrids[ThisGrid].GSglobal_z;

  MyGrids[ThisGrid].GSstart_k_x = 0;
  MyGrids[ThisGrid].GSstart_k_y = local_1_start;
  MyGrids[ThisGrid].GSstart_k_z = 0;

  MyGrids[ThisGrid].total_local_size_fft = 2*alloc_local;
  MyGrids[ThisGrid].total_local_size = MyGrids[ThisGrid].GSlocal_x * MyGrids[ThisGrid].GSlocal_y * MyGrids[ThisGrid].GSlocal_z;

  /* this gives is the number of unused x-grid points in rvector_fft */
  if (MyGrids[0].GSglobal_x%2)
    MyGrids[0].off=1;
  else
    MyGrids[0].off=2;

  return 0;
}

int initialize_fft()
{
  ptrdiff_t L, M, N;
  int igrid;

  for (igrid=0; igrid<Ngrids; igrid++)
    {

      L=MyGrids[igrid].GSglobal_x;
      M=MyGrids[igrid].GSglobal_y;
      N=MyGrids[igrid].GSglobal_z;

      /* create plan for out-of-place r2c DFT */
      MyGrids[igrid].forward_plan = fftw_mpi_plan_dft_r2c_3d(L, M, N, rvector_fft[igrid], cvector_fft[igrid],
							     MPI_COMM_WORLD, FFTW_MEASURE + FFTW_MPI_TRANSPOSED_OUT);

      MyGrids[igrid].reverse_plan = fftw_mpi_plan_dft_c2r_3d(L, M, N, cvector_fft[igrid], rvector_fft[igrid],
							     MPI_COMM_WORLD, FFTW_MEASURE + FFTW_MPI_TRANSPOSED_IN);
    }

  return 0;
}


double forward_transform(int ThisGrid)
{
  double time;

  time=MPI_Wtime();

  fftw_execute(MyGrids[ThisGrid].forward_plan);

  return MPI_Wtime()-time;
}


double reverse_transform(int ThisGrid)
{

  int i;
  double time;
  // double ave, var; //levare

  time=MPI_Wtime();

  fftw_execute(MyGrids[ThisGrid].reverse_plan);

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
      rvector_fft[ThisGrid][i]*=MyGrids[ThisGrid].norm;

#ifdef DEBUG
  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    fprintf(results, " %d  %g \n", i, rvector_fft[ThisGrid][i]);
  fclose(results);
#endif

  return MPI_Wtime()-time;
}


int finalize_fftw()
{
  int igrid;

  for (igrid=Ngrids-1; igrid>=0; igrid--)
    {
      fftw_destroy_plan(MyGrids[igrid].forward_plan);
      fftw_destroy_plan(MyGrids[igrid].reverse_plan);
    }

  for (igrid=Ngrids-1; igrid>=0; igrid--)
    if (deallocate_fft_vectors(igrid))
      return 1;

  fftw_cleanup();

  return 0;
}



int compute_derivative(int ThisGrid, int first_derivative, int second_derivative)
{
  int swap, local_x,local_y,local_z,ixx,iyy,izz,index,nxhalf,nyhalf,nzhalf;
  double kx,ky,kz,kxnorm,kynorm,kznorm,green,k_squared,smoothing,tmp,time;
  double diff_comp[4];


#ifdef DEBUG
  sprintf(filename,"results.%d-%d.%d",first_derivative,second_derivative,ThisTask);
  results=fopen(filename,"w");
#endif

  /* k vectors */
  kxnorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_x;
  kynorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_y;
  kznorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_z;

  /* Nyquist frequencies */
  nxhalf = MyGrids[ThisGrid].GSglobal_x/2;
  nyhalf = MyGrids[ThisGrid].GSglobal_y/2;
  nzhalf = MyGrids[ThisGrid].GSglobal_z/2;

  /* for first derivatives the real and imaginary parts must be swapped */
  swap=((first_derivative==0 && second_derivative> 0) ||
	(first_derivative> 0 && second_derivative==0));

/*
  NB: ix, iy and iz are global coordinates of the box (from 0 to N-1)
      ix: [0,N/2]
      iy, iz: [0,N-1]
      iyy and izz are unfolded, they run from -N/2+1 to N/2 (ixx = ix)
      local_x, local_y and local_z are local coordinates of the slice
      (for the x and z coordinates they are the same as ix and iz)
*/


/* loop over k-space indices */
/* This loop is correct for transposed order in FFTW */

  for (local_z = 0; local_z < MyGrids[ThisGrid].GSlocal_k_z; local_z++)
    {
      izz = local_z + MyGrids[ThisGrid].GSstart_k_z;
      if (local_z > nzhalf)
	izz -= MyGrids[ThisGrid].GSglobal_z;
      kz  = kznorm*izz;

      for (local_y = 0; local_y < MyGrids[ThisGrid].GSlocal_k_y; local_y++)
	{
	  iyy = local_y + MyGrids[ThisGrid].GSstart_k_y;
	  if (iyy > nyhalf)
	    iyy -= MyGrids[ThisGrid].GSglobal_y;
	  ky  = kynorm*iyy;

	  for (local_x = 0; local_x <= nxhalf; local_x++)
	    {
	      ixx = local_x;
	      kx  = kxnorm*ixx;

              k_squared  = kx*kx + ky*ky + kz*kz;

	      /* corresponding index of real part in vector (imaginary in index + 1) */
	      index = 1 + 2*local_x + (MyGrids[ThisGrid].GSglobal_x+MyGrids[ThisGrid].off)
		*(local_z + local_y* MyGrids[ThisGrid].GSglobal_z);

              if (k_squared != 0.)
		{

		  /* Gaussian smoothing window */
		  smoothing = exp(-0.5 * k_squared * Rsmooth * Rsmooth);

		  /* the components are stored in the vectors that control the differentiation */
		  diff_comp[0] = 1.0;
		  diff_comp[1] = kx;
		  diff_comp[2] = ky;
		  diff_comp[3] = kz;

		  green = greens_function(diff_comp, k_squared, first_derivative, second_derivative);

		  (cvector_fft[ThisGrid][index/2])[0] *=  green * smoothing;
		  (cvector_fft[ThisGrid][index/2])[1] *=  green * smoothing;
		}

              if (swap)
		{
		  tmp                                 = (cvector_fft[ThisGrid][index/2])[1];
		  (cvector_fft[ThisGrid][index/2])[1] = (cvector_fft[ThisGrid][index/2])[0];
		  (cvector_fft[ThisGrid][index/2])[0] = -tmp;
		}
	    }
	}
    }

  if (!ThisTask)
    printf("[%s] compute_derivative: starting fft\n",fdate());

  time=reverse_transform(ThisGrid);

  if (!ThisTask)
    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);

  cputime.fft+=time;

  return 0;
}


double greens_function(double *diff_comp, double k_squared, int first_derivative, int second_derivative)
{

  /* in this case the greens_function is simply 1 */
  if (first_derivative == -1 && second_derivative == -1)
    return 1.0;

  if (first_derivative==0 && second_derivative==0)
    return -diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;
  else
    return diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;

}


void write_in_cvector(int ThisGrid, double *vector)
{
  int i;

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    *((double*)cvector_fft[ThisGrid]+i)=*(vector+i);

}


void write_from_cvector(int ThisGrid, double *vector)
{
  int i;

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    *(vector+i)=*((double*)cvector_fft[ThisGrid]+i);

}


void write_in_rvector(int ThisGrid, double *vector)
{
  int lx,ly,lz;

  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) =
	    *(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y));

  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=MyGrids[ThisGrid].GSlocal_x; lx<MyGrids[ThisGrid].GSlocal_x+MyGrids[ThisGrid].off; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) = 0.0;

}


void write_from_rvector(int ThisGrid, double *vector)
{
  int lx,ly,lz;

  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y)) =
          *(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y));

}
```


----- FILE: src/fmax-pfft.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"


// -------------------------
// defines and variables
// -------------------------

//#define DEBUG

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif

#define GRID MyGrids[ThisGrid]

typedef ptrdiff_t point[3];
point *starts;

// -------------------------
// functions' prototypes
// -------------------------


double greens_function                  (double *, double, int, int);
int    cubes_order                      (const void *, const void *);

// -------------------------
// code segment
// -------------------------


int cubes_order(const void *A, const void *B)
{
  float a = *(starts[*(int*)A]);
  float b = *(starts[*(int*)B]);

  return (a - b);

  /*
  if( a > b)
    return 1;
  if(a < b)
    return -1;
  return 0;
  */
}




int set_one_grid(int ThisGrid)
{
  ptrdiff_t    alloc_local;
  unsigned int pfft_flags;

  GRID.norm = (double)1.0 /
    ((double)GRID.Ntotal);

  GRID.CellSize = (double)GRID.BoxSize / GRID.GSglobal[_x_];

  //pfft_flags  = PFFT_MEASURE;
  pfft_flags  = 0;
  if(params.use_transposed_fft)
    pfft_flags  |= PFFT_TRANSPOSED_OUT;

  alloc_local = pfft_local_size_dft_r2c_3d(GRID.GSglobal, FFT_Comm, pfft_flags,
					   GRID.GSlocal, GRID.GSstart,
					   GRID.GSlocal_k, GRID.GSstart_k);


  dprintf(VDBG, ThisTask, "[set grid %02d] task %d %ld "
	  "i: %ld %ld %ld - i start: %ld %ld %ld - "
	  "o: %ld %ld %ld - o start: %ld %ld %ld\n",
	  ThisGrid, ThisTask, alloc_local,
	  GRID.GSlocal[_x_], GRID.GSlocal[_y_], GRID.GSlocal[_z_],
	  GRID.GSstart[_x_], GRID.GSstart[_y_], GRID.GSstart[_z_],
	  GRID.GSlocal_k[_x_], GRID.GSlocal_k[_y_], GRID.GSlocal_k[_z_],
	  GRID.GSstart_k[_x_], GRID.GSstart_k[_y_], GRID.GSstart_k[_z_]);


  GRID.total_local_size_fft = 2 * alloc_local;
  GRID.total_local_size     = GRID.GSlocal[_x_] * GRID.GSlocal[_y_] * GRID.GSlocal[_z_];

  MyGrids[0].off = 0;
  // order the sub-blocks by row-major order (i, j, k), k first, then j, then i

  /* int           i; */
  /* unsigned int index; */

  /* starts         = (point*)malloc(sizeof(point) * NTasks); */

  /* MPI_Allgather(GRID.GSlocal, sizeof(point), MPI_BYTE, starts, sizeof(point), MPI_BYTE, MPI_COMM_WORLD); */
  /* for(i = 0; i < NTasks; i++) */
  /*   { */
  /*     index            = (starts[i][_x_]*GRID.GSglobal[_y_] + starts[i][_y_])*GRID.GSglobal[_z_] + starts[i][_z_]; */
  /*     starts[i][_x_]   = index; */
  /*     cubes_ordering[i] = i; */
  /*   } */

  /* qsort(cubes_ordering, NTasks, sizeof(int), cubes_order); */

  /* free(starts); */

  return 0;
}




int compute_fft_plans()
{
  ptrdiff_t DIM[3];
  int       pfft_flags;
  int       ThisGrid;

  dprintf(VMSG, 0, "[%s] Computing fft plans\n",fdate());

  for (ThisGrid = 0; ThisGrid < Ngrids; ThisGrid++)
    {

      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];

      /* create plan for out-of-place DFT */

      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE ;
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_OUT;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
#endif
      GRID.forward_plan = pfft_plan_dft_r2c_3d(DIM, rvector_fft[ThisGrid], cvector_fft[ThisGrid],
					       FFT_Comm, PFFT_FORWARD, pfft_flags);


      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];

      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE;
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_IN;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
#endif
      GRID.reverse_plan = pfft_plan_dft_c2r_3d(DIM, cvector_fft[ThisGrid], rvector_fft[ThisGrid],
					       FFT_Comm, PFFT_BACKWARD, pfft_flags);
    }


  dprintf(VMSG, 0, "[%s] fft plans done\n",fdate());

  return 0;
}


double forward_transform(int ThisGrid)
{
  double time;

  time=MPI_Wtime();

  pfft_execute(GRID.forward_plan);

  return MPI_Wtime()-time;
}


double reverse_transform(int ThisGrid)
{

  int i;
  double time;

  time=MPI_Wtime();

  pfft_execute(GRID.reverse_plan);

//  dvec         NORM   = {GRID.norm, GRID.norm, GRID.norm, GRID.norm};
//  unsigned int mysize = GRID.total_local_size_fft / 4;

// #pragma GCC ivdep
//   for (i = 0; i < mysize; i++)
//     rvector_fft[ThisGrid][i] *= NORM;

  for (i = GRID.total_local_size_fft - GRID.total_local_size_fft%4 ; i < GRID.total_local_size_fft; i++)
    rvector_fft[ThisGrid][i] *= GRID.norm;

  // non-vector code
  for (i = 0 ; i < GRID.total_local_size_fft; i++)
      rvector_fft[ThisGrid][i] *= GRID.norm;

  return MPI_Wtime() - time;
}


int finalize_fft()
{
#ifndef RECOMPUTE_DISPLACEMENTS
  int ThisGrid;

  for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--)
    {
      pfft_destroy_plan(GRID.forward_plan);
    }

  /* for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--) */
  /*   { */
  /*     pfft_free(GRID.forward_plan); */
  /*     pfft_free(GRID.reverse_plan); */
  /*   } */

  pfft_cleanup();
  MPI_Comm_free(&FFT_Comm);
#endif

  return 0;
}


int compute_derivative(int ThisGrid, int first_derivative, int second_derivative)
{
  int    swap, local[3], start[3], C[3], N[3], Nhalf[3];
  double knorm[3];

#ifdef DEBUG
  sprintf(filename,"results.%d-%d.%d",first_derivative,second_derivative,ThisTask);
  results=fopen(filename,"w");
#endif

  // note: for PFFT the transposed order is different when
  // the subdivision is done either in 1 or higher dimensions.
  // In 1D the memory order is [ y, x, z ], while in 2D and 3D
  // it is [ y, z, x].


  if(!params.use_transposed_fft)
    // this is the non-transposed order x, y, z
    for(int i = 0; i < 3; i++)
      C[i] = i;
  else
    {
      if(internal.tasks_subdivision_dim > 1)
	C[0] = _y_, C[1] = _z_, C[2] = _x_;
      else
	C[0] = _y_, C[1] = _x_, C[2] = _z_;
    }

  for(int i = 0; i < 3; i++)
    {
      // grid number
      N[i]     = GRID.GSglobal[C[i]];
      // Nyquist frequencies
      Nhalf[i] = N[i] / 2;
      // k vectors
      knorm[i] = 2.*PI / (double)GRID.GSglobal[ C[i] ];

      local[i] = GRID.GSlocal_k[C[i]];
      start[i] = GRID.GSstart_k[C[i]];

    }

  // for first derivatives the real and imaginary parts must be swapped
  swap=((first_derivative==0 && second_derivative> 0) ||
	(first_derivative> 0 && second_derivative==0));

  double growth_rate=1.0;

  // loop over k-space indices

  int idx;
  for (idx = 0; idx < local[_x_]; idx++)
    {
      int    ii[3];

      ii[_x_] = idx + start[_x_]; // this is the kx index

      if (ii[_x_] > Nhalf[_x_])   // kx-N if kx > Nx/2
	ii[_x_] -= N[_x_];

      double k_x  = knorm[_x_] * ii[_x_];
      double k2_0 = k_x * k_x;

      int idy;
      for (idy = 0; idy < local[_y_]; idy++)
	{
	  ii[_y_] = idy + start[_y_]; // this is the ky index

	  if (ii[_y_] > Nhalf[_y_])
	    ii[_y_] -= N[_y_];        // ky-N if ky > Ny/2

	  double k_y  = knorm[_y_] * ii[_y_];
	  double k2_1 = k2_0 + k_y * k_y;

	  int idz;
	  for (idz = 0; idz < local[_z_]; idz++)
	    {
	      ii[_z_] = idz + start[_z_]; // this is the kz index

	      if (ii[_z_] > Nhalf[_z_])
		ii[_z_] -= N[_z_];        // kz-N if ky > Nz/2

	      double k_z  = knorm[_z_] * ii[_z_];

              double k_squared  = k2_1 + k_z * k_z;
	      double k_module = sqrt(k_squared);

	      /* In the scale-dependent case the delta(k) must be multiplied
		 by the relevant growth rate */
	      switch (ScaleDep.order)
		{
		case 0:
		  growth_rate = 1.0;
		  break;
		case 1:
		  growth_rate = GrowingMode(ScaleDep.redshift,k_module);
		  break;
		case 2:
		  growth_rate = GrowingMode_2LPT(ScaleDep.redshift,k_module);
		  break;
		case 3:
		  growth_rate = GrowingMode_3LPT_1(ScaleDep.redshift,k_module);
		  break;
		case 4:
		  growth_rate = GrowingMode_3LPT_2(ScaleDep.redshift,k_module);
		  break;
		default:
		  growth_rate = 1.0;
		  break;
		}

	      int index = 2*(( idx * local[_y_] + idy ) * local[_z_] + idz);

              if (k_squared != 0.)
		{

		  // Gaussian smoothing window
		  double smoothing = exp(-0.5 * k_squared * Rsmooth * Rsmooth);

		  double diff_comp[4];

		  // the components are stored in the vectors that control the differentiation

		  diff_comp[0] = 1.0;
		  diff_comp[1] = k_x;
		  diff_comp[2] = k_y;
		  diff_comp[3] = k_z;

		  double green = greens_function(diff_comp, k_squared, first_derivative, second_derivative);

		  (cvector_fft[ThisGrid][index/2])[0] *=  green * smoothing * growth_rate;
		  (cvector_fft[ThisGrid][index/2])[1] *=  green * smoothing * growth_rate;
		}

              if (swap)
		{
		  double tmp                          = (cvector_fft[ThisGrid][index/2])[1];
		  (cvector_fft[ThisGrid][index/2])[1] = (cvector_fft[ThisGrid][index/2])[0];
		  (cvector_fft[ThisGrid][index/2])[0] = -tmp;
		}
	    }
	}
    }

  /* int idy, idz, index; */
  /* for (idx = 0; idx < local[_x_]; idx++) */
```


----- FILE: src/fmax.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

int compute_second_derivatives(double, int);
int Fmax_PDF(void);


/* Computation of collapse times and displacements */
int compute_fmax(void)
{
  int ismooth,ThisGrid;
  double cputmp,cpusm;

  cputime.fmax=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] First part: computation of collapse times\n",fdate());

  ThisGrid=0;

  /*
    TO EXTEND TO MULTIPLE GRIDS here we should

    1) initialize fft for the large-scale grid (the same as the small-scale one?)
    2) compute first derivatives of the large grid
    3) compute second derivatives of the large grid
    4) distribute to all processors and store the few used grid points
    5) re-initialize fft for the small-scale grid
  */

  ScaleDep.order=0; /* collapse times must be computed as for LambdaCDM */
  ScaleDep.redshift=0.0;


  /****************************
   * CYCLE ON SMOOTHING RADII *
   ****************************/

  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {
      if (!ThisTask)
	printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n",
	       fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
	       sqrt(Smoothing.Variance[ismooth]) );

      cpusm=MPI_Wtime();

      /*
	 Part 1:
	 Compute second derivatives of the potential
      */

      if (!ThisTask)
	printf("[%s] Computing second derivatives\n",fdate());

      cputmp=MPI_Wtime();

      if (compute_second_derivatives(Smoothing.Radius[ismooth] ,ThisGrid))
	return 1;

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done second derivatives, cpu time = %f s\n",fdate(),cputmp);
      cputime.deriv += cputmp;

      /*
	 Part 2:
	 Compute collapse times
      */

      if (!ThisTask)
	printf("[%s] Computing collapse times\n",fdate());

      cputmp=MPI_Wtime();

#ifdef TABULATED_CT
      /* initialize spline for interpolating collapse times */
      if (initialize_collapse_times(ismooth,0))
	return 1;

      cputmp=MPI_Wtime()-cputmp;

      if (!ThisTask)
	{
	  if (strcmp(params.CTtableFile,"none"))
	    printf("[%s] Collapse times read from file %s\n",fdate(),params.CTtableFile);
	  else
	    printf("[%s] Collapse times computed for interpolation, cpu time =%f s\n",fdate(),cputmp);
	}

      cputime.coll+=cputmp;
      cputmp=MPI_Wtime();
#endif

      if (compute_collapse_times(ismooth))
	return 1;

#ifdef TABULATED_CT
      /* this is needed only for debug options, to be removed in the official code */
      if (reset_collapse_times(ismooth))
	return 1;
#endif

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing collapse times, cpu time = %f s\n",fdate(),cputmp);
      cputime.coll+=cputmp;

      /*
	End of cycle on smoothing radii
      */

      cpusm = MPI_Wtime()-cpusm;

      if (!ThisTask)
	printf("[%s] Completed, R=%6.3f, expected sigma: %7.4f, computed sigma: %7.4f, cpu time = %f s\n",fdate(),
	       Smoothing.Radius[ismooth], sqrt(Smoothing.Variance[ismooth]),
	       sqrt(Smoothing.TrueVariance[ismooth]),
	       cpusm );

      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  /********************************
   * COMPUTATION OF DISPLACEMENTS *
   ********************************/

  /* computes the displacements for all particles at zero smoothing,
     assuming second derivatives are already in place
     displacements are computed to the redshift of the first (or only) segment
  */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("\n[%s] Computing displacements  for redshift %f\n",fdate(),ScaleDep.z[0]);

  if (compute_displacements(1,0,ScaleDep.z[0]))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done computing displacements, cpu time = %f s\n",fdate(),cputmp);

  /***************/
  /* END OF FMAX */
  /***************/

  if (finalize_fft())
    return 1;

  if (Fmax_PDF())
    return 1;

  cputime.fmax = MPI_Wtime() - cputime.fmax;
  if (!ThisTask)
    printf("[%s] Finishing fmax, total fmax cpu time = %14.6f\n"
	   "\t\t IO       : %14.6f (%14.6f total time without I/O)\n"
	   "\t\t FFT      : %14.6f\n"
	   "\t\t COLLAPSE : %14.6f\n",
	   fdate(), cputime.fmax, cputime.io, cputime.fmax-cputime.io, cputime.fft, cputime.coll);

  return 0;
}


int compute_first_derivatives(double R, int ThisGrid, int order, double* vector)
{
  /* computes second derivatives of the potential */

  double timetmp = 0;

  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for (int ia=1; ia<=3; ia++)
    {

      if (!ThisTask)
	printf("[%s] Computing 1st derivative: %d\n",fdate(),ia);

      double tmp = MPI_Wtime();
      write_in_cvector(ThisGrid, vector);
      timetmp += MPI_Wtime() - tmp;

      if (compute_derivative(ThisGrid,ia,0))
	return 1;

      tmp = MPI_Wtime();
      write_from_rvector_to_products(ThisGrid, ia-1, order);
      timetmp += MPI_Wtime() - tmp;
    }

  cputime.mem_transf += timetmp;
  return 0;
}


int compute_second_derivatives(double R, int ThisGrid)
{

  /* computes second derivatives of the potential */

  double timetmp = 0;

  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for ( int ia = 1; ia <= 3; ia++ )
    for ( int ib = ia; ib <= 3; ib++ )
      {

	int ider=( ia == ib ? ia : ia+ib+1 );

	if (!ThisTask)
	  printf("[%s] Computing 2nd derivative: %d\n",fdate(),ider);

	double tmp = MPI_Wtime();
	write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	timetmp += MPI_Wtime() - tmp;

	if (compute_derivative(ThisGrid,ia,ib))
	  return 1;

	tmp = MPI_Wtime();
	write_from_rvector(ThisGrid, second_derivatives[ThisGrid][ider-1]);
	timetmp += MPI_Wtime() - tmp;
      }

  cputime.mem_transf += timetmp;
  return 0;
}


char *fdate()
{
  /* returns a 24-char string with the full date and time
     (the year is moved at the middle of the string) */

  time_t current_time;
  char *string;
  int n;

  /* format from ctime:
    0123456789
              0123456789
	                0123
    Www Mmm dd hh:mm:ss yyyy
  */

  current_time=time(NULL);
  string=ctime(&current_time);

  for (n=0; n<10; n++)
    *(date_string+n)=*(string+n);
  for (n=10; n<15; n++)
    *(date_string+n)=*(string+n+9);
  for (n=10; n<19; n++)
    *(date_string+n+5)=*(string+n);
  *(date_string+24)='\0';

  return date_string;
}


int compute_displacements(int compute_sources, int recompute_sd, double redshift)
{
  /* This routine computes the displacement fields at a specified redshift
     if recompute_sd is true, it recomputes the second derivatives */

  double cputmp;

#ifdef TWO_LPT

  if (recompute_sd)
    {
      /* second derivatives are needed for LPT displacements, if not already in place they must be recomputed */
      cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Computing second derivatives\n",fdate());

      ScaleDep.order=0;  /* sources are computed as for LambdaCDM at z=0 */
      ScaleDep.redshift=0.0;

      if (compute_second_derivatives(0.0, 0))
	return 1;

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done second derivatives, cpu time = %f s\n",fdate(),cputmp);
      cputime.deriv+=cputmp;
    }

  /* computes the 2LPT and 3LPT displacement fields and stores them in the products */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("\n[%s] Computing LPT displacements\n",fdate());

  if (compute_LPT_displacements(compute_sources, redshift))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done LPT displacements, cpu time = %f s\n",fdate(),cputmp);
  cputime.lpt+=cputmp;

#endif

  /* computes Zeldovich displacements */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("[%s] Computing first derivatives\n",fdate());

  /* for scale-dependent growth we compute the displacements for the first redshift segment
     in case there is only one segment, this implies to compute them at the final redshift */
  ScaleDep.redshift=redshift;
  ScaleDep.order=1;   /* here we need the first-order growth */

  if (compute_first_derivatives(0.0, 0, 1, kdensity[0]))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done first derivatives, cpu time = %f s\n",fdate(),cputmp);
  cputime.deriv+=cputmp;

  /* /\* Store Zeldovich displacements in the products *\/ */
  /* cputmp=MPI_Wtime(); */
  /* if (!ThisTask) */
  /*   printf("[%s] Storing velocities\n",fdate()); */

  /* if (store_velocities()) */
  /*   return 1; */

  /* cputmp=MPI_Wtime()-cputmp; */
  /* if (!ThisTask) */
  /*   printf("[%s] Done storing velocities, cpu time = %f s\n",fdate(),cputmp); */
  /* cputime.vel+=cputmp; */

  return 0;
}


#include <sys/stat.h>

int dump_products()
{
  /* dumps fmax products on files to reuse them */

  struct stat dr;
  FILE *file;
  char fname[LBLENGTH];

  /* Task 0 checks that the DumpProducts directory exists, or creates it
     and writes a summary file to check that one does not start from a different run */
  if (!ThisTask)
    {
      if(stat(params.DumpDir,&dr))
	  {
	    printf("Creating directory %s\n",params.DumpDir);
	    if (mkdir(params.DumpDir,0755))
	      {
		printf("ERROR IN CREATING DIRECTORY %s (task 0)\n",params.DumpDir);
		return 1;
	      }
	  }

      sprintf(fname,"%ssummary",params.DumpDir);
      file=fopen(fname,"w");
      fprintf(file,"%d   # NTasks\n",NTasks);
      fprintf(file,"%d   # random seed\n",params.RandomSeed);
      fprintf(file,"%d   # grid size\n",params.GridSize[0]);
      fprintf(file,"%d   # length of product_data\n",(int)sizeof(product_data));
      fclose(file);
```


----- FILE: src/fragment.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"
#include <sys/types.h>
#include <sys/stat.h>

// #define DEBUG

int fragment(void);
int count_peaks(int *);
void set_fragment_parameters(int);
void write_map(int);
int create_map(void);
void reorder(int *, int);
void reorder_nofrag(int *, int);
void sort_and_organize(void);
#ifdef RECOMPUTE_DISPLACEMENTS
int shift_all_displacements(void);
int recompute_group_velocities(void);
#endif

void set_fragment_parameters(int order)
{

#ifndef TWO_LPT
  order = 1;
#else
#ifndef THREE_LPT
  if (order > 2)
    order = 2;
#else
  if (order > 3)
    order = 3;
#endif
#endif

  /* Parameters for the simulation */

  f_200 = 0.171;
  switch (order)
  {
  case 1:
#ifdef USE_SIM_PARAMS
    f_m = f_a = 0.495;
    f_rm = -0.075;
    espo = 0.852;
    f_ra = 0.500;
#else
    f_m = f_a = 0.505;
    f_rm = 0.000;
    espo = 0.820;
    f_ra = 0.300;
#endif
    sigmaD0 = 1.7;
    break;

  case 2:
#ifdef USE_SIM_PARAMS
    f_m = f_a = 0.475;
    f_rm = -0.020;
    espo = 0.780;
    f_ra = 0.650;
#else
    f_m = f_a = 0.501;
    f_rm = 0.052;
    espo = 0.745;
    f_ra = 0.334;
#endif
    sigmaD0 = 1.5;
    break;

  case 3:
#ifdef USE_SIM_PARAMS
    f_m = f_a = 0.455;
    f_rm = 0.000;
    espo = 0.755;
    f_ra = 0.700;
#else
    f_m = f_a = 0.5024;
    f_rm = 0.1475;
    espo = 0.6852;
    f_ra = 0.4584;
#endif
    sigmaD0 = 1.2;
    break;

  default:
    break;
  }
}

static int index_compare_F(const void *a, const void *b)
{
  if (frag[*((int *)a)].Fmax == frag[*((int *)b)].Fmax)
    return 0;
  else if (frag[*((int *)a)].Fmax > frag[*((int *)b)].Fmax)
    return -1;
  else
    return 1;
}

static int index_compare_P(const void *a, const void *b)
{
  if (frag_pos[*((int *)a)] == frag_pos[*((int *)b)])
    return 0;
  else if (frag_pos[*((int *)a)] < frag_pos[*((int *)b)])
    return -1;
  else
    return 1;
}

int fragment_driver()
{
  /* This routine can be edited to call fragmentation for many combinations
     of fragment parameters */

  if (!ThisTask)
    printf("\n[%s] Second part: fragmentation of the collapsed medium\n", fdate());

  /* parameters are assigned in the standard way */
  set_fragment_parameters(ORDER_FOR_GROUPS);

  /* reallocates the memory */
  if (reallocate_memory_for_fragmentation())
    return 1;

  if (fragment())
    return 1;

  return 0;
}

int fragment()
{
  /* This function is the driver for fragmentation of collapsed medium
     and construction of halo catalogs */

  int Npeaks, Ngood;
#ifndef CLASSIC_FRAGMENTATION
  unsigned int nadd[2];
#endif
  double BestPredPeakFactor;
  unsigned long long mynadd[2], nadd_all[2];
  double tmp;

  /* timing */
  cputime.frag = MPI_Wtime();
  cputime.group = cputime.sort = cputime.distr = 0.0;
#ifdef PLC
  cputime.plc = 0.0;
#endif

  /* In the updated code, fragmentation is performed twice;
     in the first turn a quick fragmentation is performed to locate
     halos and determine the interesting particles; velocities are not updated;
     in the second turn the full fragmentation is performed, segmentation is applied
     and velocities are updated.
     For the classic fragmentation the first turn is skipped */

#ifdef CLASSIC_FRAGMENTATION
#define START_TURN 1
#else
#define START_TURN 0
  int all_pbc = (subbox.pbc[_x_] && subbox.pbc[_y_] && subbox.pbc[_z_]);
#endif

  for (int turn = START_TURN; turn < 2; turn++)
  {

#ifndef CLASSIC_FRAGMENTATION
    if (!turn)
    {

      /* it creates the map by setting to 1 the well-resolved region */
      if (!ThisTask)
        printf("[%s] Creating map of needed particles\n", fdate());

      if (create_map())
        return 1;

      /* the rest of the first turn is skipped if there are PBCs in all directions */
      if (all_pbc)
        continue;
    }
    else if (!all_pbc)
    {

      /* it updates the map by adding spheres around halos in the boundary layer */
      if (!ThisTask)
        printf("[%s] Updating map of needed particles\n", fdate());

      update_map(nadd);

      mynadd[0] = nadd[0];
      mynadd[1] = nadd[1];
      MPI_Reduce(mynadd, nadd_all, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

      if (!ThisTask)
      {
        printf("[%s] Requesting %Ld particles from the boundary layer\n", fdate(), nadd_all[0]);
        if (nadd_all[1])
          printf("WARNING: %Ld requested particles lie beyond the boundary layer, some halos may be inaccurate\n", nadd_all[1]);
      }
    }

    /* redistribution of products for queried particles */
    if (!ThisTask)
      printf("[%s] Starting %s re-distribution of products\n", fdate(), (turn ? "second" : "first"));
#else

    if (!ThisTask)
      printf("[%s] Starting re-distribution of products\n", fdate());

#endif

    tmp = MPI_Wtime();

    /* distribution of products from the fft space to the subbox space */
    map_to_be_used = 0;
    if (!turn)
      subbox.Nneeded = 0;
    if (distribute())
      return 1;

    tmp = MPI_Wtime() - tmp;
    cputime.distr += tmp;

    if (!ThisTask)
      printf("[%s] Re-distribution of products done, cputime = %14.6f\n", fdate(), tmp);

#ifndef CLASSIC_FRAGMENTATION
    if (subbox.Nneeded > subbox.Nalloc)
    {
      if (turn)
      {
        printf("CRITICAL WARNING in Task %d: it allocated %d but received %d particles\n",
               ThisTask, subbox.Nalloc, subbox.Nneeded);
        printf("The required overhead is %f\n", (float)subbox.Nneeded / (float)MyGrids[0].ParticlesPerTask);
        printf("Please increase MaxMemPerParticle by at least %d and start again\n",
               (int)((float)(subbox.Nneeded - subbox.Nalloc) / (float)MyGrids[0].ParticlesPerTask *
                     (sizeof(product_data) + FRAGFIELDS * sizeof(int))) +
                   1 + (int)params.MaxMemPerParticle);
        if (params.ExitIfExtraParticles)
          return 1;
      }
      else
      {
        printf("ERROR in Task %d: the number of allocated particles (%d) is WAY too small! I need at least %d\n",
               ThisTask, subbox.Nalloc, subbox.Nneeded);
        printf("The required overhead is %f\n", (float)subbox.Nneeded / (float)MyGrids[0].ParticlesPerTask);
        printf("Please increase MaxMemPerParticle to at least %d and start again\n",
               (int)((float)(subbox.Nneeded - subbox.Nalloc) / (float)MyGrids[0].ParticlesPerTask *
                     (sizeof(product_data) + FRAGFIELDS * sizeof(int))) +
                   1 + (int)params.MaxMemPerParticle);
        return 1;
      }
    }

    /* memory checks on the redistributed particles */
    mynadd[0] = subbox.Nstored;
    MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!ThisTask)
      printf("[%s] %s re-distribution of Fmax done, %Ld particles stored by all tasks, average overhead: %f, cputime = %14.6f\n",
             fdate(), (turn ? "Second" : "First"), nadd_all[0],
             (float)nadd_all[0] / (float)MyGrids[0].Ntotal, tmp);

    mynadd[0] = subbox.Nstored;
    MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(mynadd, nadd_all + 1, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!ThisTask)
      printf("[%s] Smallest and largest overhead: %f, %f\n", fdate(),
             (float)nadd_all[0] / (float)MyGrids[0].ParticlesPerTask,
             (float)nadd_all[1] / (float)MyGrids[0].ParticlesPerTask);

    /* sets or updates the map of potentially loaded particles */
    if (!turn)
      for (int i = 0; i < subbox.maplength; i++)
        frag_map[i] = frag_map_update[i];
    else
      for (int i = 0; i < subbox.maplength; i++)
        frag_map[i] |= frag_map_update[i];

#ifdef DEBUG
    write_map(turn);
#endif

    /* this part sorts particles and reorganizes index arrays */
    sort_and_organize();
#else

    /* sorting of particles according to their collapse time */
    tmp = MPI_Wtime();
    if (!ThisTask)
      printf("[%s] Starting sorting\n", fdate());

    for (int i = 0; i < subbox.Npart; i++)
      *(indices + i) = i;
    qsort((void *)indices, subbox.Npart, sizeof(int), index_compare_F);

    tmp = MPI_Wtime() - tmp;
    cputime.sort += tmp;
    if (!ThisTask)
      printf("[%s] Sorting done, total cputime = %14.6f\n", fdate(), tmp);

#endif

    /* counts the number of peaks in the sub-volume to know the number of groups */
    Npeaks = count_peaks(&Ngood);

    mynadd[0] = Ngood;
    mynadd[1] = Npeaks;
    MPI_Reduce(mynadd, nadd_all, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!ThisTask)
    {
      printf("[%s] Task 0 found %d peaks, %d in the well resolved region. Total number of peaks: %Ld\n",
             fdate(), Npeaks, Ngood, nadd_all[0]);
    }

    /* the number of peaks was supposed to be at most 1/10 of the number of particles
(we need space for Npeaks+2 groups) */
    if (Npeaks + 2 > subbox.PredNpeaks)
    {
      printf("ERROR on task %d: the number of peaks %d exceeds the predicted one (%d)\n",
             ThisTask, Npeaks, subbox.PredNpeaks);
      printf("      Please increase PredPeakFactor and restart\n");
      fflush(stdout);
      return 1;
    }

    if (turn)
    {
      /* Gives the minimal PredPeakFactor */
      BestPredPeakFactor = (double)nadd_all[1] / (double)MyGrids[0].ParticlesPerTask * 6.0;
    }

    /* sets arrays to zero */
    memset(group_ID, 0, subbox.Nalloc * sizeof(int));
    memset(linking_list, 0, subbox.Nalloc * sizeof(int));
    memset(groups, 0, subbox.PredNpeaks * sizeof(group_data));
    memset(wheretoplace_mycat, 0, subbox.PredNpeaks * sizeof(histories_data));
#ifdef PLC
    memset(plcgroups, 0, plc.Nmax * sizeof(plcgroup_data));
#endif

#ifndef CLASSIC_FRAGMENTATION
    if (!turn)
    {

      /* quick fragmentation is done in a shot */
      tmp = MPI_Wtime();

      /* quickly generates the group catalogue */
      if (quick_build_groups(Npeaks))
        return 1;

      tmp = MPI_Wtime() - tmp;
      cputime.group += tmp;
      if (!ThisTask)
        printf("[%s] Quick fragmentation done, cputime = %14.6f\n", fdate(), tmp);
    }
    else
#endif
    {
      /* full fragmentation is done segmenting the redshift interval */
      for (int mysegment = 0; mysegment < ScaleDep.nseg; mysegment++)
      {
        ScaleDep.myseg = mysegment;

#ifdef RECOMPUTE_DISPLACEMENTS
        /* The first segment uses only the first computed velocity
     The second segment has the two velocities already loaded */
```


----- FILE: src/fragment.h -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define NV 6
#define FILAMENT 1
#define SHIFT 0.5

#define ORDER_FOR_GROUPS 2
#define ORDER_FOR_CATALOG 3

extern int ngroups;

typedef struct
{
  int M,i;
  double q[3],v[3],D,z;
#ifdef TWO_LPT
  double D2,v2[3];
#ifdef THREE_LPT
  double D31,v31[3],D32,v32[3];
#endif
#endif
} pos_data;
pos_data obj, obj1, obj2;

typedef struct
{
  unsigned long long int name;
  double M,q[3],x[3],v[3];
  int n, pad;
}  catalog_data;

void condition_for_accretion(int, int, int, int, double, int, double *, double *);
void condition_for_merging(double, int, int, int *);
void set_obj(int, double, pos_data *);
void set_point(int, int, int, int, double, pos_data *);
void set_group(int, pos_data *);
double q2x(int, pos_data *, int);
double vel(int, pos_data *);
double distance(int, pos_data *, pos_data *);
void clean_list(int *);
double virial(int, double, int);
void merge_groups(int, int, double);
void update_history(int, int, double);
void accretion(int, int, int, int, int, double);
void update(pos_data *, pos_data *);
int write_catalog(int);
int write_histories(void);
int compute_mf(int);
#ifdef PLC
int write_PLC();
#endif

#ifdef MASS_MAPS
int write_mass_maps(double z_start, double z_end);
#endif
```


----- FILE: src/GenIC.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*

  The code contained in this file has been adapted from the 2LPTic code
  by Roman Scoccimarro, downloadable from:

  http://cosmo.nyu.edu/roman/2LPT/

  The 2LPTic code is a 2LPT extension of the N-GenIC code by V. Springel,
  downloadable from:

  http://www.gadgetcode.org/right.html#ICcode

  and distributed under GNU/GPL license.  The part of the 2LPTic code
  we use here is entirely contained in N-GenIC.

*/


#include "pinocchio.h"

#define GRID MyGrids[ThisGrid]

#define BOTTOM 0
#define TOP    1

typedef long long int map_index;
typedef int map_point[2];             // NOTE: this "limits" the plane to a side length of 2^31 points

unsigned int *OLDSEEDTABLE;


unsigned int *SEEDTABLE;
unsigned int *SEED_k0[2], *SEED_x0, *SEED_y0;
map_point     SEED_k0_BL[2];         // bottom-left corner of the seed regions for symmetries on k=0 plane
int           SEED_k0_nx[2];         // x-dimension of the seed regions for symmetries on k=0 plane
int           SEED_k0_ny[2];         // y-dimension of the seed regions for symmetries on k=0 plane


int  generate_seeds_plane(int, map_point, map_point, unsigned, int);
int  quadrant            (map_point, int);
void dump_seed_plane     (int, map_point, map_point, int[3], int, int, int);



int GenIC_large(int ThisGrid)
{
  int           Nmesh, Nsample, VectorLength;
  int           C[3], local_n[3], local_start[3];
  double        Box, fac;

  VectorLength   = GRID.total_local_size_fft;
  Nmesh          = GRID.GSglobal[_x_];
  Nsample        = GRID.GSglobal[_x_];
  Box            = GRID.BoxSize;

  fac = pow(1./Box,1.5);


  dprintf(VDBG, ThisTask, "==== %d in k-space:: (%td, %td) - (%td, %td) - (%td, %td)\n",
	  ThisTask,
	  GRID.GSstart_k[0], GRID.GSstart_k[0] + GRID.GSlocal_k[0],
	  GRID.GSstart_k[1], GRID.GSstart_k[1] + GRID.GSlocal_k[1],
	  GRID.GSstart_k[2], GRID.GSstart_k[2] + GRID.GSlocal_k[2]);


  // -----------------------------------------------
  // generate seed plane
  // bottom-left and top-right corners for the region of k-plane
  // that this Task must deal with


  // prepare coordinates
  if ( !params.use_transposed_fft )
    // this is the non-transposed order x, y, z
    for(int i = 0; i < 3; i++)
      C[i] = i;
  else
    {
      if ( internal.tasks_subdivision_dim > 1 )
  	C[0] = _y_, C[1] = _z_, C[2] = _x_;
      else
  	C[0] = _y_, C[1] = _x_, C[2] = _z_;
    }

  // set up starting points and lenghts
  for(int i = 0; i< 3; i++)
    {
      local_n    [i] = GRID.GSlocal_k[i];
      local_start[i] = GRID.GSstart_k[i];
    }

  // actually generate seeds
  {
    map_point     bottom_left, top_right;
    map_point     plane_bottom_left, plane_top_right;

    bottom_left[_x_] = GRID.GSstart_k[_x_];
    bottom_left[_y_] = GRID.GSstart_k[_y_];

    top_right[_x_]   = GRID.GSstart_k[_x_] + GRID.GSlocal_k[_x_];
    top_right[_y_]   = GRID.GSstart_k[_y_] + GRID.GSlocal_k[_y_];

    if( !internal.large_plane )
      {
	plane_bottom_left[_x_] = 0;
	plane_bottom_left[_y_] = 0;

	plane_top_right[_x_] = Nmesh;
	plane_top_right[_y_] = Nmesh;
      }
    else
      {
	plane_bottom_left[_x_] = bottom_left[_x_];
	plane_bottom_left[_y_] = bottom_left[_y_];
	plane_top_right[_x_]   = top_right[_x_];
	plane_top_right[_y_]   = top_right[_y_];
      }

    int ret = generate_seeds_plane(ThisGrid, plane_bottom_left, plane_top_right, Nmesh, (local_start[_z_] == 0));

    if( ret > 0 )
      return ret;

    if ( internal.dump_seedplane )
      {

	if ( params.use_transposed_fft && internal.tasks_subdivision_dim > 1 )
	  C[1] = _x_;
	dump_seed_plane(ThisGrid, bottom_left, top_right, C, local_n[_x_], GRID.GSstart_k[_z_], Nmesh);
      }
  }

  // -----------------------------------------------




  // the vector is initialized to zero
  memset( kdensity[ThisGrid], 0, VectorLength * sizeof(double) );

  int           Nmesh_2, Nmesh_odd;
  gsl_rng      *k0_generator;

  Nmesh_2   = Nmesh / 2;
  Nmesh_odd = (Nmesh % 2);


  // initialize the pseudo-random number chain

  gsl_rng_set(random_generator, params.RandomSeed);

  // initialize a second random generator, used on the plane k = 0

  k0_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  // --------------------------------------------------------------
  // this is the main loop on grid modes
  //--------------------------------------------------------------

  for( int i = 0; i < local_n[_x_]; i++ )
    {

      // ii is the x grid coordinate
      int ii = local_start[_x_] + i;
      if ( ii == Nmesh_2 )
	continue;

      double kvec[3];

      if ( ii < Nmesh_2 )
	kvec[0] = ii * 2 * PI / Box;
      else
	kvec[0] = -(Nmesh - ii) * 2 * PI / Box;

      double kmag2_i = kvec[0] * kvec[0];
      int iii = ii;


      for(int j = 0; j < local_n[_y_]; j++)
	{

	  // jj is the y grid coordinate
	  int jj = local_start[_y_] + j;
	  if(jj == Nmesh_2)
	    continue;
	  int jjj = jj;

	  if(jj < Nmesh_2)
	    kvec[1] = jj * 2 * PI / Box;
	  else
	    kvec[1] = -(Nmesh - jj) * 2 * PI / Box;

	  double kmag2_ij = kmag2_i + kvec[1] * kvec[1];

	  unsigned int seed;

	  // initialize the random generator with the seed found in the
	  // seedtable at the (i,j) coordinate in grid coordinates
	  if(internal.large_plane)
	    {
	      gsl_rng_set(random_generator, SEEDTABLE[j * local_n[_x_] + i]);
	      seed = SEEDTABLE[j * local_n[_x_] + i];
	    }
	  else
	    {
	      gsl_rng_set(random_generator, SEEDTABLE[jj * Nmesh + ii]);
	      seed = SEEDTABLE[jj * Nmesh + ii];
	    }

	  // since the chain of random numbers is subsequent, to
	  // reproduce all the time the same sequence, this
	  // process must generate all the random numbers on the
	  // k column that are below its starting k coordinate

	  double phase;
	  if(local_start[_z_] > 0 && local_start[_z_] < Nmesh_2)
	    for(int RR = 0; RR < local_start[_z_]; RR++)
	      {
	  	phase = gsl_rng_uniform(random_generator);
	  	do
	  	  phase = gsl_rng_uniform(random_generator);
	  	while(phase == 0);
	      }

	  // inner loop
	  for(int k = 0; (k < local_n[_z_]) && (k+local_start[_z_] < Nmesh_2); k++)
	    {
	      // kk is the z grid coordinate
	      int kk = local_start[_z_] + k;

	      // generate phase and amplitude
	      phase = gsl_rng_uniform(random_generator) * 2 * PI;
	      double ampl;
	      do
		ampl = gsl_rng_uniform(random_generator);
	      while(ampl == 0);

	      // blind points on the cube
	      if(ii == 0 && jj == 0 && kk == 0)
		continue;
	      if(kk == Nmesh_2)
		continue;

	      if(kk < Nmesh_2)
		kvec[2] = kk * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - kk) * 2 * PI / Box;

	      double kmag2_local = kmag2_ij + kvec[2] * kvec[2];
	      double kmag = sqrt(kmag2_local);

	      if(kmag * Box / (2 * PI) > NYQUIST * Nsample / 2)
		continue;

	      double p_of_k = PowerSpectrum(kmag);

	      double sign = 1.0;
	      int addr_j = j;

	      // special symmetries on the plane k = 0
	      if( kk == 0 )
		{

		  // some points are left empty
		  if( (ii == 0) && (jj == Nmesh_2 || jj == Nmesh_2 + Nmesh_odd) )
		    continue;
		  if( (ii == Nmesh_2) || (ii == Nmesh_2 + Nmesh_odd) )
		    continue;

		  // the half-plane for ii > Nmesh_2 takes the seeds
		  // from points in the other half-plane
		  if( ( ii > Nmesh_2 ) ||
		      ( ii == 0 && jj > Nmesh/2 ) )
		    {
		      // symmetries on j
		      jjj = Nmesh -jj;
		      if(jjj == Nmesh)
			jjj = 0;

		      if(Nmesh_odd && jj == Nmesh_2+1)
			{
			  jjj = Nmesh_2+1;
			  addr_j = Nmesh_2;
			}

		      // symmetries on i
		      if(ii > Nmesh_2)
			iii = Nmesh - ii;


		      //unsigned int seed;
		      if(internal.large_plane)
			// each task owns only a region of the seeds' plane
			{
			  if(jjj == 0)
			    // the line below would read
			    //   seed = SEED_y0[iii];
			    // however, since SEED_y0 starts from 0 we
			    // have to displace accordingly to i position
			    // reversed because of the requested symmetry
			    {
			      seed = SEED_y0[local_n[_x_] - i];
			    }
			  else if(iii == 0)
			    // the line below would read
			    //   seed = SEED_x0[jjj];
			    // however, since SEED_y0 starts from 0 we
			    // have to displace accordingly to j position
			    // reversed because of the requested symmetry
			    {
			      seed = SEED_x0[local_n[_y_] - j];
			    }
			  else
			    {
			      map_point V;
			      V[_x_] = ii;
			      V[_y_] = jj;
			      int Q = quadrant(V , Nmesh/2) / 2;

			      seed = SEED_k0[Q][(jjj - SEED_k0_BL[Q][_y_]) * SEED_k0_nx[Q] + (iii - SEED_k0_BL[Q][_x_]) ];
			    }
			}
		      else
			{
			  // each task owns the whole seeds plane
			  seed = SEEDTABLE[jjj * Nmesh + iii];
			}


		      sign = -1.0;

		      // re-initialize the spare random generator to the right seed
		      gsl_rng_set(k0_generator, seed);
		      phase = gsl_rng_uniform(k0_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(k0_generator);
		      while(ampl == 0);

		    }
		}

	      // Paired initial conditions
	      if (params.PairedIC)
		phase += PI;

	      // Fixed initial conditions
	      if (!params.FixedIC)
		p_of_k *= -log(ampl);

	      double delta = fac * sqrt(p_of_k);

	      int addr;

	      // calculate the storing address in local coordinates
	      if(!params.use_transposed_fft)
		addr    = 2*((i * local_n[_y_] + addr_j) * local_n[_z_] + k);
	      else
		{
		  // the addressing here is obviously awful and inefficient because
		  //   when the transposed fft are used the memory ordering is
		  //   row-major in z, y, x instead of x, y, z.
		  //
		  // however the nested loops have been designed for x, y, z order
		  //   that is "natural" because of the z=0 symmetries.
		  //
		  // due to this fact, for the sake of a clearer coding, we prefer
		  //   not to write a different version of the loops, taking into
		  //   account that this part of the code is not expèected to be
		  //   dominant

		  if(internal.tasks_subdivision_dim > 1)
		    addr    = 2*((addr_j * local_n[_z_] + k) * local_n[_x_] + i);
```


----- FILE: src/initialization.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

// #define BL_GRANDISSIMA

#ifdef PRECISE_TIMING // SI PUO` TOGLIERE?

#define SET_WTIME cputime.partial = MPI_Wtime();
#define ASSIGN_WTIME(INIT, ACC)       \
  do                                  \
  {                                   \
    double ttt = MPI_Wtime();         \
    cputime.ACC = ttt - cputime.INIT; \
  } while (0)
#define ACCUMULATE_WTIME(INIT, ACC)    \
  do                                   \
  {                                    \
    double ttt = MPI_Wtime();          \
    cputime.ACC += ttt - cputime.INIT; \
  } while (0)
#else
#define SET_WTIME
#define ASSIGN_WTIME(INIT, ACC)
#define ACCUMULATE_WTIME(INIT, ACC)
#endif

int initialize_fft(void);
int init_cosmology(void);
int generate_densities(void);
int set_plc(void);
#ifdef SCALE_DEPENDENT
int set_scaledep_GM(void);
#endif
unsigned int gcd(unsigned int, unsigned int);
int set_fft_decomposition(void);

int initialization()
{

  /* timing */
  cputime.init = MPI_Wtime();

  /* this is for gsl integration */
  workspace = gsl_integration_workspace_alloc(NWINT);
  /* this is the initialization of the random number generator */
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* reading of parameters from file and few other parameter initializations */
  if (set_parameters())
    return 1;

  /* call to cosmo.c: initialization of cosmological functions */
  if (initialize_cosmology())
    return 1;

  /* initialize pfft and fftw functions */
  if (initialize_fft())
    return 1;

  /* set the properties of grids and initialize FFTW quantities, including vectors */
  if (set_grids())
    return 1;

  /* computes the smoothing radii */
  if (set_smoothing())
    return 1;

  /* now it re-initializes the variance with a top-hat filter */
  WindowFunctionType = 2;
  if (initialize_MassVariance())
    return 1;

  /* initializes quantities needed for the on-the-fly reconstruction of PLC */
  SET_WTIME;
  if (set_plc())
    return 1;
  ASSIGN_WTIME(partial, set_plc);

  /* computes the number of sub-boxes for fragmentation */

  SET_WTIME;
  if (set_subboxes())
    return 1;
  ASSIGN_WTIME(partial, set_subboxes);

#ifdef SCALE_DEPENDENT
  /* computes the growth rates for displacements */
  if (set_scaledep_GM())
    return 1;
#endif

  /* checks that parameters and directives are coherent */
  if (check_parameters_and_directives())
    return 1;

  /* estimates the size of output file */
  if (estimate_file_size())
    return 1;

  /* this barrier is set to have correct stdout in case the run cannot start */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* allocations of memory for fmax and memory tests */
  SET_WTIME;
  if (allocate_main_memory())
    return 1;
  ASSIGN_WTIME(partial, memory_allocation);

  /* initialization of fft plans */
  SET_WTIME;
  if (compute_fft_plans())
    return 1;
  ASSIGN_WTIME(partial, fft_initialization);

  /* generation of initial density field */
  if (!params.ReadProductsFromDumps) /* no generation if products are read from dumps */
    if (generate_densities())
      return 1;

  cputime.init = MPI_Wtime() - cputime.init;

  if (!ThisTask)
  {
    dprintf(VMSG, ThisTask, "[%s] initialization done, initialization cpu time = %14.6f\n", fdate(), cputime.init);
    dprintf(VMSG, ThisTask, "\t\t set subboxes time = %14.6f\n"
                            "\t\t set plc time = %14.6f\n"
                            "\t\t memory allocation time = %14.6f\n"
                            "\t\t fft initialization time = %14.6f\n"
                            "\t\t density generation time = %14.6f\n",
            cputime.set_subboxes, cputime.set_plc, cputime.memory_allocation,
            cputime.fft_initialization, cputime.dens);
  }

  return 0;
}

int initialize_fft(void)
{

#ifdef USE_FFT_THREADS
  // if ( internal.nthreads_fft < 0 )
  internal.nthreads_fft = internal.nthreads_omp;
  if (internal.nthreads_fft > 1)
    dprintf(VMSG, 0, "Using %d threads for FFTs\n", internal.nthreads_fft);
#endif

  /* Initialize pfft */
  pfft_init();

  /* Inititalize fftw */
#ifdef USE_FFT_THREADS
  fftw_init_threads();
#endif
  fftw_mpi_init();

  if (set_fft_decomposition())
    return 1;

  if (!ThisTask)
    dprintf(VMSG, ThisTask, "cube subdivision [%d dim]: %d x %d x %d = %d processes\n",
            internal.tasks_subdivision_dim,
            internal.tasks_subdivision_3D[0],
            internal.tasks_subdivision_3D[1],
            internal.tasks_subdivision_3D[2],
            internal.tasks_subdivision_3D[0] *
                internal.tasks_subdivision_3D[1] *
                internal.tasks_subdivision_3D[2]);

  if (pfft_create_procmesh(internal.tasks_subdivision_dim, MPI_COMM_WORLD, internal.tasks_subdivision_3D, &FFT_Comm))
  {
    int all = 1;
    for (int iii = 0; iii < internal.tasks_subdivision_dim; iii++)
      all *= internal.tasks_subdivision_3D[iii];

    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error while creating communicator and mesh with %d processes\n", all);
    return 1;
  }

  return 0;
}

int set_parameters()
{
  int i;

  /* set default internal parameters */
  internal.verbose_level = VDIAG;
  internal.dump_seedplane = 0;
  internal.dump_kdensity = 0;
  internal.large_plane = 1;
  internal.mimic_original_seedtable = 0;
  // internal.dump_vectors                    = 0;
  internal.constrain_task_decomposition[0] = 0;
  internal.constrain_task_decomposition[1] = 0;
  internal.constrain_task_decomposition[2] = 0;
  internal.tasks_subdivision_3D[0] = 0;
  internal.tasks_subdivision_3D[1] = 0;
  internal.tasks_subdivision_3D[2] = 0;

  if (read_parameter_file())
    return 1;

  /* the smallest legitimate value of MinHaloMass is 1 */
  if (params.MinHaloMass <= 0)
    params.MinHaloMass = 1;

  if (params.BoxInH100)
  {
    params.BoxSize_h100 = params.BoxSize;
    params.BoxSize_htrue = params.BoxSize / params.Hubble100;
  }
  else
  {
    params.BoxSize_h100 = params.BoxSize * params.Hubble100;
    params.BoxSize_htrue = params.BoxSize;
  }
  params.InterPartDist = params.BoxSize_htrue / params.GridSize[0];

  params.ParticleMass = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0 * pow(params.InterPartDist, 3.);
  strcpy(params.DumpDir, "DumpProducts/");

  /* The Nyquist wavenumber is used in generic calls of scale-dependent growth rates */
  params.k_for_GM = PI / params.InterPartDist;

  if (!params.NumFiles)
    params.NumFiles = 1;

  /* The number of files must be a divisor of the number of tasks */
  if (NTasks % params.NumFiles != 0)
  {
    while (NTasks % params.NumFiles != 0)
      params.NumFiles--;

    if (!ThisTask)
      printf("Warning: NumFiles must be a divisor of NTasks, it has been fixed to %d\n",
             params.NumFiles);
  }

  /* inverse collapse times for the required outputs */
  for (i = 0; i < outputs.n; i++)
    outputs.F[i] = 1. + outputs.z[i];
  outputs.Flast = outputs.F[outputs.n - 1];

  if (!ThisTask)
  {
    dprintf(VMSG, 0, "Flag for this run: %s\n\n", params.RunFlag);
    dprintf(VMSG, 0, "PARAMETER VALUES from file %s:\n", params.ParameterFile);
    dprintf(VMSG, 0, "Omega0                      %f\n", params.Omega0);
    dprintf(VMSG, 0, "OmegaLambda                 %f\n", params.OmegaLambda);
    dprintf(VMSG, 0, "OmegaBaryon                 %f\n", params.OmegaBaryon);
    if (strcmp(params.TabulatedEoSfile, "no"))
    {
      dprintf(VMSG, 0, "Dark Energy EoS will be read from file %s\n", params.TabulatedEoSfile);
    }
    else
    {
      dprintf(VMSG, 0, "DE EoS parameters           %f %f\n", params.DEw0, params.DEwa);
    }

    dprintf(VMSG, 0, "Hubble100                   %f\n", params.Hubble100);
    dprintf(VMSG, 0, "Sigma8                      %f\n", params.Sigma8);
    dprintf(VMSG, 0, "PrimordialIndex             %f\n", params.PrimordialIndex);
    dprintf(VMSG, 0, "RandomSeed                  %d\n", params.RandomSeed);
    dprintf(VMSG, 0, "PairedIC                    %d\n", params.PairedIC);
    dprintf(VMSG, 0, "FixedIC                     %d\n", params.FixedIC);
    dprintf(VMSG, 0, "OutputList                  %s\n", params.OutputList);
    dprintf(VMSG, 0, "Number of outputs           %d\n", outputs.n);
    dprintf(VMSG, 0, "Output redshifts           ");
    for (i = 0; i < outputs.n; i++)
      dprintf(VMSG, 0, " %f ", outputs.z[i]);
    dprintf(VMSG, 0, "\n");
    dprintf(VMSG, 0, "GridSize                    %d %d %d\n", params.GridSize[0], params.GridSize[1], params.GridSize[2]);
    dprintf(VMSG, 0, "BoxSize (true Mpc)          %f\n", params.BoxSize_htrue);
    dprintf(VMSG, 0, "BoxSize (Mpc/h)             %f\n", params.BoxSize_h100);
    dprintf(VMSG, 0, "Particle Mass (true Msun)   %g\n", params.ParticleMass);
    dprintf(VMSG, 0, "Particle Mass (Msun/h)      %g\n", params.ParticleMass * params.Hubble100);
    dprintf(VMSG, 0, "Inter-part dist (true Mpc)  %f\n", params.InterPartDist);
    dprintf(VMSG, 0, "Inter-part dist (Mpc/h)     %f\n", params.InterPartDist * params.Hubble100);
    dprintf(VMSG, 0, "MinHaloMass (particles)     %d\n", params.MinHaloMass);
    dprintf(VMSG, 0, "MinHaloMass (Msun/h)        %g\n", params.MinHaloMass * params.ParticleMass * params.Hubble100);
    dprintf(VMSG, 0, "BoundaryLayerFactor         %f\n", params.BoundaryLayerFactor);
    dprintf(VMSG, 0, "MaxMem per task (Mb)        %d\n", params.MaxMem);
    dprintf(VMSG, 0, "MaxMem per particle (b)     %f\n", params.MaxMemPerParticle);
    dprintf(VMSG, 0, "PredPeakFactor              %f\n", params.PredPeakFactor);
    dprintf(VMSG, 0, "CatalogInAscii              %d\n", params.CatalogInAscii);
    dprintf(VMSG, 0, "NumFiles                    %d\n", params.NumFiles);
    dprintf(VMSG, 0, "DoNotWriteCatalogs          %d\n", params.DoNotWriteCatalogs);
    dprintf(VMSG, 0, "DoNotWriteHistories         %d\n", params.DoNotWriteHistories);
    dprintf(VMSG, 0, "WriteTimelessSnapshot       %d\n", params.WriteTimelessSnapshot);
    dprintf(VMSG, 0, "OutputInH100                %d\n", params.OutputInH100);
    dprintf(VMSG, 0, "DumpProducts                %d\n", params.DumpProducts);
    dprintf(VMSG, 0, "ReadProductsFromDumps       %d\n", params.ReadProductsFromDumps);
    dprintf(VMSG, 0, "ExitIfExtraParticles        %d\n", params.ExitIfExtraParticles);

    switch (params.AnalyticMassFunction)
    {
    case 0:
      dprintf(VMSG, 0, "Using Press & Schechter (1974) for the analytic mass function\n");
      break;
    case 1:
      dprintf(VMSG, 0, "Using Sheth & Tormen (2001) for the analytic mass function\n");
      break;
    case 2:
      dprintf(VMSG, 0, "Using Jenkins et al. (2001) for the analytic mass function\n");
      break;
    case 3:
      dprintf(VMSG, 0, "Using Warren et al. (2006) for the analytic mass function\n");
      break;
    case 4:
      dprintf(VMSG, 0, "Using Reed et al. (2007) for the analytic mass function\n");
      break;
    case 5:
      dprintf(VMSG, 0, "Using Crocce et al. (2010) for the analytic mass function\n");
      break;
    case 6:
      dprintf(VMSG, 0, "Using Tinker et al. (2008) for the analytic mass function\n");
      break;
    case 7:
      dprintf(VMSG, 0, "Using Courtin et al. (2010) for the analytic mass function\n");
      break;
    case 8:
      dprintf(VMSG, 0, "Using Angulo et al. (2012) for the analytic mass function\n");
      break;
    case 9:
      dprintf(VMSG, 0, "Using Watson et al. (2013) for the analytic mass function\n");
      break;
    case 10:
      dprintf(VMSG, 0, "Using Crocce et al. (2010) with forced universality for the analytic mass function\n");
      break;
    default:
      dprintf(VMSG, 0, "Unknown value for AnalyticMassFunction, Using Watson et al. (2013)\n");
      params.AnalyticMassFunction = 9;
      break;
    }
    dprintf(VMSG, 0, "\n");

    dprintf(VMSG, 0, "\n");
    dprintf(VMSG, 0, "GENIC parameters:\n");
    dprintf(VMSG, 0, "InputSpectrum_UnitLength_in_cm %f\n", params.InputSpectrum_UnitLength_in_cm);
    dprintf(VMSG, 0, "FileWithInputSpectrum          %s\n", params.FileWithInputSpectrum);
    dprintf(VMSG, 0, "WDM_PartMass_in_kev            %f\n", params.WDM_PartMass_in_kev);
#ifdef TABULATED_CT
    dprintf(VMSG, 0, "CTtableFile                    %s\n", params.CTtableFile);
#endif
#ifdef READ_PK_TABLE
    dprintf(VMSG, 0, "CAMBMatterFile                 %s\n", params.camb.MatterFile);
    dprintf(VMSG, 0, "CAMBRedshiftsFile              %s\n", params.camb.RedshiftsFile);
#endif
    dprintf(VMSG, 0, "\n");
  }

  /* Task 0 may have changed the value of this parameter */
  MPI_Bcast(&params.AnalyticMassFunction, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return 0;
}

int set_smoothing()
{
  int ismooth;
  double var_min, var_max, rmin;

  var_min = pow(1.686 / NSIGMA / GrowingMode(outputs.zlast, params.k_for_GM), 2.0);
  rmin = params.InterPartDist / 6.;
  var_max = MassVariance(rmin);
  Smoothing.Nsmooth = (log10(var_max) - log10(var_min)) / STEP_VAR + 2;

  if (Smoothing.Nsmooth <= 0)
  {
    if (!ThisTask)
      dprintf(VERR, 0, "I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
    Smoothing.Nsmooth = 1;
```


----- FILE: src/LPT.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"
#ifdef TWO_LPT

int compute_LPT_displacements(int compute_sources, double redshift)
{
  double time;

  /*
      second derivative components are sorted this way
       0 -> (1,1)
       1 -> (2,2)
       2 -> (3,3)
       3 -> (1,2)
       4 -> (1,3)
       5 -> (2,3)
  */

  if (compute_sources)
    {

      /**************************/
      /* Computation of sources */
      /**************************/

      ScaleDep.order=0; /* sources are computed as for LambdaCDM */
      ScaleDep.redshift=0.0;

  /* Computes the source for 2 LPT, 3LPT_1 and 3LPT_2 term */
  /* loop on all local particles */
#ifdef _OPENMP
#pragma omp parallel
      {

#pragma omp for nowait
#endif
	for (int index=0; index<MyGrids[0].total_local_size; index++)
	  {

	    /* HERE IT WILL BE NECESSARY TO SUM DIFFERENT SOURCES IF MULTIPLE GRIDS ARE USED */

	    /* source term for 2LPT */
	    source_2LPT[index] =
	      second_derivatives[0][0][index] * second_derivatives[0][1][index] +
	      second_derivatives[0][0][index] * second_derivatives[0][2][index] +
	      second_derivatives[0][1][index] * second_derivatives[0][2][index] -
	      second_derivatives[0][3][index] * second_derivatives[0][3][index] -
	      second_derivatives[0][4][index] * second_derivatives[0][4][index] -
	      second_derivatives[0][5][index] * second_derivatives[0][5][index];

#ifdef THREE_LPT
	    source_3LPT_1[index] = 3.0 *(second_derivatives[0][0][index]*
					 (second_derivatives[0][1][index]*second_derivatives[0][2][index]
					  -second_derivatives[0][5][index]*second_derivatives[0][5][index])
					 - second_derivatives[0][3][index]*
					 (second_derivatives[0][3][index]*second_derivatives[0][2][index]
					  -second_derivatives[0][4][index]*second_derivatives[0][5][index])
					 + second_derivatives[0][4][index]*
					 (second_derivatives[0][3][index]*second_derivatives[0][5][index]
					  -second_derivatives[0][4][index]*second_derivatives[0][1][index]));

	    source_3LPT_2[index] = 2.0 *     /* this factor is needed because nabla2phi is half the theoretical one */
	      (second_derivatives[0][0][index] + second_derivatives[0][1][index] + second_derivatives[0][2][index]) *
	      source_2LPT[index];
#endif
	  }
#ifdef _OPENMP
      }
#endif

      /* forward FFT for 2LPT source */
      write_in_rvector(0, source_2LPT);

      if (!ThisTask)
	printf("[%s] source for 2LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] source for 2LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft += time;

      write_from_cvector(0, kvector_2LPT);

#ifdef THREE_LPT

      /* second derivatives of second-order potential */
      for ( int ia = 1; ia <= 3; ia++ )
	for ( int ib = ia; ib <= 3; ib++ )
	  {
	    int ider = ( ia==ib? ia : ia+ib+1 );

	    if (!ThisTask)
	      printf("[%s] Computing 2nd derivative of 2LPT source: %d\n",fdate(),ia);

	    write_in_cvector(0, kvector_2LPT);

	    if (compute_derivative(0,ia,ib))
	      return 1;

#ifdef _OPENMP
#pragma omp parallel
	    {
#pragma omp for nowait
#endif
	      for (int index=0; index<MyGrids[0].total_local_size; index++)
		/* the first 2 factor is needed because nabla2phi is half the theoretical one */
		source_3LPT_2[index] -= 2.0 * (ider<=3? 1.0 : 2.0) *
		  rvector_fft[0][index] * second_derivatives[0][ider-1][index];
#ifdef _OPENMP
	    }
#endif
	  }

      /* forward FFT for 3LPT_1 source */
      write_in_rvector(0, source_3LPT_1);

      if (!ThisTask)
	printf("[%s] source for 3LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] sources for 3LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft+=time;

      write_from_cvector(0, kvector_3LPT_1);


      /* forward FFT for 3LPT_2 source */
      write_in_rvector(0, source_3LPT_2);

      if (!ThisTask)
	printf("[%s] source for 3LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] source for 3LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft+=time;

      write_from_cvector(0, kvector_3LPT_2);

#endif
    }

  /* displacements for the three terms */
  if (!ThisTask)
    printf("[%s] Computing displacements\n",fdate());

  ScaleDep.order=2; /* here we need the 2LPT  */
  ScaleDep.redshift=redshift;

  compute_first_derivatives(0., 0, 2, kvector_2LPT);

/*   for ( int ia = 1; ia <= 3; ia++ ) */
/*     { */
/*       if (!ThisTask) */
/* 	printf("[%s] Computing 1st derivative of 2LPT source: %d\n",fdate(),ia); */

/*       write_in_cvector(0, kvector_2LPT); */

/*       if (compute_derivative(0,ia,0)) */
/* 	return 1; */

/*       write_from_rvector(0, first_derivatives[0][ia-1]); */
/*     } */

/* #ifdef _OPENMP */
/* #pragma omp parallel */
/*   { */
/*     /\* assigns displacement to particles *\/ */

/* #pragma omp for nowait */
/* #endif */
/*     for (int index=0; index<MyGrids[0].total_local_size; index++) */
/*       for ( int ia = 0; ia < 3; ia++ ) */
/* 	products[index].Vel_2LPT[ia] = first_derivatives[0][ia][index]; */

/* #ifdef _OPENMP */
/*   } */
/* #endif */

#ifdef THREE_LPT

  if (!ThisTask)
    printf("[%s] Computing 3LPT_1 displacements\n",fdate());

  ScaleDep.order=3; /* here we need the 3LPT_1  */

  compute_first_derivatives(0., 0, 3, kvector_3LPT_1);

  if (!ThisTask)
    printf("[%s] Computing 3LPT_2 displacements\n",fdate());

  ScaleDep.order=4; /* here we need the 3LPT_2  */

  compute_first_derivatives(0., 0, 4, kvector_3LPT_2);

#endif


  /* bye! */
  return 0;
}

#endif
```


----- FILE: src/pinocchio.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "pinocchio.h"
#include "def_splines.h"

void abort_code(void);
void write_cputimes(void);

int main(int argc, char **argv, char **envp)
{

  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level); // Hybrid MPI and OPENMP parallel
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); // gives you your rank i.e. ID
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks); // size of your pool

#ifdef _OPENMP
  /* initialization of OpemMP */
#pragma omp parallel
  {
#pragma omp master
    internal.nthreads_omp = omp_get_num_threads();
  }
#endif

  /* timing of the code */
  cputime.total=MPI_Wtime();
  greetings();

  /* checks that the parameter file is given in the command line */
  if (argc<2)
    {
      if (!ThisTask)
	printf("Usage: pinocchio.x parameterfile\n");
      MPI_Finalize();
      return 0;
    }

  /* exit now if a snapshot is required and SNAPSHOT is not set */
#ifndef SNAPSHOT
  if (argc>=3 && atoi(argv[2])>1)
    {
      if (!ThisTask)
	printf("Sorry but you have to use the SNAPSHOT directive in compilation to write a snapshot\n");
      MPI_Finalize();
      return 0;
    }
#endif

  /* exit now if a snapshot is required and SNAPSHOT is not set */
#ifndef TABULATED_CT
  if (argc>=3 && atoi(argv[2])==1)
    {
      if (!ThisTask)
	printf("Sorry but you have to use the TABULATED_CT directive in compilation to write the collapse time table\n");
      MPI_Finalize();
      return 0;
    }
#endif


  /* initialization */
  memset(&params, 0, sizeof(param_data));
  strcpy(params.ParameterFile,argv[1]);
  if (initialization())
    abort_code();

  /*****************************************/
  /*********** Special behaviour ***********/
  /*****************************************/
  /* called as "pinocchio.x parameterfile 1" it computes and writes collapse time table, then exit */
  if (argc>=3 && atoi(argv[2])==1)
    {
#ifdef TABULATED_CT
      if (!ThisTask)
	{
	  printf("In this configuration pinocchio will only compute a table of collapse times\n");
	}

      /* CYCLE ON SMOOTHING RADII */
      for (int ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
	{
	  double cputmp=MPI_Wtime();

	  if (!ThisTask)
	    printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n",
		   fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
		   sqrt(Smoothing.Variance[ismooth]) );

	  if (initialize_collapse_times(ismooth,1))
	    return 1;

	  if (!ThisTask)
	    printf("[%s] Collapse times computed, cpu time =%f s\n",fdate(),cputmp);

	}

      if (!ThisTask)
	printf("Pinocchio done!\n");

      MPI_Finalize();

#endif
      return 0;
    }

  /* On request, it writes the density field */
  if (argc>=3 && atoi(argv[2])==2)
    {
#ifdef SNAPSHOT
      /* called as "pinocchio.x parameterfile 2" it writes the density field
	 in configuration space and exits */
      if (!ThisTask)
	printf("In this configuration pinocchio only writes the linear density field in a snapshot\n");

      for (int ThisGrid=0; ThisGrid<Ngrids; ThisGrid++)
	{
	  write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	  double time=reverse_transform(ThisGrid);
	  if (!ThisTask)
	    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);
	  write_from_rvector(ThisGrid, density[ThisGrid]);
	  if (write_density(ThisGrid))
	    abort_code();
	}
#else
      if (!ThisTask)
	printf("Please compile the code with SNAPSHOT directive to use this option\n");

#endif
      if (!ThisTask)
	{
	  write_cputimes();
	  printf("Pinocchio done!\n");
	}

      MPI_Finalize();
      return 0;
    }

  /* called as "pinocchio.x parameterfile 3" it computes displacements and writes a standard snapshot,
     then exit */
  if (argc>=3 && atoi(argv[2])==3)
    {
#ifdef SNAPSHOT

      if (!ThisTask)
	      {
	        printf("In this configuration pinocchio will only produce a GADGET snapshot\n");
	        printf("at the first redshift specified by the %s file (z=%f)\n", params.OutputList,outputs.z[0]);
        }

      /* compute displacements for all particles */
      double cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Computing displacements\n",fdate());

      if (compute_displacements(1,1,outputs.z[0]))
          abort_code();

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing displacements, cpu time = %f s\n",fdate(),cputmp);

      /* write the snapshot */
      cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Writing the snapshot\n",fdate());

      if (write_LPT_snapshot())
          abort_code();

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done snapshot, cpu time = %f s\n",fdate(),cputmp);
#endif

      if (!ThisTask)
        printf("Pinocchio done!\n");

      MPI_Finalize();

      return 0;
    }


  /******************************************/
  /*********** Standard behaviour ***********/
  /******************************************/


  if (params.ReadProductsFromDumps)
    {
      /* on request, read products from dump files and skip the first part */
      if (read_dumps())
	abort_code();
    }
  else
    {
      /* computation of collapse times and displacements */
      if (compute_fmax())
	abort_code();

      /* on request, dump products to files for skipping fmax */
      if (params.DumpProducts)
	if (dump_products())
	  abort_code();
    }

  /* fragmentation of the collapsed medium */
  if (fragment_driver())
    abort_code();

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* output detailed cpu times */
  cputime.total=MPI_Wtime()-cputime.total;
  if (!ThisTask)
    write_cputimes();

  /* done */
  if (!ThisTask)
    printf("Pinocchio done!\n");
  MPI_Finalize();

  return 0;
}


void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}


void write_cputimes()
{
  printf("Total:            %14.6f\n", cputime.total);
  printf("Initialization:   %14.6f (%5.2f%%)\n", cputime.init, 100.*cputime.init/cputime.total);
  printf("  Density in PS:  %14.6f (%5.2f%%)\n", cputime.dens, 100.*cputime.dens/cputime.total);
  printf("fmax:             %14.6f (%5.2f%%)\n", cputime.fmax, 100.*cputime.fmax /cputime.total);
#ifdef TWO_LPT
  printf("  LPT:            %14.6f (%5.2f%%)\n", cputime.lpt,  100.*cputime.lpt  /cputime.total);
#endif
  printf("  Derivatives:    %14.6f (%5.2f%%)\n", cputime.deriv,  100.*cputime.deriv  /cputime.total);
  printf("    Mem transfer: %14.6f (%5.2f%%)\n", cputime.mem_transf, 100.*cputime.mem_transf  /cputime.total);
  printf("    FFTs:         %14.6f (%5.2f%%)\n", cputime.fft,  100.*cputime.fft  /cputime.total);
  printf("  Collapse times: %14.6f (%5.2f%%)\n", cputime.coll, 100.*cputime.coll /cputime.total);
  printf("    inv.collapse: %14.6f (%5.2f%%)\n", cputime.invcoll, 100.*cputime.invcoll /cputime.total);
  printf("    ellipsoid:    %14.6f (%5.2f%%)\n", cputime.ell, 100.*cputime.ell /cputime.total);
  printf("  Velocities:     %14.6f (%5.2f%%)\n", cputime.vel,  100.*cputime.vel  /cputime.total);
  printf("Fragmentation:    %14.6f (%5.2f%%)\n", cputime.frag, 100.*cputime.frag /cputime.total);
  printf("  Redistribution: %14.6f (%5.2f%%)\n", cputime.distr,100.*cputime.distr/cputime.total);
  printf("  Sorting:        %14.6f (%5.2f%%)\n", cputime.sort, 100.*cputime.sort /cputime.total);
#ifdef PLC
  printf("  Groups total:   %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
  printf("  Groups PLC:     %14.6f (%5.2f%%)\n", cputime.plc,100.*cputime.plc/cputime.total);
#else
  printf("  Groups:         %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
#endif
  printf("Total I/O:        %14.6f (%5.2f%%)\n", cputime.io,   100.*cputime.io   /cputime.total);
}

```


----- FILE: src/pinocchio.h -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <fftw3-mpi.h>
#include <pfft.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* this library is used to vectorize the computation of collapse times */
/* #if !(defined(__aarch64__) || defined(__arm__)) */
/* #include <immintrin.h> */
/* #endif */

/* Defines */
#define NYQUIST 1.
#define PI 3.14159265358979323846
#define LBLENGTH 400
#define SBLENGTH 100
#define GBYTE 1073741824.0
#define MBYTE 1048576.0
#define alloc_verbose 0
#define MAXOUTPUTS 100
#define SPEEDOFLIGHT ((double)299792.458) /* km/s */
#define GRAVITY ((double)4.30200e-9)      /*  (M_sun^-1 (km/s)^2 Mpc)  */
#define NBINS 210                         /* number of time bins in cosmological quantities */
#define FRAGFIELDS 6

#define NSIGMA ((double)6.0)
#define STEP_VAR ((double)0.3) // 0.2)   /* this sets the spacing for smoothing radii */

#define NV 6
#define FILAMENT 1
#define SHIFT 0.5

#define ORDER_FOR_GROUPS 2
#define ORDER_FOR_CATALOG 3

#define ALIGN 32   /* for memory alignment */
#define UINTLEN 32 /* 8*sizeof(unsigned int) */

// #define ADD_RMAX_TO_SNAPSHOT

/* these templates define how to pass from coordinates to indices */
#define INDEX_TO_COORD(I, X, Y, Z, L) ({Z=(I)%L[_z_]; int _KK_=(I)/L[_z_]; Y=_KK_%L[_y_]; X=_KK_/L[_y_]; })
#define COORD_TO_INDEX(X, Y, Z, L) ((Z) + L[_z_] * ((Y) + L[_y_] * (X)))

/* coordinates */
#define _x_ 0
#define _y_ 1
#define _z_ 2

#define DECOMPOSITION_LIMIT_FACTOR_2D 1 /* smallest allowed side lenght of rectangular pencils in */
                                        /* 2D decomposition of FFT */

/* debug levels */
#define dprintf(LEVEL, TASK, ...)                                    \
  do                                                                 \
  {                                                                  \
    if (((LEVEL) <= internal.verbose_level) && (ThisTask == (TASK))) \
      fprintf(stdout, __VA_ARGS__);                                  \
  } while (1 == 0)
#define VDBG 4   // verbose level for debug
#define VDIAG 2  // verbose level for diagnostics
#define VMSG 1   // verbose level for flow messages
#define VXX 0    // essential messages
#define VERR VXX // non letal errors
#define VXERR -1 // letal errors

#define SWAP_INT(A, B) (A) ^= (B), (B) ^= (A), (A) ^= (B);

/*
  Global diagnostic wrapper for MPI_Reduce
  - Logs rank and file:line for each Reduce
  - Adds barriers before/after to catch divergence
  - Prints MPI error string and aborts on failure
  Enable by default; define DISABLE_MPI_REDUCE_DIAG to turn off.
*/
#define DISABLE_MPI_REDUCE_DIAG
#ifndef DISABLE_MPI_REDUCE_DIAG
static inline int pin_MPI_Reduce_diag(const void *sendbuf, void *recvbuf, int count,
                                      MPI_Datatype datatype, MPI_Op op, int root,
                                      MPI_Comm comm, const char *file, int line)
{
  int rank = -1, size = -1;
  PMPI_Comm_rank(comm, &rank);
  PMPI_Comm_size(comm, &size);
  fprintf(stderr, "[diag] rank %d/%d entering MPI_Reduce at %s:%d (count=%d, root=%d)\n",
          rank, size, file, line, count, root);
  fflush(stderr);
  PMPI_Barrier(comm);

  int rc = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  if (rc != MPI_SUCCESS)
  {
    char errstr[MPI_MAX_ERROR_STRING];
    int elen = 0;
    PMPI_Error_string(rc, errstr, &elen);
    fprintf(stderr, "[diag] MPI_Reduce FAILED at %s:%d on rank %d: %.*s\n",
            file, line, rank, elen, errstr);
    fflush(stderr);
    PMPI_Abort(comm, rc);
  }

  PMPI_Barrier(comm);
  if (rank == root)
  {
    fprintf(stderr, "[diag] MPI_Reduce OK at %s:%d (root=%d)\n", file, line, root);
    fflush(stderr);
  }
  return rc;
}

#define MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm) \
  pin_MPI_Reduce_diag((sendbuf), (recvbuf), (count), (datatype), (op), (root), (comm), __FILE__, __LINE__)
#endif /* DISABLE_MPI_REDUCE_DIAG */

/* checks of compiler flags */
#if defined(THREE_LPT) && !defined(TWO_LPT)
#define TWO_LPT
#endif

#if !defined(ELL_SNG) && !defined(ELL_CLASSIC)
#define ELL_CLASSIC
#endif

#if (defined(READ_PK_TABLE) || defined(MOD_GRAV_FR)) && !defined(SCALE_DEPENDENT)
#define SCALE_DEPENDENT
#endif

#if defined(READ_PK_TABLE) && defined(MOD_GRAV_FR)
#error REAK_PK_TABLE and MOD_GRAV_FR cannot be chosen together
#endif

#if defined(MOD_GRAV_FR) && !defined(FR0)
#error Please set a value to FR0 when you choose MOD_GRAV_FR
#endif

#if defined(MOD_GRAV_FR) && defined(FR0)
#warning "You have correctly compiled the code for the modified gravity scenario. However, please keep in mind that the modified gravity run (MOD_GRAV_FR) is still under development, and this mode should be used with extreme caution as it may not be fully stable. If you are unsure about its usage, please contact the developers for guidance."
#endif

/* vectorialization */
#define DVEC_SIZE 4

typedef double dvec __attribute__((vector_size(DVEC_SIZE * sizeof(double))));
typedef long int ivec __attribute__((vector_size(DVEC_SIZE * sizeof(long int))));

typedef union
{
  dvec V;
  double v[DVEC_SIZE];
} dvec_u;

typedef union
{
  ivec V;
  int v[DVEC_SIZE];
} ivec_u;

/* variables and type definitions */
extern int ThisTask, NTasks;
/* pfft-related variables */
/* extern int pfft_flags_c2r, pfft_flags_r2c; */
extern MPI_Comm FFT_Comm;

typedef struct
{
  int tasks_subdivision_dim;           /* 1, 2 or 3 to divide in slabs, pencils and volumes */
  int tasks_subdivision_3D[4];         /* ??? */
  int constrain_task_decomposition[3]; /* constraints on the number of subdivisions for each dimension */
  int verbose_level;                   /* for dprintf */
  int mimic_original_seedtable;        /* logical, set to 1 to reproduce exactly GenIC */
  // int dump_vectors;                     /* logical, dump vectors to files */
  int dump_seedplane; /* logical, dump seedplane to files */
  int dump_kdensity;  /* logical, dump Fourier-space density to files */
  int large_plane;    /* select the new generation of ICs */
  int nthreads_omp;   /* number of OMP threads */
  int nthreads_fft;   /* number of FFT threads */
} internal_data;
extern internal_data internal;

typedef unsigned int uint;
// typedef unsigned long long int UL;  // mi pare non ci sia

#ifdef DOUBLE_PRECISION_PRODUCTS
#define MPI_PRODFLOAT MPI_DOUBLE
typedef double PRODFLOAT;
#else
#define MPI_PRODFLOAT MPI_FLOAT
typedef float PRODFLOAT;
#endif

typedef struct // RIALLINEARE?
{
  int Rmax;
  PRODFLOAT Fmax, Vel[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1[3], Vel_3LPT_2[3];
#endif
#endif

#ifdef SNAPSHOT
  PRODFLOAT zacc;
  int group_ID;
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT Vel_prev[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT_prev[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1_prev[3], Vel_3LPT_2_prev[3];
#endif
#endif
#endif

} product_data __attribute__((aligned(ALIGN))); // VERIFICARE

extern char *main_memory, *wheretoplace_mycat;

extern product_data *products, *frag;

extern unsigned int *cubes_ordering;

extern unsigned int **seedtable;

extern double **kdensity;
extern double **density;
extern double ***first_derivatives;
extern double ***second_derivatives;

#ifdef TWO_LPT
extern double *kvector_2LPT;
extern double *source_2LPT;
#ifdef THREE_LPT
extern double *kvector_3LPT_1, *kvector_3LPT_2;
extern double *source_3LPT_1, *source_3LPT_2;
#endif
#endif

extern double Rsmooth;
typedef struct
{
  int Nsmooth;
  double *Radius, *Variance, *TrueVariance;
#ifdef SCALE_DEPENDENT
  double *Rad_GM, *k_GM_dens, *k_GM_displ, *k_GM_vel;
#endif
} smoothing_data;
extern smoothing_data Smoothing;

extern int Ngrids;
typedef struct
{
  unsigned int total_local_size, total_local_size_fft;
  unsigned int off, ParticlesPerTask;
  ptrdiff_t GSglobal[3];
  ptrdiff_t GSlocal[3];
  ptrdiff_t GSstart[3];
  ptrdiff_t GSlocal_k[3];
  ptrdiff_t GSstart_k[3];
  double lower_k_cutoff, upper_k_cutoff, norm, BoxSize, CellSize;
  pfft_plan forward_plan, reverse_plan;
  unsigned long long Ntotal;
} grid_data;
extern grid_data *MyGrids;

extern pfft_complex **cvector_fft;
extern double **rvector_fft;

#ifdef READ_PK_TABLE
typedef struct
{
  int Nkbins, NCAMB;
  char MatterFile[SBLENGTH], RedshiftsFile[LBLENGTH];
  double *Logk, *LogPkref, D2ref, *Scalef, *RefGM;
} camb_data;
#endif

typedef struct
{
  double Omega0, OmegaLambda, Hubble100, Sigma8, OmegaBaryon, DEw0, DEwa,
      PrimordialIndex, InterPartDist, BoxSize, BoxSize_htrue, BoxSize_h100, ParticleMass,
      StartingzForPLC, LastzForPLC, InputSpectrum_UnitLength_in_cm, WDM_PartMass_in_kev,
      BoundaryLayerFactor, Largest, MaxMemPerParticle, k_for_GM, PredPeakFactor, PLCAperture,
      PLCCenter[3], PLCAxis[3];
  char RunFlag[SBLENGTH], DumpDir[SBLENGTH], TabulatedEoSfile[LBLENGTH], ParameterFile[LBLENGTH],
      OutputList[LBLENGTH], FileWithInputSpectrum[LBLENGTH], CTtableFile[LBLENGTH];
  int GridSize[3], DumpProducts, ReadProductsFromDumps,
      CatalogInAscii, DoNotWriteCatalogs, DoNotWriteHistories, WriteTimelessSnapshot,
      OutputInH100, RandomSeed, MaxMem, NumFiles,
      BoxInH100, simpleLambda, AnalyticMassFunction, MinHaloMass, PLCProvideConeData, ExitIfExtraParticles,
      use_transposed_fft, FixedIC, PairedIC,
      NumMassPlanes,         /* number of mass planes for MASS_MAPS feature (0 disables) */
      MassMapNSIDE;          /* HEALPix NSIDE for MASS_MAPS (0 disables) */
  double MassMapMasterMaxGB; /* Max GB of memory rank 0 may use for one HEALPix plane (counts array) */
#ifdef READ_HUBBLE_TABLE
  char HubbleTableFile[LBLENGTH];
#endif
#ifdef READ_PK_TABLE
  camb_data camb;
#endif
} param_data;
extern param_data params;

typedef struct
{
  int n;
  double F[MAXOUTPUTS], z[MAXOUTPUTS], zlast, Flast;
} output_data;
extern output_data outputs;

typedef struct
{
  unsigned int Npart, Ngood, Nstored, PredNpeaks, maplength;
  unsigned int Nalloc, Nneeded;
  int nbox[3];
  int mybox[3];
  int Lgrid[3];
  int Lgwbl[3];
  int start[3];
  int stabl[3];
  int safe[3];
  int pbc[3];
  double SafetyBorder, overhead;
} subbox_data;
extern subbox_data subbox;

typedef struct
{
  double init, total, dens, fft, coll, invcoll, ell, vel, lpt, fmax, distr, sort, group, frag, io,
      deriv, mem_transf, partial, set_subboxes, set_plc, memory_allocation, fft_initialization
#ifdef PLC
      ,
      plc
#endif
      ;
} cputime_data;
extern cputime_data cputime;

extern int WindowFunctionType;

typedef struct
{
  int Mass;
  PRODFLOAT Pos[3], Vel[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1[3], Vel_3LPT_2[3];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT Vel_prev[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT_prev[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1_prev[3], Vel_3LPT_2_prev[3];
```


----- FILE: src/Pk_from_CAMB.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
   Refactored note:
   This file expects CAMB to output directly the CDM+baryon power spectrum P_cb(k,z)
   with k in h Mpc^-1 and P in (Mpc/h)^3. No transfer-function correction is applied
   and no unit conversion by Hubble100 is performed when reading the spectra.
*/

#ifdef SCALE_DEPENDENT_GROWTH

#include "pinocchio.h"

/* Optional diagnostics */
/* #define OUTPUT_GM */

static double *StoredLogK, *CAMBRedshifts, *CAMBScalefac, *RefGrowingMode,
    *StoredLogTotalPowerSpectrum;

static gsl_interp_accel *accelGrow = 0x0, *accel = 0x0;
static gsl_spline **splineGrowMatter, **splineGrowVel, **splineInvGrowMatter, **splineFomega, **splineFomega2, *spline = 0x0;

int ThisOutput, ThisRadius;

double IntegrandForSDMassVariance(double, void *);
double IntegrandForSDVelVariance(double, void *);
double PowerFromCAMB(double);

int read_power_table_from_CAMB()
{
    /* This routine reads P_cb(k) at various redshifts from a series of CAMB outputs
       and stores them in memory. It assumes inputs are already in h-units. */

    int i, j, dummy, ff;
    double kappa, Pk;
    char filename[BLENGTH], buffer[BLENGTH], *ugo;
    FILE *fd;

    if (!ThisTask)
    {
        /* count the number of CAMB files and, from the first, the number of data lines */
        params.camb.NCAMB = 0;
        sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, params.camb.NCAMB);
        while ((fd = fopen(filename, "r")) != 0x0)
        {
            if (!params.camb.NCAMB)
            {
                params.camb.Nkbins = 0;
                while (!feof(fd))
                {
                    ugo = fgets(buffer, BLENGTH, fd);
                    if (ugo && sscanf(buffer, "%lf", &kappa) == 1)
                        params.camb.Nkbins++;
                }
            }
            fclose(fd);
            params.camb.NCAMB++;
            sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, params.camb.NCAMB);
        }

        if (!params.camb.NCAMB)
        {
            printf("Error on Task 0: CAMB file %s not found\n", filename);
            return 1;
        }
        else if (!params.camb.Nkbins)
        {
            sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, 0);
            printf("Error on Task 0: problem in reading CAMB file %s\n", filename);
            return 1;
        }
        else if (params.camb.ReferenceOutput > params.camb.NCAMB - 1)
        {
            printf("Error on Task 0: ReferenceOutput is larger than NCAMB-1\n");
            return 1;
        }

        printf("Found %d CAMB matter power files with %d lines each (assuming P_cb in h-units)\n",
               params.camb.NCAMB, params.camb.Nkbins);
    }

    /* Broadcast NCAMB and Nkbins */
    MPI_Bcast(&params.camb.NCAMB, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.camb.Nkbins, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

    /* Allocate memory for spectra */
    StoredLogTotalPowerSpectrum = (double *)malloc(params.camb.NCAMB * params.camb.Nkbins * sizeof(double));
    StoredLogK = (double *)malloc(params.camb.Nkbins * sizeof(double));
    CAMBRedshifts = (double *)malloc(params.camb.NCAMB * sizeof(double));
    CAMBScalefac = (double *)malloc(params.camb.NCAMB * sizeof(double));
    RefGrowingMode = (double *)malloc(params.camb.NCAMB * sizeof(double));

    /* Task 0 reads all the CAMB files and stores the power spectra */
    if (!ThisTask)
    {
        for (i = 0; i < params.camb.NCAMB; i++)
        {
            sprintf(filename, "%s_%03d.dat", params.camb.MatterFile, i);
            fd = fopen(filename, "r");
            for (j = 0; j < params.camb.Nkbins; j++)
            {
                ff = fscanf(fd, "%lf %lf", &kappa, &Pk);
                (void)ff; /* suppress unused warning if not checked */
                /* Input already in h-units: k in h/Mpc and P in (Mpc/h)^3 */
                StoredLogTotalPowerSpectrum[i * params.camb.Nkbins + j] = log(Pk / pow(params.Hubble100, 3.));
                if (!i)
                    StoredLogK[j] = log(kappa * params.Hubble100);
            }
            fclose(fd);
        }

        if ((fd = fopen(params.camb.RedshiftsFile, "r")) == 0x0)
        {
            printf("Error: Redshift file %s not found\n", params.camb.RedshiftsFile);
            return 1;
        }
        for (i = 0; i < params.camb.NCAMB; i++)
            ff = fscanf(fd, "%d %lf", &dummy, CAMBRedshifts + i);
        fclose(fd);

        /* The last redshift MUST be z=0 */
        if (CAMBRedshifts[params.camb.NCAMB - 1] != 0.0)
        {
            printf("ERROR on Task 0: last CAMB redsbift must be 0.0\n");
            return 1;
        }
    }

    /* broadcast of loaded and computed quantities */
    MPI_Bcast(StoredLogK, params.camb.Nkbins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(StoredLogTotalPowerSpectrum, params.camb.NCAMB * params.camb.Nkbins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(CAMBRedshifts, params.camb.NCAMB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (i = 0; i < params.camb.NCAMB; i++)
    {
        CAMBScalefac[i] = 1. / (CAMBRedshifts[i] + 1.);
        RefGrowingMode[i] = 0.5 * (StoredLogTotalPowerSpectrum[i * params.camb.Nkbins + params.camb.ReferenceScale] -
                                   StoredLogTotalPowerSpectrum[params.camb.ReferenceOutput * params.camb.Nkbins + params.camb.ReferenceScale]);
    }
    params.camb.ReferenceRedshift = CAMBRedshifts[params.camb.ReferenceOutput];

    /* make P(k) at the reference redshift available to cosmo.c */
    params.camb.Logk = StoredLogK;
    params.camb.LogPkref = StoredLogTotalPowerSpectrum +
                           params.camb.ReferenceOutput * params.camb.Nkbins;
    params.camb.D2ref = exp(StoredLogTotalPowerSpectrum[(params.camb.NCAMB - 1) * params.camb.Nkbins + params.camb.ReferenceScale] -
                            StoredLogTotalPowerSpectrum[params.camb.ReferenceOutput * params.camb.Nkbins + params.camb.ReferenceScale]);
    params.camb.Scalef = CAMBScalefac;
    params.camb.RefGM = RefGrowingMode;

    return 0;
}

int initialize_ScaleDependentGrowth(void)
{

    /* Here it computes the scale-dependent growing modes at first and second order,
       for density and velocity, and the relative fomegas */

    int i, i1, i2;
    double *VarMatter, *VarVel, *ThisMatter, *ThisVel, *ThisfO, *ThisfO2;
    double result, error;
    gsl_function FuncMatter, FuncVel;
    double tolerance = 1.e-4;

    /* functions to integrate */
    FuncMatter.function = &IntegrandForSDMassVariance;
    FuncVel.function = &IntegrandForSDVelVariance;

    /* initialization of splines for growing modes */
    accelGrow = gsl_interp_accel_alloc();
    splineGrowMatter = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
    splineGrowVel = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
    splineInvGrowMatter = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
    splineFomega = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
    splineFomega2 = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
    for (i = 0; i < Smoothing.Nsmooth; i++)
    {
        splineGrowMatter[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
        splineGrowVel[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
        splineInvGrowMatter[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
        splineFomega[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
        splineFomega2[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
    }

    /* computes the needed growing modes and stores them */
    VarMatter = (double *)malloc(Smoothing.Nsmooth * sizeof(double));
    VarVel = (double *)malloc(Smoothing.Nsmooth * sizeof(double));
    ThisMatter = (double *)malloc(params.camb.NCAMB * sizeof(double));
    ThisVel = (double *)malloc(params.camb.NCAMB * sizeof(double));
    ThisfO = (double *)malloc(params.camb.NCAMB * sizeof(double));
    ThisfO2 = (double *)malloc(params.camb.NCAMB * sizeof(double));

    /* This is needed by PowerFromCAMB */
    accel = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, params.camb.Nkbins);

    /* values at the reference redshift */
    ThisOutput = params.camb.ReferenceOutput;
    for (ThisRadius = 0; ThisRadius < Smoothing.Nsmooth; ThisRadius++)
    {
        gsl_integration_qags(&FuncMatter, StoredLogK[0], StoredLogK[params.camb.Nkbins - 1], 0.0, tolerance, NWINT, workspace, &result, &error);
        VarMatter[ThisRadius] = result;
        gsl_integration_qags(&FuncVel, StoredLogK[0], StoredLogK[params.camb.Nkbins - 1], 0.0, tolerance, NWINT, workspace, &result, &error);
        VarVel[ThisRadius] = result;
    }

#ifdef OUTPUT_GM
    SDGM.flag = -1;
    FILE *fd;
    if (!ThisTask)
        fd = fopen("GrowingModes", "w");
#endif

    /* computation of growing modes for matter and velocity as a function of z,
       for all smoothing radii */
    for (ThisRadius = 0; ThisRadius < Smoothing.Nsmooth; ThisRadius++)
    {
        for (ThisOutput = 0; ThisOutput < params.camb.NCAMB; ThisOutput++)
        {
            gsl_integration_qags(&FuncMatter, StoredLogK[0], StoredLogK[params.camb.Nkbins - 1], 0.0, tolerance, NWINT, workspace, &result, &error);
            ThisMatter[ThisOutput] = sqrt(result / VarMatter[ThisRadius]);
            gsl_integration_qags(&FuncVel, StoredLogK[0], StoredLogK[params.camb.Nkbins - 1], 0.0, tolerance, NWINT, workspace, &result, &error);
            ThisVel[ThisOutput] = sqrt(result / VarVel[ThisRadius]);
        }
        /* initialization of splines */
        gsl_spline_init(splineGrowMatter[ThisRadius], CAMBScalefac, ThisMatter, params.camb.NCAMB);
        gsl_spline_init(splineGrowVel[ThisRadius], CAMBScalefac, ThisVel, params.camb.NCAMB);
        gsl_spline_init(splineInvGrowMatter[ThisRadius], ThisMatter, CAMBScalefac, params.camb.NCAMB);

        /* fomega */
        for (ThisOutput = 0; ThisOutput < params.camb.NCAMB; ThisOutput++)
        {

            if (ThisOutput < params.camb.NCAMB - 1)
            {
                i1 = (ThisOutput > 0 ? ThisOutput - 1 : 0);
                i2 = ThisOutput + 1; //(ThisOutput<params.camb.NCAMB-1 ? ThisOutput+1 : params.camb.NCAMB-1);

                ThisfO[ThisOutput] = (ThisVel[i2] - ThisVel[i1]) / (CAMBScalefac[i2] - CAMBScalefac[i1]) * CAMBScalefac[ThisOutput] / ThisVel[ThisOutput];
                ThisfO2[ThisOutput] = ((3. / 7. * pow(ThisVel[i2], 2.0) * pow(Omega(CAMBRedshifts[i2]), -0.007)) -
                                       (3. / 7. * pow(ThisVel[i1], 2.0) * pow(Omega(CAMBRedshifts[i1]), -0.007))) /
                                      (CAMBScalefac[i2] - CAMBScalefac[i1]) * CAMBScalefac[ThisOutput] /
                                      (3. / 7. * pow(ThisVel[ThisOutput], 2.0) * pow(Omega(CAMBRedshifts[ThisOutput]), -0.007));
            }
            else
            {
                ThisfO[ThisOutput] = -ThisfO[ThisOutput - 2] * (CAMBScalefac[ThisOutput] - CAMBScalefac[ThisOutput - 1]) / (CAMBScalefac[ThisOutput - 1] - CAMBScalefac[ThisOutput - 2]) + ThisfO[ThisOutput - 1] * (CAMBScalefac[ThisOutput] - CAMBScalefac[ThisOutput - 2]) / (CAMBScalefac[ThisOutput - 1] - CAMBScalefac[ThisOutput - 2]);
                ThisfO2[ThisOutput] = -ThisfO2[ThisOutput - 2] * (CAMBScalefac[ThisOutput] - CAMBScalefac[ThisOutput - 1]) / (CAMBScalefac[ThisOutput - 1] - CAMBScalefac[ThisOutput - 2]) + ThisfO2[ThisOutput - 1] * (CAMBScalefac[ThisOutput] - CAMBScalefac[ThisOutput - 2]) / (CAMBScalefac[ThisOutput - 1] - CAMBScalefac[ThisOutput - 2]);
            }

#ifdef OUTPUT_GM
            if (!ThisTask)
                fprintf(fd, " %d %d   %g %g   %g %g %g   %g %g %g   %g %g %f  %g %g\n",
                        ThisRadius, ThisOutput,
                        CAMBRedshifts[ThisOutput], Smoothing.Radius[ThisRadius],
                        ThisMatter[ThisOutput], ThisVel[ThisOutput], GrowingMode(CAMBRedshifts[ThisOutput]),
                        3. / 7. * pow(ThisMatter[ThisOutput], 2.0) * pow(Omega(CAMBRedshifts[ThisOutput]), -0.007),
                        3. / 7. * pow(ThisVel[ThisOutput], 2.0) * pow(Omega(CAMBRedshifts[ThisOutput]), -0.007),
                        GrowingMode_2LPT(CAMBRedshifts[ThisOutput]),
                        ThisfO[ThisOutput], fomega(CAMBRedshifts[ThisOutput]), Omega(CAMBRedshifts[ThisOutput]),
                        ThisfO2[ThisOutput], fomega_2LPT(CAMBRedshifts[ThisOutput]));
#endif
        }
        /* initialization of splines */
        gsl_spline_init(splineFomega[ThisRadius], CAMBScalefac, ThisfO, params.camb.NCAMB);
        gsl_spline_init(splineFomega2[ThisRadius], CAMBScalefac, ThisfO2, params.camb.NCAMB);
    }

#ifdef OUTPUT_GM
    if (!ThisTask)
        fclose(fd);
#endif

    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);

    free(ThisfO2);
    free(ThisfO);
    free(ThisVel);
    free(ThisMatter);
    free(VarVel);
    free(VarMatter);

    return 0;
}

double IntegrandForSDMassVariance(double logk, void *param)
{
    double k, R;

    k = exp(logk);
    R = (ThisRadius < Smoothing.Nsmooth - 1 ? Smoothing.Radius[ThisRadius] : params.InterPartDist / 6.);
    return PowerFromCAMB(k) * exp(-k * k * R * R) * k * k * k / (2. * PI * PI);
}

double IntegrandForSDVelVariance(double logk, void *radius)
{
    double k, R;

    k = exp(logk);
    R = (ThisRadius < Smoothing.Nsmooth - 1 ? Smoothing.Radius[ThisRadius] : params.InterPartDist / 6.);
    return PowerFromCAMB(k) * exp(-k * k * R * R) * k / (2. * PI * PI);
}

double PowerFromCAMB(double k)
{
    /* This function gives the interpolated power spectrum from the CAMB outputs.
       The power spectrum is interpolated from */

    static int LastOutput = -99;

    /* initialize gsl spline the first time the function is called, or when output changes */
    if (spline == 0x0 || ThisOutput != LastOutput)
    {
        gsl_spline_init(spline, StoredLogK,
                        StoredLogTotalPowerSpectrum + ThisOutput * params.camb.Nkbins,
                        params.camb.Nkbins);
    }

    LastOutput = ThisOutput;
    double result = exp(my_spline_eval(spline, log(k), accel));
    return result;
}

double MatterGrowingMode(double z)
{

    int mysmooth, interpolate;
    double w;

    /* computes the smoothing radius to be used (including the case of interpolation) */
    if (SDGM.ismooth >= 0 && SDGM.ismooth < Smoothing.Nsmooth)
    {
        mysmooth = SDGM.ismooth;
        interpolate = 0;
    }
    else
    {
        if (SDGM.radius > Smoothing.Radius[0])
        {
            mysmooth = 0;
            interpolate = 0;
        }
        else if (SDGM.radius < Smoothing.Radius[Smoothing.Nsmooth - 1])
        {
            mysmooth = Smoothing.Nsmooth - 1;
            interpolate = 0;
        }
        else
        {
            for (mysmooth = 1; mysmooth < Smoothing.Nsmooth && SDGM.radius <= Smoothing.Radius[mysmooth]; mysmooth++)
                ;
            interpolate = 1;
        }
    }

    if (interpolate)
    {
        w = (SDGM.radius - Smoothing.Radius[mysmooth]) / (Smoothing.Radius[mysmooth - 1] - Smoothing.Radius[mysmooth]);
        return ((1. - w) * my_spline_eval(splineGrowMatter[mysmooth], 1. / (1. + z), accelGrow) + w * my_spline_eval(splineGrowMatter[mysmooth - 1], 1. / (1. + z), accelGrow));
    }
    else
        return my_spline_eval(splineGrowMatter[mysmooth], 1. / (1. + z), accelGrow);
}

double VelGrowingMode(double z)
{

    int mysmooth, interpolate;
    double w;

    /* computes the smoothing radius to be used (including the case of interpolation) */
    if (SDGM.ismooth >= 0 && SDGM.ismooth < Smoothing.Nsmooth)
    {
```


----- FILE: src/ReadParamfile.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
   The similarity of this code with the read_parameter_file() function in the Gadget code,
   by V. Springel (http://www.gadgetcode.org), is not a coincidence.
   Many thanks to Volker for the inspiration.
*/

#include "pinocchio.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define LOGICAL 4
#define INT3 5
#define DOUBLE3 6
#define INT_SKIP_DEF 98 /* this is previously set to the default value */
#define INT_SKIP 99     /* this is set to 0 if not present in the parameter file */
#define MAXTAGS 100

int read_parameter_file()
{
  FILE *fd;
  char buf[SBLENGTH], buf1[SBLENGTH], buf2[SBLENGTH], buf3[SBLENGTH], buf4[SBLENGTH];
  int i, j, nt, number_of_fields;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][SBLENGTH];
  double z;

  if (!ThisTask)
  {

    nt = 0;
    int idx_NumMassPlanes = -1; /* track presence to inform user if missing */
    int idx_MassMapNSIDE = -1;
    int idx_MassMapMasterMaxGB = -1;

    /* list of requested parameters */
    strcpy(tag[nt], "RunFlag");
    addr[nt] = params.RunFlag;
    id[nt++] = STRING;

    strcpy(tag[nt], "OutputList");
    addr[nt] = params.OutputList;
    id[nt++] = STRING;

    strcpy(tag[nt], "BoxSize");
    addr[nt] = &params.BoxSize;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "BoxInH100");
    addr[nt] = &params.BoxInH100;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "GridSize");
    addr[nt] = &(params.GridSize[0]);
    id[nt++] = INT;

    strcpy(tag[nt], "RandomSeed");
    addr[nt] = &params.RandomSeed;
    id[nt++] = INT;

    strcpy(tag[nt], "Omega0");
    addr[nt] = &params.Omega0;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "OmegaLambda");
    addr[nt] = &params.OmegaLambda;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "OmegaBaryon");
    addr[nt] = &params.OmegaBaryon;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "Hubble100");
    addr[nt] = &params.Hubble100;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "Sigma8");
    addr[nt] = &params.Sigma8;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PrimordialIndex");
    addr[nt] = &params.PrimordialIndex;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "DEw0");
    addr[nt] = &params.DEw0;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "DEwa");
    addr[nt] = &params.DEwa;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "TabulatedEoSfile");
    addr[nt] = params.TabulatedEoSfile;
    id[nt++] = STRING;

    strcpy(tag[nt], "FileWithInputSpectrum");
    addr[nt] = params.FileWithInputSpectrum;
    id[nt++] = STRING;

    strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
    addr[nt] = &(params.InputSpectrum_UnitLength_in_cm);
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "WDM_PartMass_in_kev");
    addr[nt] = &(params.WDM_PartMass_in_kev);
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "BoundaryLayerFactor");
    addr[nt] = &params.BoundaryLayerFactor;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "MaxMem");
    addr[nt] = &params.MaxMem;
    id[nt++] = INT;

    strcpy(tag[nt], "MaxMemPerParticle");
    addr[nt] = &(params.MaxMemPerParticle);
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PredPeakFactor");
    addr[nt] = &(params.PredPeakFactor);
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "CatalogInAscii");
    addr[nt] = &params.CatalogInAscii;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "OutputInH100");
    addr[nt] = &params.OutputInH100;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "NumFiles");
    addr[nt] = &params.NumFiles;
    id[nt++] = INT_SKIP;

    strcpy(tag[nt], "MinHaloMass");
    addr[nt] = &params.MinHaloMass;
    id[nt++] = INT;

    strcpy(tag[nt], "AnalyticMassFunction");
    addr[nt] = &params.AnalyticMassFunction;
    id[nt++] = INT;

    strcpy(tag[nt], "WriteTimelessSnapshot");
    addr[nt] = &params.WriteTimelessSnapshot;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "DoNotWriteCatalogs");
    addr[nt] = &params.DoNotWriteCatalogs;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "DoNotWriteHistories");
    addr[nt] = &params.DoNotWriteHistories;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "StartingzForPLC");
    addr[nt] = &params.StartingzForPLC;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "LastzForPLC");
    addr[nt] = &params.LastzForPLC;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PLCAperture");
    addr[nt] = &params.PLCAperture;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PLCProvideConeData");
    addr[nt] = &params.PLCProvideConeData;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "PLCCenter");
    addr[nt] = &(params.PLCCenter);
    id[nt++] = DOUBLE3;

    strcpy(tag[nt], "PLCAxis");
    addr[nt] = &(params.PLCAxis);
    id[nt++] = DOUBLE3;

    strcpy(tag[nt], "FixedIC");
    addr[nt] = &(params.FixedIC);
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "PairedIC");
    addr[nt] = &(params.PairedIC);
    id[nt++] = LOGICAL;

#ifdef TABULATED_CT
    strcpy(tag[nt], "CTtableFile");
    addr[nt] = params.CTtableFile;
    id[nt++] = STRING;
#endif

#ifdef READ_PK_TABLE
    strcpy(tag[nt], "CAMBMatterFile");
    addr[nt] = params.camb.MatterFile;
    id[nt++] = STRING;

    strcpy(tag[nt], "CAMBRedshiftsFile");
    addr[nt] = params.camb.RedshiftsFile;
    id[nt++] = STRING;
#endif

#ifdef READ_HUBBLE_TABLE
    strcpy(tag[nt], "HubbleTableFile");
    addr[nt] = params.HubbleTableFile;
    id[nt++] = STRING;
#endif

    strcpy(tag[nt], "DumpProducts");
    addr[nt] = &params.DumpProducts;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "ReadProductsFromDumps");
    addr[nt] = &params.ReadProductsFromDumps;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "ExitIfExtraParticles");
    addr[nt] = &params.ExitIfExtraParticles;
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "UseTransposedFFT");
    addr[nt] = &(params.use_transposed_fft);
    id[nt++] = LOGICAL;

    /* Optional: number of mass planes for MASS_MAPS feature */
    strcpy(tag[nt], "NumMassPlanes");
    addr[nt] = &params.NumMassPlanes;
    id[nt++] = INT_SKIP_DEF;

    /* HEALPix mass map parameters */
    strcpy(tag[nt], "MassMapNSIDE");
    addr[nt] = &params.MassMapNSIDE;
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "MassMapMasterMaxGB");
    addr[nt] = &params.MassMapMasterMaxGB;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "MimicOldSeed");
    addr[nt] = &(internal.mimic_original_seedtable);
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "DumpSeedPlane");
    addr[nt] = &(internal.dump_seedplane);
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "DumpKDensity");
    addr[nt] = &(internal.dump_kdensity);
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "VerboseLevel");
    addr[nt] = &(internal.verbose_level);
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "LargePlane");
    addr[nt] = &(internal.large_plane);
    id[nt++] = LOGICAL;

    strcpy(tag[nt], "Constrain_dim0");
    addr[nt] = &(internal.constrain_task_decomposition[0]);
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "Constrain_dim1");
    addr[nt] = &(internal.constrain_task_decomposition[1]);
    id[nt++] = INT_SKIP_DEF;

    strcpy(tag[nt], "Constrain_dim2");
    addr[nt] = &(internal.constrain_task_decomposition[2]);
    id[nt++] = INT_SKIP_DEF;

    for (j = 0; j < nt; j++) /* All logical tags are FALSE by default */
      if (id[j] == LOGICAL)
        *((int *)addr[j]) = 0;

    printf("Reading parameters from file %s\n", params.ParameterFile);
    fflush(stdout);

    if ((fd = fopen(params.ParameterFile, "r")))
    {
      while (!feof(fd))
      {
        *buf = 0;
        (void)fgets(buf, SBLENGTH, fd);
        number_of_fields = sscanf(buf, "%s %s %s %s", buf1, buf2, buf3, buf4);
        if (number_of_fields < 1)
          continue;

        if (buf1[0] == '#' || buf1[0] == '%')
          continue;

        /* searches for a match with the tag list */
        for (i = 0, j = -1; i < nt; i++)
          if (strcmp(buf1, tag[i]) == 0)
          {
            j = i;
            tag[i][0] = 0;
            break;
          }

        if (j >= 0)
        {
          switch (id[j])
          {
          case DOUBLE:
            if (number_of_fields < 2)
              j = -10;
            else
              *((double *)addr[j]) = atof(buf2);
            break;

          case STRING:
            if (number_of_fields < 2)
              j = -10;
            else
              strcpy(addr[j], buf2);
            break;

          case INT:
            if (number_of_fields < 2)
              j = -10;
            else
              *((int *)addr[j]) = atoi(buf2);
            break;

          case INT_SKIP:
            if (number_of_fields < 2)
              j = -10;
            else
              *((int *)addr[j]) = atoi(buf2);
            break;

          case INT_SKIP_DEF:
            if (number_of_fields < 2)
              j = -10;
            else
              *((int *)addr[j]) = atoi(buf2);
            break;

          case INT3:
            if (number_of_fields < 4)
              j = -10;
            else
            {
              *((int *)addr[j]) = atoi(buf2);
              *((int *)addr[j] + 1) = atoi(buf3);
              *((int *)addr[j] + 2) = atoi(buf4);
            }
            break;

          case DOUBLE3:
            if (number_of_fields < 4)
              j = -10;
            else
            {
              *((double *)addr[j]) = atof(buf2);
              *((double *)addr[j] + 1) = atof(buf3);
              *((double *)addr[j] + 2) = atof(buf4);
            }
            break;

          case LOGICAL:
            *((int *)addr[j]) = 1;
            break;

          default:
            break;
          }

          if (j < 0)
```


----- FILE: src/ReadWhiteNoise.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

/* This routine reads in a white noise field as generated by GRAFIC2
   and adds the correct power */

// TUTTO DA RIVERIFICARE!!!

#ifdef WHITENOISE

int read_white_noise()
{
  int dummy,np1,np2,np3,seed,Ninplane,gz,lx,ly,lz,Task,i;
  int local_x,local_y,local_z,ixx,iyy,izz,index,nxhalf,nyhalf,nzhalf;
  size_t planesize;
  float *plane;
  double kx,ky,kz,kxnorm,kynorm,kznorm,k_squared,k_module,k_physical,time,norm,delta;
  double ave=0.0, var=0.0, gave, gvar;
  char filename[LBLENGTH];
  FILE *file;
  MPI_Status status;

  Ninplane=MyGrids[0].GSglobal_x * MyGrids[0].GSglobal_y;
  planesize=Ninplane * sizeof(float);
  plane=(float*)malloc(planesize);

  if (!ThisTask)
    {
      strcpy(filename,"WhiteNoise");
      file=fopen(filename,"r");
      if (file==0x0)
	{
	  printf("ERROR on Task 0: cannot find file WhiteNoise required by read_white_noise\n");
	  return 1;
	}

      fread(&dummy, 1, sizeof(int), file);
      fread(&np1, 1, sizeof(int), file);
      fread(&np2, 1, sizeof(int), file);
      fread(&np3, 1, sizeof(int), file);
      fread(&seed, 1, sizeof(int), file);
      fread(&dummy, 1, sizeof(int), file);

      if (np1 != MyGrids[0].GSglobal_x || np2 != MyGrids[0].GSglobal_y || np3 != MyGrids[0].GSglobal_z)
	{
	  printf("ERROR on Task 0: WhiteNoise file has wrong dimensions: %d %d %d\n",np1,np2,np3);
	  return 1;
	}
    }

  for (gz=0; gz<MyGrids[0].GSglobal_z; gz++)
    {
      /* Task 0 reads the plane */
      if (!ThisTask)
	{
	  fread(&dummy,1,sizeof(int),file);
	  if (dummy != planesize)
	    {
	      printf("ERROR on Task 0: mis-matched record size for file WhiteNoise, I expected %d but obtained %d\n",
		     (int)planesize, dummy);
	      return 1;
	    }

	  fread(plane, Ninplane, sizeof(float), file);
	  fread(&dummy,1,sizeof(int),file);

	  printf("Task 0 has read plane %d\n",gz);
	}

      /* each task writes into an int variable either 0 or its task number if the plane belongs to it */
      if (gz>=MyGrids[0].GSstart_z && gz<MyGrids[0].GSstart_z + MyGrids[0].GSlocal_z)
	dummy=ThisTask;
      else
	dummy=0;

      /* This way the variable Task (for Task 0) contains the task that must receive the plane */
      MPI_Reduce(&dummy, &Task, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      if (!ThisTask)  /* Task 0 either stores the data or sends them to Task */
	{
	  printf("The plane will be stored by Task %d\n",Task);

	  if (Task)
	    MPI_Send(plane, planesize, MPI_BYTE, Task, 0, MPI_COMM_WORLD);
	  else
	    {
	      for (i=0; i<Ninplane; i++)
		{
		  ave+=plane[i];
		  var+=plane[i]*plane[i];
		}
	      lz=gz-MyGrids[0].GSstart_z;
	      for (ly=0; ly<MyGrids[0].GSlocal_y; ly++)
		for (lx=0; lx<MyGrids[0].GSlocal_x; lx++)
		  *(rvector_fft[0] + lx + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (ly + lz * MyGrids[0].GSlocal_y)) =
		    *(plane + lx + MyGrids[0].GSlocal_x * ly);

	    }
	}
      else     /* other tasks receive the data and store them if dummy>0 */
	{
	  if (dummy==ThisTask)
	    {
	      MPI_Recv(plane, planesize, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
	      for (i=0; i<Ninplane; i++)
		{
		  ave+=plane[i];
		  var+=plane[i]*plane[i];
		}
	      lz=gz-MyGrids[0].GSstart_z;
	      for (ly=0; ly<MyGrids[0].GSlocal_y; ly++)
		for (lx=0; lx<MyGrids[0].GSlocal_x; lx++)
		  *(rvector_fft[0] + lx + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (ly + lz * MyGrids[0].GSlocal_y)) =
		    *(plane + lx + MyGrids[0].GSlocal_x * ly);
	    }
	}
    }

  MPI_Reduce(&ave, &gave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&var, &gvar, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (!ThisTask)
    {
      gave /= (double)MyGrids[0].GSglobal_x*(double)MyGrids[0].GSglobal_y*(double)MyGrids[0].GSglobal_z;
      gvar = sqrt(gvar/((double)MyGrids[0].GSglobal_x*(double)MyGrids[0].GSglobal_y*(double)MyGrids[0].GSglobal_z));
      printf("Average and variance of white noise: %f, %f\n",gave,gvar);
    }

  free(plane);


  if (!ThisTask)
    printf("[%s] read_white_noise: starting fft\n",fdate());

  time = forward_transform(0);

  if (!ThisTask)
    printf("[%s] read_white_noise: done fft, cputime = %f\n",fdate(),time);

  cputime.fft+=time;


  /* loop on the cvector to fix the power spectrum */
  /* k vectors */
  kxnorm   = 2.*PI/(double)MyGrids[0].GSglobal_x;
  kynorm   = 2.*PI/(double)MyGrids[0].GSglobal_y;
  kznorm   = 2.*PI/(double)MyGrids[0].GSglobal_z;

  /* Nyquist frequencies */
  nxhalf = MyGrids[0].GSglobal_x/2;
  nyhalf = MyGrids[0].GSglobal_y/2;
  nzhalf = MyGrids[0].GSglobal_z/2;

  norm = 1.0/pow(params.InterPartDist, 3.0);

  for (local_z = 0; local_z < MyGrids[0].GSlocal_k_z; local_z++)
    {
      izz = local_z + MyGrids[0].GSstart_k_z;
      if (local_z > nzhalf)
	izz -= MyGrids[0].GSglobal_z;
      kz  = kznorm*izz;

      for (local_y = 0; local_y < MyGrids[0].GSlocal_k_y; local_y++)
	{
	  iyy = local_y + MyGrids[0].GSstart_k_y;
	  if (iyy > nyhalf)
	    iyy -= MyGrids[0].GSglobal_y;
	  ky  = kynorm*iyy;

	  for (local_x = 0; local_x <= nxhalf; local_x++)
	    {
	      ixx = local_x;
	      kx  = kxnorm*ixx;

              k_squared  = kx*kx + ky*ky + kz*kz;
	      k_module   = sqrt(k_squared);
	      k_physical = k_module / params.InterPartDist;

	      /* corresponding index of real part in vector (imaginary in index + 1) */
	      index = 1 + 2*local_x + (MyGrids[0].GSglobal_x+MyGrids[0].off)
		*(local_z + local_y* MyGrids[0].GSglobal_z);

	      if (k_squared > 0)
		delta = sqrt(norm * PowerSpectrum(k_physical) *
			     exp(-pow(k_module/MyGrids[0].upper_k_cutoff,16.0) ) );
	      else
		delta = 0.0;

	      (cvector_fft[0][index/2])[0] *= delta;
	      (cvector_fft[0][index/2])[1] *= delta;

	    }
	}
    }

  write_from_cvector(0,kdensity[0]);

  return 0;

}

#endif
```


----- FILE: src/variables.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "pinocchio.h"

int ThisTask,NTasks;

//int            pfft_flags_c2r, pfft_flags_r2c;
MPI_Comm        FFT_Comm;
internal_data   internal;

char *main_memory, *wheretoplace_mycat;
product_data *products, *frag;
unsigned int **seedtable;  // QUESTO RIMANE?
unsigned int   *cubes_ordering;
double **kdensity;
double **density;
double ***first_derivatives;
double ***second_derivatives;
#ifdef TWO_LPT
double *kvector_2LPT;
double *source_2LPT;
#ifdef THREE_LPT
double *kvector_3LPT_1,*kvector_3LPT_2;
double *source_3LPT_1,*source_3LPT_2;
#endif
#endif

double Rsmooth;
smoothing_data Smoothing;

grid_data *MyGrids;
int Ngrids;

pfft_complex **cvector_fft;
double **rvector_fft;

param_data params={0};
output_data outputs;
subbox_data subbox;
#ifdef PLC
plc_data plc;
plcgroup_data *plcgroups;
#endif

cputime_data cputime={0.0};

int WindowFunctionType;

group_data *groups;

char date_string[25];

int *frag_pos,*indices,*indicesY,*sorted_pos,*group_ID,*linking_list;
unsigned int *frag_map, *frag_map_update;
int map_to_be_used;
double f_m, f_rm, espo, f_a, f_ra, f_200, sigmaD0;

gsl_integration_workspace * workspace;
gsl_rng *random_generator;

mf_data mf;

gsl_spline **SPLINE;
gsl_interp_accel **ACCEL;
#if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC)
gsl_spline **SPLINE_INVGROW;
gsl_interp_accel **ACCEL_INVGROW;
#endif

#ifdef MOD_GRAV_FR
double H_over_c;
#endif

memory_data memory;

int ngroups;
extern pos_data obj, obj1, obj2;
ScaleDep_data ScaleDep;
```


----- FILE: src/write_halos.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

#define FTOZ(A) (A > 0.0 ? A - 1 : A);

int NSlices = 1;

int compute_mf(int iout)
{
	/* computes the mass function of the groups and the analytic mass function */

	int i, ibin;
	char filename[LBLENGTH], lab1[10], lab2[10];
	double amass, x, dm, a, a1, a2, a3, m, D, mx, r, massvar, sigma, ni;
	FILE *file;

	D = GrowingMode(outputs.z[iout], params.k_for_GM);

	/* sets counters to zero */
	for (i = 0; i < mf.NBIN; i++)
	{
		mf.ninbin_local[i] = 0;
		mf.massinbin_local[i] = 0.0;
		mf.ninbin[i] = 0;
		mf.massinbin[i] = 0.0;
	}

	/* constructs the histograms */
	for (i = FILAMENT + 1; i <= ngroups; i++)
		if (groups[i].point >= 0 && groups[i].good && groups[i].Mass >= params.MinHaloMass)
		{
			amass = groups[i].Mass * params.ParticleMass;
			ibin = (int)((log10(amass) - mf.mmin) / DELTAM);
			if (ibin < 0 || ibin >= mf.NBIN)
				continue;
			mf.ninbin_local[ibin]++;
			mf.massinbin_local[ibin] += amass;
		}

	/* sums counters over all tasks */
	MPI_Reduce(mf.ninbin_local, mf.ninbin, mf.NBIN, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(mf.massinbin_local, mf.massinbin, mf.NBIN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Writes result */
	if (!ThisTask)
	{
		sprintf(filename, "pinocchio.%6.4f.%s.mf.out", outputs.z[iout], params.RunFlag);
		printf("[%s] Writing mass function into file  %s\n", fdate(), filename);
		file = fopen(filename, "w");

		fprintf(file, "# Mass function for redshift %f\n", outputs.z[iout]);
		if (params.OutputInH100)
		{
			strcpy(lab1, "/h");
			strcpy(lab2, "h^4");
		}
		else
		{
			strcpy(lab1, "");
			strcpy(lab2, "");
		}

		fprintf(file, "# 1) mass (Msun%s)\n", lab1);
		fprintf(file, "# 2) n(m) (Mpc^-3 Msun^-1 %s)\n", lab2);
		fprintf(file, "# 3) upper +1-sigma limit for n(m) (Mpc^-3 Msun^-1 %s)\n", lab2);
		fprintf(file, "# 4) lower -1-sigma limit for n(m) (Mpc^-3 Msun^-1 %s)\n", lab2);
		fprintf(file, "# 5) number of halos in the bin\n");
		switch (params.AnalyticMassFunction)
		{
		case 0:
			fprintf(file, "# 6) analytical n(m) from Press & Schechter 1974\n");
			break;
		case 1:
			fprintf(file, "# 6) analytical n(m) from Sheth & Tormen 2001\n");
			break;
		case 2:
			fprintf(file, "# 6) analytical n(m) from Jenkins et al. 2001\n");
			break;
		case 3:
			fprintf(file, "# 6) analytical n(m) from Warren et al. 2006\n");
			break;
		case 4:
			fprintf(file, "# 6) analytical n(m) from Reed et al. 2007\n");
			break;
		case 5:
			fprintf(file, "# 6) analytical n(m) from Crocce et al. 2010\n");
			break;
		case 6:
			fprintf(file, "# 6) analytical n(m) from Tinker et al. 2010\n");
			break;
		case 7:
			fprintf(file, "# 6) analytical n(m) from Courtin et al. 2010\n");
			break;
		case 8:
			fprintf(file, "# 6) analytical n(m) from Angulo et al. 2012\n");
			break;
		case 9:
			fprintf(file, "# 6) analytical n(m) from Watson et al. 2013\n");
			break;
		case 10:
			fprintf(file, "# 6) analytical n(m) from Crocce et al. 2010, universal\n");
			break;
		default:
			break;
		}
		fprintf(file, "#\n");

		for (i = 0; i < mf.NBIN; i++)
		{
			x = mf.mmin + (i + 0.5) * DELTAM;
			m = pow(10., x);
			dm = params.ParticleMass *
				 (double)((int)(pow(10., mf.mmin + (i + 1) * DELTAM) / params.ParticleMass) -
						  (int)(pow(10., mf.mmin + i * DELTAM) / params.ParticleMass));
			if (dm > 0.0)
			{
				a = (double)mf.ninbin[i] / mf.vol / dm;
				a1 = ((double)mf.ninbin[i] + sqrt((double)mf.ninbin[i])) / mf.vol / dm;
				a2 = ((double)mf.ninbin[i] - sqrt((double)mf.ninbin[i])) / mf.vol / dm;
			}
			else
				a = a1 = a2 = 0.0;

			if (mf.ninbin[i] > 1)
				mx = mf.massinbin[i] / (double)mf.ninbin[i];
			else
				mx = m;
			a3 = AnalyticMassFunction(mx, outputs.z[iout]);

			r = SizeForMass(mx);
			massvar = MassVariance(r) * D * D;
			sigma = sqrt(massvar);
			ni = 1.686 / sigma;

			fprintf(file,
					" %15.8g %15.8g %15.8g %15.8g   %10d  %15.8g    %15.8g\n",
					mx * mf.hfactor,
					a / mf.hfactor4,
					a1 / mf.hfactor4,
					a2 / mf.hfactor4,
					mf.ninbin[i],
					a3 / mf.hfactor4,
					ni);
		}
		fclose(file);
	}

	// LEVARE
	/*   if (iout==1) */
	/*     { */
	/*       FILE *ff=fopen("mydump","w"); */
	/*       for (int u=0; u<subbox.Nstored; u++) */
	/* #ifdef RECOMPUTE_DISPLACEMENTS */
	/* 	fprintf(ff," %d %d %d   %f %f %f\n", */
	/* 		u,frag_pos[u],sorted_pos[u],frag[u].Fmax,frag[u].Vel[0],frag[u].Vel_prev[0]); */
	/* #else */
	/* 	fprintf(ff," %d %d %d   %f %f %f\n", */
	/* 		u,frag_pos[u],sorted_pos[u],frag[u].Fmax,frag[u].Vel[0],0.0); */
	/* #endif	 */
	/*       fclose(ff); */
	/*     } */

	/* Bye! */
	return 0;
}

/*
   GENERAL COMMUNICATION SCHEME

   The code wants to write a catalog splitted in params.NumFiles
   Each file will be written by this number of tasks
	  NTasksPerFile=NTasks/params.NumFiles;

   Task ThisTask will write on file
	  ThisFile=ThisTask/NTasksPerFile;

   Among the tasks that access the file, one task will collect information
   and write it on the file. This is:
	  collector=ThisFile*NTasksPerFile;

   1: initialization of whatever

   2: each task constructs its own catalog (whatever it is)

   3: collector task opens the file and writes the header

   4: collector task writes the catalog on the file

   5: loop on all other tasks that write on the same file

	 5a: the task sends its catalog to collector task
	 5b: collector task receives the catalog

	 5c: collector task writes the catalog on the file

   6: collector task closes the file

*/

int write_catalog(int iout)
{
	/* Writes the group catalogues */

	int igood, i, ngood, nhalos, j;
	double hfactor, GGrid[3], SGrid[3];
	char filename[LBLENGTH], labh[3];
	int NTasksPerFile, collector, itask, next, ThisFile;
	catalog_data *mycat;
	FILE *file;
	MPI_Status status;
	int idummy;
	pos_data obj1;

	/* ordering of coordinates to accomodate for rotation caused by fft ordering */

	if (params.DoNotWriteCatalogs)
	{
		if (!ThisTask)
			printf("Halo catalog at z=%f will not be written\n", outputs.z[iout]);
		return 0;
	}

	GGrid[0] = (double)MyGrids[0].GSglobal[_x_];
	GGrid[1] = (double)MyGrids[0].GSglobal[_y_];
	GGrid[2] = (double)MyGrids[0].GSglobal[_z_];
	SGrid[0] = (double)subbox.stabl[_x_];
	SGrid[1] = (double)subbox.stabl[_y_];
	SGrid[2] = (double)subbox.stabl[_z_];

	NTasksPerFile = NTasks / params.NumFiles;
	ThisFile = ThisTask / NTasksPerFile;
	collector = ThisFile * NTasksPerFile;

	/* output in H100 or in Htrue */
	if (params.OutputInH100)
		hfactor = params.Hubble100;
	else
		hfactor = 1.0;

	/* each processor builds the catalogue */
	nhalos = 0;
	for (i = FILAMENT + 1, ngood = 0; i <= ngroups; i++)
		if (groups[i].point >= 0 && groups[i].good &&
			groups[i].Mass >= params.MinHaloMass)
			ngood++;

	/* space to store catalogs */
	mycat = (catalog_data *)wheretoplace_mycat;
	if (ngood * sizeof(catalog_data) > subbox.PredNpeaks * sizeof(histories_data))
	{
		printf("ERROR on task %d: surprisingly, memory reserved to mycat is insufficient in write_catalog\n", ThisTask);
		fflush(stdout);
		return 1;
	}

	if (ngood)
	{
		for (i = FILAMENT + 1, igood = 0; i <= ngroups; i++)
			if (groups[i].point >= 0 && groups[i].good &&
				groups[i].Mass >= params.MinHaloMass)
			{
				set_obj(i, outputs.F[iout], &obj1);
				set_obj_vel(i, outputs.F[iout], &obj1);
				mycat[igood].name = groups[i].name;
#ifndef LIGHT_OUTPUT
				mycat[igood].n = groups[i].Mass;
#endif
				mycat[igood].M = groups[i].Mass * params.ParticleMass * hfactor;
				for (j = 0; j < 3; j++)
				{
					/* q and x are in sub-box coordinates, they are transformed
					   to global box coordinates (forcing PBCs) */
					mycat[igood].q[j] = groups[i].Pos[j] + SGrid[j];
					if (mycat[igood].q[j] < 0)
						mycat[igood].q[j] += GGrid[j];
					if (mycat[igood].q[j] >= GGrid[j])
						mycat[igood].q[j] -= GGrid[j];
					mycat[igood].q[j] *= params.InterPartDist * hfactor;

					/* displacement is done up to ORDER_FOR_CATALOG */
					mycat[igood].x[j] = q2x(j, &obj1, subbox.pbc[j], (double)subbox.Lgwbl[j], ORDER_FOR_CATALOG) + SGrid[j];
					if (mycat[igood].x[j] < 0)
						mycat[igood].x[j] += GGrid[j];
					if (mycat[igood].x[j] >= GGrid[j])
						mycat[igood].x[j] -= GGrid[j];
					mycat[igood].x[j] *= params.InterPartDist * hfactor;
					mycat[igood].v[j] = vel(j, &obj1);
				}
				igood++;
			}
	}

	/* The collector task opens the file and writes its catalogue */
	if (ThisTask == collector)
	{
		if (params.NumFiles > 1)
			sprintf(filename, "pinocchio.%6.4f.%s.catalog.out.%d",
					outputs.z[iout], params.RunFlag, ThisFile);
		else
			sprintf(filename, "pinocchio.%6.4f.%s.catalog.out",
					outputs.z[iout], params.RunFlag);

		if (!ThisTask)
			printf("[%s] Opening file %s\n", fdate(), filename);

		file = fopen(filename, "w");

		if (params.CatalogInAscii)
		{
			fprintf(file, "# Group catalog for redshift %f and minimal mass of %d particle%s\n",
					outputs.z[iout], params.MinHaloMass, (params.MinHaloMass == 1 ? "" : "s"));

			if (params.OutputInH100)
				strcpy(labh, "/h");
			else
				strcpy(labh, "");

			fprintf(file, "#    1) group ID\n");
			fprintf(file, "#    2) group mass (Msun%s)\n", labh);
			fprintf(file, "# 3- 5) initial position (Mpc%s)\n", labh);
			fprintf(file, "# 6- 8) final position (Mpc%s)\n", labh);
			fprintf(file, "# 9-11) velocity (km/s)\n");
#ifndef LIGHT_OUTPUT
			fprintf(file, "#   12) number of particles\n");
#endif
			fprintf(file, "#\n");

			if (ngood)
				for (igood = 0; igood < ngood; igood++)
#ifndef LIGHT_OUTPUT
					fprintf(file, " %12Lu %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
#else
					fprintf(file, " %12Lu %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
#endif
							mycat[igood].name,
							mycat[igood].M,
							mycat[igood].q[0], mycat[igood].q[1], mycat[igood].q[2],
							mycat[igood].x[0], mycat[igood].x[1], mycat[igood].x[2],
							mycat[igood].v[0], mycat[igood].v[1], mycat[igood].v[2]
#ifndef LIGHT_OUTPUT
							,
							mycat[igood].n
#endif
					);
		}
		else
		{
			idummy = 2 * sizeof(int);
			fwrite(&idummy, sizeof(int), 1, file);
			fwrite(&NTasksPerFile, sizeof(int), 1, file);
			idummy = sizeof(catalog_data);
			fwrite(&idummy, sizeof(int), 1, file);
			idummy = 2 * sizeof(int);
			fwrite(&idummy, sizeof(int), 1, file);

			idummy = sizeof(int);
			fwrite(&idummy, sizeof(int), 1, file);
			fwrite(&ngood, sizeof(int), 1, file);
			fwrite(&idummy, sizeof(int), 1, file);

			if (ngood)
			{
				idummy = ngood * sizeof(catalog_data);
				fwrite(&idummy, sizeof(int), 1, file);
				for (igood = 0; igood < ngood; igood++)
					fwrite(mycat + igood, sizeof(catalog_data), 1, file);
				fwrite(&idummy, sizeof(int), 1, file);
			}
		}

		nhalos += ngood;
	}

```


----- FILE: src/write_mass_maps.c -----
```text
/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

int write_mass_maps(double z_start, double z_end)
{
  (void)z_start;
  (void)z_end;
  if (ThisTask == 0)
    printf("[%s] MASS_MAPS: diagnostics-only skeleton (positions + summaries).\n", fdate());
  return 0;
}

#endif /* MASS_MAPS */
```


----- FILE: src/write_snapshot.c -----
```text
/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

//#define POS_IN_KPC   /* with this directive on, positions will be in kpc/h */


/* This file contains functions to handle gadget-like snapshots
   (format 2 with INFO block) with information on all particles
   Possible snapshots are:
   - timeless snapshot -- gives ID and LPT fields for all particles, plus Fmax, Rmax (if required)
     and the ZACC time at which the particle enters a halo
   - LPT snapshot -- gives ID, position and velocity of each particle, according
     to LPT at a given order
   - density -- it writes the ID and the linear density field
*/


#ifdef SNAPSHOT

static unsigned long long largest32 = (unsigned)1<<31;

#ifdef LONGIDS
#define MYIDTYPE unsigned long long int
#else
#define MYIDTYPE unsigned int
#endif

/* gadget header */
typedef struct
{
  unsigned NPart[6];
  double   Mass[6];
  double   Time;
  double   RedShift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned NPartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  int      flag_stellarage;
  int      flag_metals;
  unsigned npartTotalHighWord[6];
  int      flag_entropy_instead_u;
  int      flag_metalcooling;
  int      flag_stellarevolution;
  char     fill[52];  /* filling to 256 Bytes */
} SnapshotHeader_data;

/* information for the INFO block */
typedef struct
{
  char name[4], type[8];
  int ndim, active[6];
  size_t sizeof_type;
  void *data;
} Block_data;

Block_data *InfoBlock;
int NBlocks, NextBlock;

/* for vector blocks */
typedef struct
{
  float axis[3];
} AuxStruct;

MPI_Status status;

void WriteBlockName(FILE *,unsigned long long, char*);
void my_strcpy(char *,char *, int);
int write_header();
int write_block(Block_data);
void free_block(Block_data);
int add_to_info(Block_data);
int write_info_block(void);
int initialize_ID(Block_data*);
int initialize_FMAX(Block_data*);
int initialize_RMAX(Block_data*);
int initialize_ZEL(Block_data*);
#ifdef TWO_LPT
int initialize_2LPT(Block_data*);
#ifdef THREE_LPT
int initialize_3LPT_1(Block_data*);
int initialize_3LPT_2(Block_data*);
#endif
#endif
int initialize_ZACC(Block_data*);
int initialize_GRUP(Block_data*);
int initialize_POS(Block_data*);
int initialize_VEL(Block_data*);
int initialize_density(int, Block_data*);
void set_point_timedep(double);

char filename[LBLENGTH];
FILE *file;
int NTasksPerFile,collector,ThisFile,NPartInFile,myNpart,myiout;
int *Npart_array;

/* this is the redshift at which an LPT snapshot is written */
int myiout;

/* this is used to shift particles to the final position */
//pos_data myobj;

double redshift;

int write_LPT_snapshot()
{
  /* write positions of all particles obtained with LPT */

  Block_data block;
  //myiout=0;
  //set_point_timedep(outputs.z[myiout]);
  redshift = outputs.z[0];

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%6.4f.%s.LPT_snapshot.out",outputs.z[myiout],params.RunFlag);

  /* allocates structure to handle the INFO block */
  NBlocks=4;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of POS */
  if (initialize_POS(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_VEL(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(Npart_array);
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;

}


int write_timeless_snapshot()
{
  /* writes the fmax products as in a snapshot format */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.t_snapshot.out",params.RunFlag);

  myiout=outputs.n-1;

  /* allocates structure to handle the INFO block */
#ifdef TWO_LPT
#ifdef THREE_LPT
  NBlocks=8;
#else
  NBlocks=6;
#endif
#else
  NBlocks=4;
#endif
#ifdef ADD_RMAX_TO_SNAPSHOT
  ++NBlocks;
#endif

  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing products in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef ADD_RMAX_TO_SNAPSHOT
  /* writing of RMAX */
  if (initialize_RMAX(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);
#endif

  /* writing of FMAX */
  if (initialize_FMAX(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_ZEL(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef TWO_LPT
  /* writing of VEL2 */
  if (initialize_2LPT(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef THREE_LPT
  /* writing of VL31 */
  if (initialize_3LPT_1(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VL32 */
  if (initialize_3LPT_2(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);
#endif
#endif

  /* writing of ZACC */
  if (initialize_ZACC(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of particle GROUP_ID */
  if (initialize_GRUP(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

int write_density(int ThisGrid)
{
  /* writes the fmax products as in a snapshot format */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.density%d.out",params.RunFlag,ThisGrid);

  /* allocates structure to handle the INFO block */
  NBlocks=2;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing density in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of density */
  if (initialize_density(ThisGrid, &block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;

}


int write_header()
```