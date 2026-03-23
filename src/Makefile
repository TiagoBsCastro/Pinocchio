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


# =============================================================
#  User Configuration
# =============================================================

RUNDIR  = ./
EXEC    = pinocchio.x
DEBUG   = NO
OMP     = YES
SYSTYPE ?= "LeonardoDCGP"


# =============================================================
#  Feature Toggles  (preprocessor definitions)
# =============================================================
#
# Uncomment / comment the OPTIONS lines below to enable / disable
# compile-time features.  Only enabled flags produce code in the
# binary.
#
# --- LPT order -----------------------------------------------------------
#  -DTWO_LPT              Second-order Lagrangian perturbation theory.
#  -DTHREE_LPT            Third-order LPT (adds computational cost).
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# --- Collapse model (choose one) ----------------------------------------
#  -DELL_CLASSIC           Classic ellipsoidal collapse (default).
#  -DELL_SNG               Spherical-no-gravity approximation.
#  -DTABULATED_CT          Use tabulated collapse times from file.
OPTIONS += -DELL_CLASSIC

# --- Output and products ------------------------------------------------
#  -DSNAPSHOT              Populate products[] (z_acc, group_ID, etc.).
#                          Required by MASS_MAPS_FILTER_UNCOLLAPSED.
#  -DLONGIDS               Use 64-bit particle IDs.
#  -DLIGHT_OUTPUT          Reduce output footprint.
#  -DDOUBLE_PRECISION_PRODUCTS  Store products in double precision.
OPTIONS += -DLONGIDS

# --- Past Light-Cone & Mass Maps ----------------------------------------
#  -DPLC                   Enable past light-cone geometry.
#  -DMASS_MAPS             Build per-segment HEALPix maps.
#  -DHAVE_CFITSIO          Write FITS files (requires CFITSIO).
#  -DMASS_MAPS_FILTER_UNCOLLAPSED  Keep only Z_PLC > Z_ACC entries.
#  -DMASS_MAPS_THREAD_ACCUM       Per-thread HEALPix accumulators.
#  -DMASS_MAPS_CAP_DIAG    Cap prefix diagnostics (aperture QA).
#  -DMASS_MAPS_MATCH_DIAG  Extra cross-matching diagnostics (debug).
#  -DDISP_HISTOGRAM        Displacement magnitude histograms (QA).

# --- Scale-dependent growth & beyond LCDM --------------------------------
#  -DSCALE_DEPENDENT       Enable k-dependent growth (FFT per segment).
#  -DREAD_PK_TABLE         Read external P(k) table.
#  -DRECOMPUTE_DISPLACEMENTS  Maintain prev/curr segment fields.
#  -DREAD_HUBBLE_TABLE     Read H(z) table for precise distances.
OPTIONS += -DSCALE_DEPENDENT
OPTIONS += -DREAD_PK_TABLE
OPTIONS += -DRECOMPUTE_DISPLACEMENTS
OPTIONS += -DREAD_HUBBLE_TABLE

# --- f(R) modified gravity -----------------------------------------------
#  -DMOD_GRAV_FR           Enable f(R) modified gravity.
#  -DFR0=1.e-8             FR parameter value.

# --- Fragmentation -------------------------------------------------------
#  -DCLASSIC_FRAGMENTATION  Use classic fragmentation algorithm.
#OPTIONS += -DCLASSIC_FRAGMENTATION

# --- Distribution --------------------------------------------------------
#  -DUSE_DISTRIBUTE_ALLTOALL  Use MPI_Alltoall-based distribution.
OPTIONS += -DUSE_DISTRIBUTE_ALLTOALL

# --- Other ---------------------------------------------------------------
#  -DWHITENOISE            White-noise initial conditions.
#  -DNORADIATION           Neglect radiation component.


# =============================================================
#  OpenMP Configuration
# =============================================================

ifeq ($(OMP), YES)
OPTIONS  += -DUSE_FFT_THREADS
OMP_FLAG  = -fopenmp
OMP_FFTW  = -lfftw3_omp -lfftw3_threads
else
OMP_FLAG  =
OMP_FFTW  =
endif


# =============================================================
#  External Library Paths  (defaults, overridden per system)
# =============================================================

HEALPIX_ROOT ?= /usr
CFITSIO_ROOT ?= /usr

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio


# =============================================================
#  System-Specific Configuration
# =============================================================
#
# Set SYSTYPE above (or via the environment) to match your platform.
# Each block defines CC, CFLAGS, and library paths.

# ----------------------- mylaptop -------------------------

ifeq ($(SYSTYPE),"mylaptop")
CC          = mpicc
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
endif

# ----------------------- base -----------------------------

ifeq ($(SYSTYPE),"base")
CC          = mpicc
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

# ----------------------- Marconi100 -----------------------

ifeq ($(SYSTYPE),"Marconi100")
CC          = mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O2 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L$(HOME)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(HOME)/include
PFFT_LIBR   = -lpfft
PFFT_INCL   =
MPI_LIBR    = -L$(SMPI_ROOT)/lib -lmpi_ibm
MPI_INCL    = -I$(SMPI_ROOT)/include
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

# ----------------------- pleiadi --------------------------

ifeq ($(SYSTYPE),"pleiadi")
CC          = mpicc
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

# ----------------------- LeonardoBoost --------------------

ifeq ($(SYSTYPE),"LeonardoBoost")
CC          = nvc
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

# ----------------------- LeonardoDCGP --------------------

ifeq ($(SYSTYPE),"LeonardoDCGP")
CC          = mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)

FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW_INC)

GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)

PFFT_ROOT   = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/
PFFT_LIBR   = -L$(PFFT_ROOT)/lib -lpfft
PFFT_INCL   = -I$(PFFT_ROOT)/include

HEALPIX_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/Healpix_3.83/
CFITSIO_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio
endif


# =============================================================
#  Includes, Libraries & Compiler Options
# =============================================================

INC  = $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)
INC += $(HEALPIX_INCL) $(CFITSIO_INCL)

LIB  = -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR)
LIB += $(CFITSIO_LIBR) $(HEALPIX_LIBR)

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif


# =============================================================
#  Source Files
# =============================================================

OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o write_mass_maps.o cosmo.o


# =============================================================
#  Build Targets
# =============================================================

pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)

%.o: %.c Makefile pinocchio.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

clean:
	rm -f *.o *~ $(EXEC) run_planner


# =============================================================
#  Documentation  (requires doxygen)
# =============================================================
#
#   make docs         – build HTML docs under ../docs/build/html
#   make docs-clean   – remove generated docs

.PHONY: docs docs-clean

DOXYGEN ?= doxygen
DOCDIR  ?= ../docs

docs:
	@command -v $(DOXYGEN) >/dev/null 2>&1 || { \
	  echo "Error: doxygen not found. Install doxygen and re-run 'make docs'."; \
	  exit 1; \
	}
	@echo "Generating documentation with $$( $(DOXYGEN) --version )..."
	@cd $(DOCDIR) && $(DOXYGEN) Doxyfile
	@echo "Docs available at $(DOCDIR)/build/html/index.html"

docs-clean:
	@rm -rf $(DOCDIR)/build
