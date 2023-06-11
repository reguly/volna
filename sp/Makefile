#
# The following environment variables should be predefined:
#
#OP2_INSTALL_PATH = 
#HDF5_INSTALL_PATH =
#MPI_INSTALL_PATH =
#PARMETIS_INSTALL_PATH =
#PTSCOTCH_INSTALL_PATH =
#CUDA_INSTALL_PATH =
#OP2_COMPILER =
#SLOPE_INSTALL_PATH=
#METIS_INSTALL_PATH=

#
# set paths for header files and libraries
#

OP2_INC = -I$(OP2_INSTALL_PATH)/include -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB = -L$(OP2_INSTALL_PATH)/lib -L$(OP2_INSTALL_PATH)/c/lib

HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz

MPI_INC = -I$(MPI_INSTALL_PATH)/include

CUDA_INC	= -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB	= -L$(CUDA_INSTALL_PATH)/lib64

ifdef SLOPE_INSTALL_PATH
  SLOPE_INC = -I$(SLOPE_INSTALL_PATH)/sparsetiling/include
  SLOPE_LIB = -L$(SLOPE_INSTALL_PATH)/lib
endif

ifdef METIS_INSTALL_PATH
  METIS_INC = -I$(METIS_INSTALL_PATH)/include
  METIS_LIB = -L$(METIS_INSTALL_PATH)/lib
endif

ifdef PARMETIS_INSTALL_PATH
PARMETIS_INCLUDE_PATH ?= $(PARMETIS_INSTALL_PATH)
PARMETIS_LIB_PATH ?= $(PARMETIS_INSTALL_PATH)
PARMETIS_VER?=4
ifeq ($(PARMETIS_VER),4)
  PARMETIS_INC = -I$(PARMETIS_INCLUDE_PATH)/include -DHAVE_PARMETIS -DPARMETIS_VER_4
  PARMETIS_LIB = -L$(PARMETIS_LIB_PATH)/lib -lparmetis -lmetis
else
  PARMETIS_INC = -I$(PARMETIS_INCLUDE_PATH)/ -DHAVE_PARMETIS
  PARMETIS_LIB = -L$(PARMETIS_LIB_PATH)/ -lparmetis -lmetis
endif
endif

ifdef PTSCOTCH_INSTALL_PATH
PTSCOTCH_INC 	= -I$(PTSCOTCH_INSTALL_PATH)/include -DHAVE_PTSCOTCH
PTSCOTCH_INC   = -I$(PTSCOTCH_INSTALL_PATH)/include
PTSCOTCH_LIB 	= -L$(PTSCOTCH_INSTALL_PATH)/lib/ -lptscotch -lscotch \
	-L$(PTSCOTCH_INSTALL_PATH)/lib/ -lptscotcherr
endif

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
	#CPPFLAGS = -O0
	CPPFLAGS = -O2 -DUNIX -Wall #-DDEBUG
	OMPFLAGS = -fopenmp
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpic++
	MPIFLAGS = $(CPPFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
#	CPPFLAGS = -g -O0
	#CPPFLAGS = -O3 -qopt-report -xHost -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
        CPPFLAGS = -O3 -xHost -inline-forceinline -qopenmp -DMPICH_IGNORE_CXX_SEEK -restrict -fno-alias -qopt-report -parallel -std=c++11 #-qopt-zmm-usage=high
#	CPPFLAGS = -g -O0 -vec-report -xSSE4.2 -parallel
	OMPFLAGS = -qopenmp #-openmp-report2
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpicxx
	MPIFLAGS = $(CPPFLAGS)
else
ifeq ($(OP2_COMPILER),intel-sycl)
  ifdef DEBUG
    CPPFLAGS  = -g -O0
  else
    CPPFLAGS  = -g -O3
  endif
  CXX       = g++
  SYCLCXX   = dpcpp
  CXXFLAGS  = $(CCFLAGS)
  MPICXX    = $(MPICPP_PATH)
  MPICPP    = mpicxx
  MPISYCLCXX = $(MPICXX)
  MPIFLAGS  = $(CXXFLAGS)
  NVCCFLAGS = -ccbin=$(NVCC_HOST_COMPILER)
  SYCL_LIB   = #-L$(SYCL_INSTALL_PATH)/lib -lOpenCL
  SYCL_FLAGS = -std=c++17 -fsycl -fsycl-targets=spir64_x86_64 -Xs "-march=avx512" -Ofast -ffast-math -march=native -fp-model fast -fno-alias -std=c++17 -qopenmp #-I$(SYCL_INSTALL_PATH)/include -I$(SYCL_INSTALL_PATH)/include #intel sycl
  #SYCL_FLAGS = -std=c++11 -fsycl #intel sycl
else
ifeq ($(OP2_COMPILER),cray)
  ifdef DEBUG
    CCFLAGS  = -g -O0
  else
    CCFLAGS  = -Ofast
  endif
  CXX       = CC
  HIPCXX   = CC -x hip
  CXXFLAGS  = $(CCFLAGS)
  MPICXX    = $(MPICPP_PATH)
  MPIHIPCXX    = $(MPICPP_PATH)
  MPICPP    = CC
  MPIFLAGS  = $(CXXFLAGS)
  HIPFLAGS = $(CXXFLAGS) --offload-arch=$(HIP_ARCH) -munsafe-fp-atomics
else
print:
	@echo "unrecognized value for OP2_COMPILER"
endif
endif
endif
endif

ifdef PROFILE_ITT
OP2_LIB += -L$(VTUNE_PROFILER_DIR)/lib64 -littnotify
CXXFLAGS += -DPROFILE_ITT -I$(VTUNE_PROFILER_DIR)/include
endif
#VOLNA_SRC =


#
# set flags for NVCC compilation and linking
#


NVCCFLAGS	= -arch=sm_35 -Xptxas=-v -Dlcm=ca -use_fast_math -O3 -m64 -g #-G

#
# master to make all versions
#

all: clean volna volna_slope volna_openmp volna_cuda volna_mpi volna_mpi_openmp volna_mpi_cuda

volna: volna.cpp volna_util.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_writeVTK.cpp volna_simulation.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) $(PARMETIS_INC) $(PTSCOTCH_INC) volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_simulation.cpp volna_util.cpp volna_writeVTK.cpp $(PTSCOTCH_LIB) -lscotch $(PARMETIS_LIB) $(HDF5_INC) $(OP2_INC) $(HDF5_LIB) $(OP2_LIB) -lop2_seq -lop2_hdf5 -o volna
else
	$(CPP) $(CPPFLAGS) $(PARMETIS_INC) $(PTSCOTCH_INC) volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_simulation.cpp volna_util.cpp volna_writeVTK.cpp $(PTSCOTCH_LIB) -lscotch $(PARMETIS_LIB) $(HDF5_INC) $(OP2_INC) $(HDF5_LIB) $(OP2_LIB) -lop2_seq -lop2_hdf5 -o volna
endif

volna_slope: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) -DSLOPE volna_op.cpp \
	volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp \
	volna_simulation_op.cpp volna_event.cpp \
	slope/volna_kernels.cpp \
	$(HDF5_INC) $(OP2_INC) $(SLOPE_INC) $(METIS_INC) $(HDF5_LIB) $(OP2_LIB) $(SLOPE_LIB) $(METIS_LIB) \
	-lop2_openmp -lop2_hdf5 -lslope -lmetis -o volna_slope
else
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) -DSLOPE  volna_op.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp volna_event.cpp openmp/volna_kernels.cpp \
	$(HDF5_INC) $(OP2_INC) $(SLOPE_INC) $(METIS_INC) $(HDF5_LIB) $(OP2_LIB) $(SLOPE_LIB) $(METIS_LIB) \
	-lop2_openmp -lop2_hdf5 -lslope -lmetis -o volna_slope
endif

volna_openmp: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) volna_op.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp volna_event.cpp openmp/volna_kernels.cpp $(HDF5_INC) $(OP2_INC) $(HDF5_LIB) $(OP2_LIB) -lop2_openmp -lop2_hdf5 -o volna_openmp
else
	$(CPP) $(CPPFLAGS) $(OMPFLAGS)  volna_op.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp volna_event.cpp openmp/volna_kernels.cpp $(HDF5_INC) $(OP2_INC) $(HDF5_LIB) $(OP2_LIB) -lop2_openmp -lop2_hdf5 -o volna_openmp
endif


#
# CUDA version using kernel files generated by op2.m
#

volna_cuda:	volna_op.cpp cuda/volna_kernels_cu.o volna_simulation_op.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(VAR) $(CPPFLAGS) volna_op.cpp volna_simulation_op.cpp volna_event.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp cuda/volna_kernels_cu.o \
	$(PTSCOTCH_INC) $(PTSCOTCH_LIB) -lscotch $(PARMETIS_LIB) $(PARMETIS_INC) $(CUDA_INC) $(OP2_INC) $(HDF5_INC) \
	$(OP2_LIB) $(CUDA_LIB) -lcudart -lop2_cuda -lop2_hdf5 $(HDF5_LIB) -o volna_cuda
else
	$(MPICPP) $(VAR) $(CPPFLAGS) volna_op.cpp volna_simulation_op.cpp volna_event.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp cuda/volna_kernels_cu.o \
	$(PTSCOTCH_INC) $(PTSCOTCH_LIB) -lscotch $(PARMETIS_LIB) $(PARMETIS_INC) $(CUDA_INC) $(OP2_INC) $(HDF5_INC) \
	$(OP2_LIB) $(CUDA_LIB) -lcudart -lop2_cuda -lop2_hdf5 $(HDF5_LIB) -o volna_cuda
endif

volna_hip:	volna_op.cpp hip/volna_kernels.o volna_simulation_op.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp Makefile
	$(MPICPP) $(VAR) $(CPPFLAGS) volna_op.cpp volna_simulation_op.cpp volna_event.cpp volna_init_op.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp hip/volna_kernels.o \
	$(PTSCOTCH_INC) $(PTSCOTCH_LIB) -lscotch $(PARMETIS_LIB) $(PARMETIS_INC) $(OP2_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_hip -lop2_hdf5 $(HDF5_LIB) -o volna_hip


cuda/volna_kernels_cu.o:	cuda/volna_kernels.cu \
	EvolveValuesRK3_1.h EvolveValuesRK3_2.h  EvolveValuesRK3_3.h EvolveValuesRK3_4.h applyConst.h getMaxElevation.h getTotalVol.h \
	initBathymetry_formula.h initBathymetry_update.h initBore_select.h initEta_formula.h initGaussianLandslide.h \
	initU_formula.h initV_formula.h computeGradient.h computeFluxes.h limiter.h SpaceDiscretization.h NumericalFluxes.h simulation_1.h \
	values_operation2.h cuda/applyConst_kernel.cu cuda/EvolveValuesRK3_1_kernel.cu  cuda/EvolveValuesRK3_2_kernel.cu\
	cuda/EvolveValuesRK3_3_kernel.cu cuda/EvolveValuesRK3_4_kernel.cu cuda/computeGradient_kernel.cu \
	cuda/computeFluxes_kernel.cu cuda/limiter_kernel.cu cuda/SpaceDiscretization_kernel.cu cuda/NumericalFluxes_kernel.cu \
	cuda/simulation_1_kernel.cu \
	cuda/values_operation2_kernel.cu Makefile

	nvcc  $(VAR) $(INC) $(NVCCFLAGS) $(OP2_INC) $(HDF5_INC) -I$(MPI_INC) -c -o cuda/volna_kernels_cu.o cuda/volna_kernels.cu

hip/volna_kernels.o:	hip/volna_kernels.cpp \
	EvolveValuesRK3_1.h EvolveValuesRK3_2.h  EvolveValuesRK3_3.h EvolveValuesRK3_4.h applyConst.h getMaxElevation.h getTotalVol.h \
	initBathymetry_formula.h initBathymetry_update.h initBore_select.h initEta_formula.h initGaussianLandslide.h \
	initU_formula.h initV_formula.h computeGradient.h computeFluxes.h limiter.h SpaceDiscretization.h NumericalFluxes.h simulation_1.h \
	values_operation2.h hip/applyConst_kernel.cpp hip/EvolveValuesRK3_1_kernel.cpp  hip/EvolveValuesRK3_2_kernel.cpp\
	hip/EvolveValuesRK3_3_kernel.cpp hip/EvolveValuesRK3_4_kernel.cpp hip/computeGradient_kernel.cpp \
	hip/computeFluxes_kernel.cpp hip/limiter_kernel.cpp hip/SpaceDiscretization_kernel.cpp hip/NumericalFluxes_kernel.cpp \
	hip/simulation_1_kernel.cpp \
	hip/values_operation2_kernel.cpp Makefile

	$(HIPCXX)  $(VAR) $(INC) $(HIPFLAGS) $(OP2_INC) $(HDF5_INC) -I$(MPI_INC) -c -o hip/volna_kernels.o hip/volna_kernels.cpp


volna_mpi: volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_writeVTK.cpp volna_simulation.cpp volna_util.cpp Makefile
	$(MPICPP) $(MPIFLAGS) $(PARMETIS_INC) $(PTSCOTCH_INC) volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_writeVTK.cpp volna_simulation.cpp volna_util.cpp $(OP2_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o volna_mpi

volna_sycl: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
	$(MPICPP) $(SYCL_FLAGS) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(OP2_INC) $(HDF5_INC) \
	volna_op.cpp volna_init_op.cpp volna_event.cpp sycl/volna_kernels.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp -lm $(OP2_INSTALL_PATH)/c/lib/libop2_sycl.a $(OP2_LIB) -lop2_hdf5 \
	$(HDF5_LIB) -o volna_sycl

volna_mpi_sycl: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
	$(MPICPP) $(SYCL_FLAGS) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(OP2_INC) $(HDF5_INC) \
	$(PARMETIS_INC) $(PTSCOTCH_INC) \
	volna_op.cpp volna_init_op.cpp volna_event.cpp sycl/volna_kernels.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp -lm $(OP2_INSTALL_PATH)/c/lib/libop2_mpi_sycl.a \
	$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o volna_mpi_sycl

volna_mpi_vec: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
	$(MPICPP) $(VAR) -DVECTORIZE $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(OP2_INC) $(HDF5_INC) \
	$(PARMETIS_INC) $(PTSCOTCH_INC) \
	volna_op.cpp volna_init_op.cpp volna_event.cpp vec/volna_veckernels.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp -lm $(OP2_LIB) -lop2_mpi \
	$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o volna_mpi_vec

volna_mpi_genseq: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
	$(MPICPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(OP2_INC) $(HDF5_INC) \
	$(PARMETIS_INC) $(PTSCOTCH_INC) \
	volna_op.cpp volna_init_op.cpp volna_event.cpp seq/volna_seqkernels.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp -lm $(OP2_LIB) -lop2_mpi \
	$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o volna_mpi_genseq

volna_mpi_openmp: volna_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp Makefile
	$(MPICPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(OP2_INC) $(HDF5_INC) \
	$(PARMETIS_INC) $(PTSCOTCH_INC) \
	volna_op.cpp volna_init_op.cpp volna_event.cpp openmp/volna_kernels.cpp volna_output_op.cpp volna_writeVTK.cpp volna_simulation_op.cpp -lm $(OP2_LIB) -lop2_mpi \
	$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o volna_mpi_openmp

volna_mpi_cuda: volna_op.cpp volna_simulation_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp cuda/volna_kernels_mpi_cu.o Makefile
	$(MPICPP) $(MPIFLAGS) volna_op.cpp volna_simulation_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp -lm cuda/volna_kernels_mpi_cu.o \
	$(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
	$(HDF5_LIB) $(CUDA_LIB) -lcudart -o volna_mpi_cuda

volna_mpi_hip: volna_op.cpp volna_simulation_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp hip/volna_kernels_mpi.o Makefile
	$(MPICPP) $(MPIFLAGS) volna_op.cpp volna_simulation_op.cpp volna_init_op.cpp volna_event.cpp volna_output_op.cpp volna_writeVTK.cpp volna_util.cpp -lm hip/volna_kernels_mpi.o \
	$(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_mpi_hip $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
	$(HDF5_LIB) -o volna_mpi_hip

cuda/volna_kernels_mpi_cu.o: cuda/volna_kernels.cu Makefile
	nvcc  $(INC) $(NVCCFLAGS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
	-c -o cuda/volna_kernels_mpi_cu.o cuda/volna_kernels.cu

hip/volna_kernels_mpi.o: hip/volna_kernels.cpp Makefile
	$(HIPCXX)  $(INC) $(HIPFLAGS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
	-c -o hip/volna_kernels_mpi.o hip/volna_kernels.cpp

#
# cleanup
#

clean:
	rm -f volna volna_slope volna_openmp volna_cuda volna_mpi volna_mpi_openmp volna_mpi_cuda *.o cuda/*.o 
