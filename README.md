Volna-OP2
=====

This is an extension to the OP2 port of the original Volna code.

## Installation
You need to install OP2 (source in [this](https://github.com/OP2/OP2-Common) GitHub repository)
 * You will need HDF5 (preferably compiled with MPI support)
 * You may not need the partitioners PT-Scotch/ParMetis if you do not plan to run in a distributed environment

Check out the Volna-OP2 code from [this](https://github.com/reguly/Volna) GitHub repository. The code comprises of two parts, volna2hdf5 which takes configuration files from the original Volna code and dumps all the necessary information to a h5 file which will be used by the second application volna-OP2.
 * Type 'make' in sp/volna2hdf5 - note that some warnings will show because of dependencies, you can safely ignore these
 * Source the OP2 build environment, then use the shared OP2 build rules in `sp`:

   ```bash
   source /path/to/OP2-Common/scripts/source_gnuz
   cd sp
   make volna_genseq
   ```

   `volna_genseq` is the current-translator sequential build. `make volna_seq`
   builds the original sources without code generation, while `make volna_openmp`,
   `make volna_cuda`, and `make volna_c_cuda` select the corresponding generated
   backend. MPI variants use the same naming convention, for example
   `volna_mpi_genseq`, `volna_mpi_openmp`, and `volna_mpi_cuda`. `make all`
   builds every backend enabled by the OP2 environment; `make clean` removes all
   build products and the `generated/` directory.

## Use
For all details and configuration options please see the documentation.

## Running the Code

Volna requires an input HDF5 file and an output-format argument: `0` for HDF5,
`1` for ASCII VTK, or `2` for binary VTK. Append `old-format` only for a
legacy bathymetry input. The supplied mesh can be found in 
[volna_30m2](https://users.itk.ppke.hu/~regiszo/volna_30m2.h5). 
The commands below use HDF5 output and write output files
to the current directory.

```bash
cd /path/to/Volna-OP2/sp
MESH=$PWD/volna_30m2.h5

./volna_seq "$MESH" 0 old-format
./volna_genseq "$MESH" 0 old-format
OMP_NUM_THREADS=32 ./volna_openmp "$MESH" 0 old-format
./volna_cuda "$MESH" 0 old-format
./volna_c_cuda "$MESH" 0 old-format
./volna_hip "$MESH" 0 old-format
./volna_c_hip "$MESH" 0 old-format

mpirun -np 64 ./volna_mpi_genseq "$MESH" 0 old-format
OMP_NUM_THREADS=4 mpirun -np 16 --map-by ppr:16:node:PE=4 ./volna_mpi_openmp "$MESH" 0 old-format
mpirun -np 4 ./volna_mpi_cuda "$MESH" 0 old-format
mpirun -np 4 ./volna_mpi_c_cuda "$MESH" 0 old-format
mpirun -np 4 ./volna_mpi_hip "$MESH" 0 old-format
mpirun -np 4 ./volna_mpi_c_hip "$MESH" 0 old-format
```

Build the selected executable first. `volna_seq`, `volna_genseq`, and
`volna_openmp` are usable with the configured CPU toolchain; CUDA/HIP commands
require the corresponding OP2 backend and compiler, and MPI commands require
the parallel HDF5 configuration.

To use volna-OP2 with the *.vln configuration files, first you have to use volna2hdf5, e.g.
 * ./volna2hdf5 gaussian_landslide.vln which will output a gaussian_landslide.h5 file
Afterwards, call volna-op2 with the above input file, e.g.:
 * ./volna_openmp gaussian_landslide.h5
 * when using the CUDA version we suggest adding "OP_PART_SIZE=128 OP_BLOCK_SIZE=128" to the execution line
 * if you use the InitBathymetry event with files, by default it will look for v3 files, so if your string is something like "bathy%i.txt" then it will look for a bathy_init.txt, which specifies the initial bathymetry for each cell center, a bathy_geom.txt which describes the coarse mesh which will define the bathymetry deformations, this file is structured as follows: first line number of cells and number of points, next #(number of cells) number of lines listing the indices of the points of each cell, then #(number of points) number of lines, each with 3 double precision values for the x and y and z coordinates of the points. Then, it will read the actual bathymetry deformations on each point on the coarse mesh, from a sequence of files, as defined by the description of the event

## Output
The files written by volna-OP2 may be different from the original VOLNA code, due to performance optimisations.
 * OutputLocation events are bundled together, if they all use the same timings, into a gauges.h5 file, which is compressed to save space. To read the file, you can use any hdf5 tool, there are two datasets written: /dims which has 2 integer fields, the first indicating the size in the contiguous direction (number of OutputLocation events + 1 for timestaps) and the second in the strided dimension (number of timestamps). Data is stored under /gauges, a 1D float array of size dims[0]*dims[1]. To get the timestamp at iteration T and the value for gauge N (indexed 1...dims[0]-1), access gauges[T*dims[0]] for the timestamp and gauges[T*dims[0]+N] - T has to be less than dims[1]. A python script is attached (read_gauges.py) that shows an example of this.

## Recommendations, restrictions
Some restriction, constantly updated as they are fixed:
 * Refrain from using OutputSimulation too often because (compared to the actual simulation) it may take a lot of time.
 * When using large h5 files, try compressing them, e.g. h5repack -i gaussian.h5 -o gaussian_compressed.h5 -f GZIP=9
