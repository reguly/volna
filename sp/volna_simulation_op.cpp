//
// auto-generated by op2.py
//

/*Copyright 2018, Frederic Dias, Serge Guillas, Istvan Reguly

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include  "op_lib_cpp.h"

//
// op_par_loop declarations
//
#ifdef OPENACC
#ifdef __cplusplus
extern "C" {
#endif
#endif

void op_par_loop_computeGradient(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_limiter(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_computeFluxes(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_NumericalFluxes1(char const *, op_set,
  op_arg );

void op_par_loop_SpaceDiscretization(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );

void op_par_loop_computeMinTimestep(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg );
#ifdef OPENACC
#ifdef __cplusplus
}
#endif
#endif


#include "volna_common.h"
#include "computeGradient.h"
#include "limiter.h"
#include "computeFluxes.h"
#include "NumericalFluxes1.h"
#include "SpaceDiscretization.h"
#include "computeMinTimestep.h"

#ifdef SLOPE
#include "executor.h"

void op_par_loop_SpaceDiscretization_slope(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  tile_t*);

void op_par_loop_computeGradient_slope(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  tile_t*);

void op_par_loop_limiter_slope(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  tile_t* tile);

void op_par_loop_computeFluxes_slope(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  tile_t* tile);

void op_par_loop_NumericalFluxes1_slope(char const *name, op_set set,
  op_arg arg0, tile_t* tile);

#endif

#ifdef SLOPE
void spaceDiscretization(op_dat data_in, op_dat data_out, float *minTimestep,
                         op_dat bathySource, op_dat edgeFluxes, op_dat maxEdgeEigenvalues,
                         op_dat edgeNormals, op_dat edgeLength, op_dat cellVolumes, op_dat isBoundary,
                         op_set cells, op_set edges, op_map edgesToCells, op_map cellsToEdges,
                         op_map cellsToCells, op_dat edgeCenters, op_dat cellCenters, op_dat GradientatCell, op_dat q, op_dat lim, int most,
                         executor_t* exec, int nColors) {


  //for each colour
  for (int color = 0; color < nColors; color++) {
  // for all tiles of this color
    const int nTilesPerColor = exec_tiles_per_color (exec, color);

    #pragma omp parallel for
    for (int j = 0; j < nTilesPerColor; j++) {
      // execute the tile
      tile_t* tile = exec_tile_at (exec, color, j);
      int tileLoopSize;

      // loop computeGradient
      op_par_loop_computeGradient_slope("computeGradient",cells,
                  op_arg_dat(data_in,-1,OP_ID,4,"float",OP_READ),
                  op_arg_dat(data_in,0,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(data_in,1,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(data_in,2,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(cellCenters,-1,OP_ID,2,"float",OP_READ),
                  op_arg_dat(cellCenters,0,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(cellCenters,1,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(cellCenters,2,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(q,-1,OP_ID,8,"float",OP_WRITE),
                  op_arg_dat(GradientatCell,-1,OP_ID,8,"float",OP_WRITE),
                  tile);

      /*iterations_list& lc2c_0 = tile_get_local_map (tile, 0, "sl_cellsToCells");
      iterations_list& iterations_0 = tile_get_iterations (tile, 0);
      tileLoopSize = tile_loop_size (tile, 0);

      #pragma omp simd
      for (int k = 0; k < tileLoopSize; k++) {
          computeGradient(
            (float*)(data_in->data + ((iterations_0[k] * 4) * sizeof(float))),
            (float*)(data_in->data + ((lc2c_0[k * N_NODESPERCELL + 0] * 4) * sizeof(float))),
            (float*)(data_in->data + ((lc2c_0[k * N_NODESPERCELL + 1] * 4) * sizeof(float))),
            (float*)(data_in->data + ((lc2c_0[k * N_NODESPERCELL + 2] * 4) * sizeof(float))),
            (float*)(cellCenters->data + ((iterations_0[k] * 2) * sizeof(float))),
            (float*)(cellCenters->data + ((lc2c_0[k * N_NODESPERCELL + 0] * 2) * sizeof(float))),
            (float*)(cellCenters->data + ((lc2c_0[k * N_NODESPERCELL + 1] * 2) * sizeof(float))),
            (float*)(cellCenters->data + ((lc2c_0[k * N_NODESPERCELL + 2] * 2) * sizeof(float))),
            (float*)(q->data + ((iterations_0[k] * 8) * sizeof(float))),
            (float*)(GradientatCell->data + ((iterations_0[k] * 8) * sizeof(float))));
      }*/


      // *minTimestep = INFINITY;
      // loop limiter
      op_par_loop_limiter_slope("limiter",cells,
                op_arg_dat(q,-1,OP_ID,8,"float",OP_READ),
                op_arg_dat(lim,-1,OP_ID,4,"float",OP_WRITE),
                op_arg_dat(data_in,-1,OP_ID,4,"float",OP_READ),
                op_arg_dat(GradientatCell,-1,OP_ID,8,"float",OP_READ),
                op_arg_dat(edgeCenters,0,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(edgeCenters,1,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(edgeCenters,2,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(cellCenters,-1,OP_ID,2,"float",OP_READ),
                tile);

      /*iterations_list& lc2e_1 = tile_get_local_map (tile, 1, "sl_cellsToEdges");
      iterations_list& iterations_1 = tile_get_iterations (tile, 1);
      tileLoopSize = tile_loop_size (tile, 1);

      #pragma omp simd
      for (int k = 0; k < tileLoopSize; k++) {
        limiter(
          (float*)(q->data + ((iterations_1[k] * 8) * sizeof(float))),
          (float*)(lim->data + ((iterations_1[k] * 4) * sizeof(float))),
          (float*)(data_in->data + ((iterations_1[k] * 4) * sizeof(float))),
          (float*)(GradientatCell->data + ((iterations_1[k] * 8) * sizeof(float))),
          (float*)(edgeCenters->data + ((lc2e_1[k * N_NODESPERCELL + 0] * 2) * sizeof(float))),
          (float*)(edgeCenters->data + ((lc2e_1[k * N_NODESPERCELL + 1] * 2) * sizeof(float))),
          (float*)(edgeCenters->data + ((lc2e_1[k * N_NODESPERCELL + 2] * 2) * sizeof(float))),
          (float*)(cellCenters->data + ((iterations_1[k] * 2) * sizeof(float)))
        );
      }*/

      // loop computeFluxes
      op_par_loop_computeFluxes_slope("computeFluxes",edges,
                op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(lim,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(lim,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(edgeLength,-1,OP_ID,1,"float",OP_READ),
                op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(cellCenters,0,edgesToCells,2,"float",OP_READ),
                op_arg_dat(cellCenters,1,edgesToCells,2,"float",OP_READ),
                op_arg_dat(edgeCenters,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(GradientatCell,0,edgesToCells,8,"float",OP_READ),
                op_arg_dat(GradientatCell,1,edgesToCells,8,"float",OP_READ),
                op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
                op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_WRITE),
                op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_WRITE),
                op_arg_dat(maxEdgeEigenvalues,-1,OP_ID,1,"float",OP_WRITE),
                tile);

      /*iterations_list& le2c_2 = tile_get_local_map (tile, 2, "sl_edgesToCells");
      iterations_list& iterations_2 = tile_get_iterations (tile, 2);
      tileLoopSize = tile_loop_size (tile, 2);

      #pragma omp simd
      for (int k = 0; k < tileLoopSize; k++) {
        computeFluxes(
          (float*)(data_in->data + ((le2c_2[k * N_CELLSPEREDGE + 0] * 4) * sizeof(float))),
          (float*)(data_in->data + ((le2c_2[k * N_CELLSPEREDGE + 1] * 4) * sizeof(float))),
          (float*)(lim->data + ((le2c_2[k * N_CELLSPEREDGE + 0] * 4) * sizeof(float))),
          (float*)(lim->data + ((le2c_2[k * N_CELLSPEREDGE + 1] * 4) * sizeof(float))),
          (float*)(edgeLength->data + ((iterations_2[k] * 1) * sizeof(float))),
          (float*)(edgeNormals->data + ((iterations_2[k] * 2) * sizeof(float))),
          (float*)(cellCenters->data + ((le2c_2[k * N_CELLSPEREDGE + 0] * 2) * sizeof(float))),
          (float*)(cellCenters->data + ((le2c_2[k * N_CELLSPEREDGE + 1] * 2) * sizeof(float))),
          (float*)(edgeCenters->data + ((iterations_2[k] * 2) * sizeof(float))),
          (float*)(GradientatCell->data + ((le2c_2[k * N_CELLSPEREDGE + 0] * 8) * sizeof(float))),
          (float*)(GradientatCell->data + ((le2c_2[k * N_CELLSPEREDGE + 1] * 8) * sizeof(float))),
          (int*)(isBoundary->data + ((iterations_2[k] * 1) * sizeof(int))),
          (float*)(bathySource->data + ((iterations_2[k] * 4) * sizeof(float))),
          (float*)(edgeFluxes->data + ((iterations_2[k] * 3) * sizeof(float))),
          (float*)(maxEdgeEigenvalues->data + ((iterations_2[k] * 1) * sizeof(float)))
        );
      }*/

      // loop NumericalFluxes
      op_par_loop_NumericalFluxes1_slope("NumericalFluxes1",cells,
              op_arg_dat(data_out,-1,OP_ID,4,"float",OP_WRITE),tile);

      /*iterations_list& iterations_3 = tile_get_iterations (tile, 3);
      tileLoopSize = tile_loop_size (tile, 3);
      #pragma omp simd
      for (int k = 0; k < tileLoopSize; k++) {
        NumericalFluxes1(
          (float*)(data_out->data + ((iterations_3[k] * 4) * sizeof(float)))
        );
      }*/

      // loop SpaceDiscretization
      op_par_loop_SpaceDiscretization_slope("SpaceDiscretization",edges,
                op_arg_dat(data_out,0,edgesToCells,4,"float",OP_INC),
                op_arg_dat(data_out,1,edgesToCells,4,"float",OP_INC),
                op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_READ),
                op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_READ),
                op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
                op_arg_dat(cellVolumes,0,edgesToCells,1,"float",OP_READ),
                op_arg_dat(cellVolumes,1,edgesToCells,1,"float",OP_READ),
                tile);

      /*iterations_list& le2c_4 = tile_get_local_map (tile, 4, "sl_edgesToCells");
      iterations_list& iterations_4 = tile_get_iterations (tile, 4);
      tileLoopSize = tile_loop_size (tile, 4);

      for (int k = 0; k < tileLoopSize; k++) {
        SpaceDiscretization(
          (float*)(data_out->data + ((le2c_4[k * N_CELLSPEREDGE + 0] * 4) * sizeof(float))),
          (float*)(data_out->data + ((le2c_4[k * N_CELLSPEREDGE + 1] * 4) * sizeof(float))),
          (float*)(data_in->data + ((le2c_4[k * N_CELLSPEREDGE + 0] * 4) * sizeof(float))),
          (float*)(data_in->data + ((le2c_4[k * N_CELLSPEREDGE + 1] * 4) * sizeof(float))),
          (float*)(edgeFluxes->data + ((iterations_4[k] * 3) * sizeof(float))),
          (float*)(bathySource->data + ((iterations_4[k] * 4) * sizeof(float))),
          (float*)(edgeNormals->data + ((iterations_4[k] * 2) * sizeof(float))),
          (int*)(isBoundary->data + ((iterations_4[k] * 1) * sizeof(int))),
          (float*)(cellVolumes->data + ((le2c_4[k * N_CELLSPEREDGE + 0] * 1) * sizeof(float))),
          (float*)(cellVolumes->data + ((le2c_4[k * N_CELLSPEREDGE + 1] * 1) * sizeof(float)))
        );
      }*/

    }

  }

  //  op_par_loop_limiter("limiter",cells,
  //               op_arg_dat(q,-1,OP_ID,8,"float",OP_READ),
  //               op_arg_dat(lim,-1,OP_ID,4,"float",OP_WRITE),
  //               op_arg_dat(data_in,-1,OP_ID,4,"float",OP_READ),
  //               op_arg_dat(GradientatCell,-1,OP_ID,8,"float",OP_READ),
  //               op_arg_dat(edgeCenters,0,cellsToEdges,2,"float",OP_READ),
  //               op_arg_dat(edgeCenters,1,cellsToEdges,2,"float",OP_READ),
  //               op_arg_dat(edgeCenters,2,cellsToEdges,2,"float",OP_READ),
  //               op_arg_dat(cellCenters,-1,OP_ID,2,"float",OP_READ));


    

    // op_par_loop_computeFluxes("computeFluxes",edges,
    //             op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(lim,0,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(lim,1,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(edgeLength,-1,OP_ID,1,"float",OP_READ),
    //             op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
    //             op_arg_dat(cellCenters,0,edgesToCells,2,"float",OP_READ),
    //             op_arg_dat(cellCenters,1,edgesToCells,2,"float",OP_READ),
    //             op_arg_dat(edgeCenters,-1,OP_ID,2,"float",OP_READ),
    //             op_arg_dat(GradientatCell,0,edgesToCells,8,"float",OP_READ),
    //             op_arg_dat(GradientatCell,1,edgesToCells,8,"float",OP_READ),
    //             op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
    //             op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_WRITE),
    //             op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_WRITE),
    //             op_arg_dat(maxEdgeEigenvalues,-1,OP_ID,1,"float",OP_WRITE));

    

    // op_par_loop_NumericalFluxes1("NumericalFluxes1",cells,
    //             op_arg_dat(data_out,-1,OP_ID,4,"float",OP_WRITE));
    // //end NumericalFluxes
    // op_par_loop_SpaceDiscretization("SpaceDiscretization",edges,
    //             op_arg_dat(data_out,0,edgesToCells,4,"float",OP_INC),
    //             op_arg_dat(data_out,1,edgesToCells,4,"float",OP_INC),
    //             op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
    //             op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_READ),
    //             op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_READ),
    //             op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
    //             op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
    //             op_arg_dat(cellVolumes,0,edgesToCells,1,"float",OP_READ),
    //             op_arg_dat(cellVolumes,1,edgesToCells,1,"float",OP_READ));
    


}
#else
void spaceDiscretization(op_dat data_in, op_dat data_out, float *minTimestep,
                         op_dat bathySource, op_dat edgeFluxes, op_dat maxEdgeEigenvalues,
                         op_dat edgeNormals, op_dat edgeLength, op_dat cellVolumes, op_dat isBoundary,
                         op_set cells, op_set edges, op_map edgesToCells, op_map cellsToEdges,
                         op_map cellsToCells, op_dat edgeCenters, op_dat cellCenters, op_dat GradientatCell, op_dat q, op_dat lim, int most) {
  {
    { op_par_loop_computeGradient("computeGradient",cells,
                  op_arg_dat(data_in,-1,OP_ID,4,"float",OP_READ),
                  op_arg_dat(data_in,0,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(data_in,1,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(data_in,2,cellsToCells,4,"float",OP_READ),
                  op_arg_dat(cellCenters,-1,OP_ID,2,"float",OP_READ),
                  op_arg_dat(cellCenters,0,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(cellCenters,1,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(cellCenters,2,cellsToCells,2,"float",OP_READ),
                  op_arg_dat(q,-1,OP_ID,8,"float",OP_WRITE),
                  op_arg_dat(GradientatCell,-1,OP_ID,8,"float",OP_WRITE));

    }
    // *minTimestep = INFINITY;
    op_par_loop_limiter("limiter",cells,
                op_arg_dat(q,-1,OP_ID,8,"float",OP_READ),
                op_arg_dat(lim,-1,OP_ID,4,"float",OP_WRITE),
                op_arg_dat(data_in,-1,OP_ID,4,"float",OP_READ),
                op_arg_dat(GradientatCell,-1,OP_ID,8,"float",OP_READ),
                op_arg_dat(edgeCenters,0,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(edgeCenters,1,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(edgeCenters,2,cellsToEdges,2,"float",OP_READ),
                op_arg_dat(cellCenters,-1,OP_ID,2,"float",OP_READ));


    {

    op_par_loop_computeFluxes("computeFluxes",edges,
                op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(lim,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(lim,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(edgeLength,-1,OP_ID,1,"float",OP_READ),
                op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(cellCenters,0,edgesToCells,2,"float",OP_READ),
                op_arg_dat(cellCenters,1,edgesToCells,2,"float",OP_READ),
                op_arg_dat(edgeCenters,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(GradientatCell,0,edgesToCells,8,"float",OP_READ),
                op_arg_dat(GradientatCell,1,edgesToCells,8,"float",OP_READ),
                op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
                op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_WRITE),
                op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_WRITE),
                op_arg_dat(maxEdgeEigenvalues,-1,OP_ID,1,"float",OP_WRITE));

    }

    op_par_loop_NumericalFluxes1("NumericalFluxes1",cells,
                op_arg_dat(data_out,-1,OP_ID,4,"float",OP_WRITE));
    //end NumericalFluxes
    op_par_loop_SpaceDiscretization("SpaceDiscretization",edges,
                op_arg_dat(data_out,0,edgesToCells,4,"float",OP_INC),
                op_arg_dat(data_out,1,edgesToCells,4,"float",OP_INC),
                op_arg_dat(data_in,0,edgesToCells,4,"float",OP_READ),
                op_arg_dat(data_in,1,edgesToCells,4,"float",OP_READ),
                op_arg_dat(edgeFluxes,-1,OP_ID,3,"float",OP_READ),
                op_arg_dat(bathySource,-1,OP_ID,4,"float",OP_READ),
                op_arg_dat(edgeNormals,-1,OP_ID,2,"float",OP_READ),
                op_arg_dat(isBoundary,-1,OP_ID,1,"int",OP_READ),
                op_arg_dat(cellVolumes,0,edgesToCells,1,"float",OP_READ),
                op_arg_dat(cellVolumes,1,edgesToCells,1,"float",OP_READ));
    }
    // op_par_loop_computeMinTimestep("computeMinTimestep",cells,
    //             op_arg_dat(maxEdgeEigenvalues,0,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(maxEdgeEigenvalues,1,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(maxEdgeEigenvalues,2,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(edgeLength,0,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(edgeLength,1,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(edgeLength,2,cellsToEdges,1,"float",OP_READ),
    //             op_arg_dat(cellVolumes,-1,OP_ID,1,"float",OP_READ),
    //             op_arg_gbl(minTimestep,1,"float",OP_MIN));

    // printf("minTimeStep=%f\n", *minTimestep);

}
#endif