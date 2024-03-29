/*Copyright 2018, Frederic Dias, Serge Guillas, Istvan Reguly

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "op_seq.h"

#include "volna_common.h"
#include "computeGradient.h"
#include "limiter.h"
#include "computeFluxes.h"
#include "Timestep.h"
#include "NumericalFluxes.h"
#include "computeFluxes_sph.h"
#include "NumericalFluxes_sph.h"


void spaceDiscretization(op_dat data_in, op_dat data_out, float *minTimestep,
                         op_dat bathySource, op_dat edgeFluxes, op_dat maxEdgeEigenvalues,
                         op_dat edgeNormals, op_dat edgeLength, op_dat cellVolumes, op_dat isBoundary,
                         op_set cells, op_set edges, op_map edgesToCells, op_map cellsToEdges,
                         op_map cellsToCells, op_dat edgeCenters, op_dat cellCenters, op_dat GradientatCell, op_dat q, op_dat lim, float *zmin) {
  {

    {
    // TO DO: Pre calculate the geometric mesh quantities
    op_par_loop(computeGradient, "computeGradient", cells,
                  op_arg_dat(data_in, -1, OP_ID, 4, "float", OP_READ),
                  op_arg_dat(data_in , 0, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in , 1, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in , 2, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(cellCenters, -1, OP_ID , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 0, cellsToCells , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 1, cellsToCells , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 2, cellsToCells , 2, "float", OP_RW),
                  op_arg_dat(q, -1, OP_ID, 8, "float", OP_WRITE),
                  op_arg_dat(GradientatCell, -1, OP_ID, 8, "float", OP_WRITE));
    }
   op_par_loop(limiter, "limiter", cells,
                op_arg_dat(q, -1, OP_ID, 8, "float", OP_READ),
                op_arg_dat(lim, -1, OP_ID, 4, "float", OP_WRITE),
                op_arg_dat(data_in, -1, OP_ID, 4, "float", OP_READ),
                op_arg_dat(GradientatCell, -1, OP_ID, 8, "float", OP_READ),
                op_arg_dat(edgeCenters, 0, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(edgeCenters, 1, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(edgeCenters, 2, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(data_out, -1, OP_ID, 4, "float", OP_WRITE),
                op_arg_dat(cellCenters, -1, OP_ID , 2, "float", OP_READ));

    {
    op_par_loop(computeFluxes, "computeFluxes", edges,
                  op_arg_dat(data_in, 0, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in, 1, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(lim, 0, edgesToCells,  4, "float", OP_READ),
                  op_arg_dat(lim, 1, edgesToCells,  4, "float", OP_READ),
                  op_arg_dat(edgeLength, -1, OP_ID, 1, "float", OP_READ),
                  op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 0, edgesToCells, 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 1, edgesToCells, 2, "float", OP_READ),
                  op_arg_dat(edgeCenters, -1, OP_ID, 2, "float", OP_READ),
                  op_arg_dat(GradientatCell, 0, edgesToCells, 8, "float", OP_READ),
                  op_arg_dat(GradientatCell, 1, edgesToCells, 8, "float", OP_READ),
                  op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(bathySource, -1, OP_ID, 4, "float", OP_WRITE),
                  op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_WRITE),
                  op_arg_dat(maxEdgeEigenvalues, -1, OP_ID, 1, "float", OP_WRITE),
                  op_arg_gbl(zmin, 1,"float", OP_READ));

    }

    if (*minTimestep >= 0.0){
    op_par_loop(Timestep, "Timestep", cells,
                op_arg_dat(maxEdgeEigenvalues, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, -1, OP_ID, 1, "float", OP_READ),
               op_arg_gbl(minTimestep,1,"float", OP_MIN));
    }

    op_par_loop(NumericalFluxes, "NumericalFluxes", edges,
                op_arg_dat(data_out, 0, edgesToCells, 4, "float", OP_INC),
                op_arg_dat(data_out, 1, edgesToCells, 4, "float", OP_INC),
                op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_READ),
                op_arg_dat(bathySource, -1, OP_ID, 4, "float", OP_READ),
                op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(cellVolumes, 0, edgesToCells, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, 1, edgesToCells, 1, "float", OP_READ));
    }
}

void spaceDiscretization_sph(op_dat data_in, op_dat data_out, float *minTimestep,
                         op_dat bathySource, op_dat edgeFluxes, op_dat maxEdgeEigenvalues,
                         op_dat edgeNormals, op_dat edgeLength, op_dat cellVolumes, op_dat isBoundary,
                         op_set cells, op_set edges, op_map edgesToCells, op_map cellsToEdges,
                         op_map cellsToCells, op_dat edgeCenters, op_dat cellCenters, op_dat GradientatCell, op_dat q, op_dat lim, float *zmin) {
  {

        {
    // TO DO: Pre calculate the geometric mesh quantities
    op_par_loop(computeGradient, "computeGradient", cells,
                  op_arg_dat(data_in, -1, OP_ID, 4, "float", OP_READ),
                  op_arg_dat(data_in , 0, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in , 1, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in , 2, cellsToCells, 4, "float", OP_READ),
                  op_arg_dat(cellCenters, -1, OP_ID , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 0, cellsToCells , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 1, cellsToCells , 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 2, cellsToCells , 2, "float", OP_RW),
                  op_arg_dat(q, -1, OP_ID, 8, "float", OP_WRITE),
                  op_arg_dat(GradientatCell, -1, OP_ID, 8, "float", OP_WRITE));
    }
   op_par_loop(limiter, "limiter", cells,
                op_arg_dat(q, -1, OP_ID, 8, "float", OP_READ),
                op_arg_dat(lim, -1, OP_ID, 4, "float", OP_WRITE),
                op_arg_dat(data_in, -1, OP_ID, 4, "float", OP_READ),
                op_arg_dat(GradientatCell, -1, OP_ID, 8, "float", OP_READ),
                op_arg_dat(edgeCenters, 0, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(edgeCenters, 1, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(edgeCenters, 2, cellsToEdges, 2, "float", OP_READ),
                op_arg_dat(data_out, -1, OP_ID, 4, "float", OP_WRITE),
                op_arg_dat(cellCenters, -1, OP_ID , 2, "float", OP_READ));

    {
    op_par_loop(computeFluxes_sph, "computeFluxes_sph", edges,
                  op_arg_dat(data_in, 0, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in, 1, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(lim, 0, edgesToCells,  4, "float", OP_READ),
                  op_arg_dat(lim, 1, edgesToCells,  4, "float", OP_READ),
                  op_arg_dat(edgeLength, -1, OP_ID, 1, "float", OP_READ),
                  op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 0, edgesToCells, 2, "float", OP_READ),
                  op_arg_dat(cellCenters, 1, edgesToCells, 2, "float", OP_READ),
                  op_arg_dat(edgeCenters, -1, OP_ID, 2, "float", OP_READ),
                  op_arg_dat(GradientatCell, 0, edgesToCells, 8, "float", OP_READ),
                  op_arg_dat(GradientatCell, 1, edgesToCells, 8, "float", OP_READ),
                  op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(bathySource, -1, OP_ID, 4, "float", OP_WRITE),
                  op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_WRITE),
                  op_arg_dat(maxEdgeEigenvalues, -1, OP_ID, 1, "float", OP_WRITE),
                  op_arg_gbl(zmin, 1,"float", OP_READ));

    }

    if (*minTimestep >= 0.0){
    op_par_loop(Timestep, "Timestep", cells,
                op_arg_dat(maxEdgeEigenvalues, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, -1, OP_ID, 1, "float", OP_READ),
                op_arg_gbl(minTimestep,1,"float", OP_MIN));
    }

    op_par_loop(NumericalFluxes_sph, "NumericalFluxes_sph", edges,
                op_arg_dat(data_out, 0, edgesToCells, 4, "float", OP_INC),
                op_arg_dat(data_out, 1, edgesToCells, 4, "float", OP_INC),
                op_arg_dat(cellCenters, 0, edgesToCells, 2, "float", OP_READ),
                op_arg_dat(cellCenters, 1, edgesToCells, 2, "float", OP_READ),
                op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_READ),
                op_arg_dat(bathySource, -1, OP_ID, 4, "float", OP_READ),
                op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(cellVolumes, 0, edgesToCells, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, 1, edgesToCells, 1, "float", OP_READ));
    }
}
