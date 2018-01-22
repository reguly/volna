/*Copyright 2018, Frederic Dias, Serge Guillas, Istvan Reguly

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "op_seq.h"

#include "volna_common.h"
#include "computeFluxes.h"
#include "NumericalFluxes.h"
#include "SpaceDiscretization.h"


void spaceDiscretization(op_dat data_in, op_dat data_out, float *minTimestep,
                         op_dat bathySource, op_dat edgeFluxes, op_dat maxEdgeEigenvalues,
                         op_dat edgeNormals, op_dat edgeLength, op_dat cellVolumes, op_dat isBoundary,
                         op_set cells, op_set edges, op_map edgesToCells, op_map cellsToEdges, int most) {
  {
    *minTimestep = INFINITY;
    { //Following loops merged:
      //FacetsValuesFromCellValues
      //FacetsValuesFromCellValues
      //spaceDiscretisation_1
      //NumericalFluxes_1
      //SpaceDiscretization
      op_par_loop(computeFluxes, "computeFluxes", edges,
                  op_arg_dat(data_in, 0, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(data_in, 1, edgesToCells, 4, "float", OP_READ),
                  op_arg_dat(edgeLength, -1, OP_ID, 1, "float", OP_READ),
                  op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                  op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(bathySource, -1, OP_ID, 2, "float", OP_WRITE),
                  op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_WRITE),
                  op_arg_dat(maxEdgeEigenvalues, -1, OP_ID, 1, "float", OP_WRITE));

    }
#ifdef DEBUG
    printf("maxFacetEigenvalues %g edgeLen %g cellVol %g\n", normcomp(maxEdgeEigenvalues, 0), normcomp(edgeLength, 0), normcomp(cellVolumes, 0));
#endif
    op_par_loop(NumericalFluxes, "NumericalFluxes", cells,
                op_arg_dat(maxEdgeEigenvalues, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(maxEdgeEigenvalues, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 0, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 1, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(edgeLength, 2, cellsToEdges, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, -1, OP_ID, 1, "float", OP_READ),
                op_arg_dat(data_out, -1, OP_ID, 4, "float", OP_WRITE),
                op_arg_gbl(minTimestep,1,"float", OP_MIN));

    //end NumericalFluxes
    op_par_loop(SpaceDiscretization, "SpaceDiscretization", edges,
                op_arg_dat(data_out, 0, edgesToCells, 4, "float", OP_INC), //again, Zb is not needed
                op_arg_dat(data_out, 1, edgesToCells, 4, "float", OP_INC),
                op_arg_dat(edgeFluxes, -1, OP_ID, 3, "float", OP_READ),
                op_arg_dat(bathySource, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(edgeNormals, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(isBoundary, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(cellVolumes, 0, edgesToCells, 1, "float", OP_READ),
                op_arg_dat(cellVolumes, 1, edgesToCells, 1, "float", OP_READ));
    }
}
