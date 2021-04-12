/*Copyright 2018, Frederic Dias, Serge Guillas, Istvan Reguly

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "volna_common.h"
#include "applyConst.h"
#include "incConst.h"
#include "initBathymetry_formula.h"
#include "initBathymetry_large.h"
#include "initBathymetry_update.h"
#include "initBathyRelative_formula.h"
#include "initBore_select.h"
#include "initEta_formula.h"
#include "initGaussianLandslide.h"
#include "initU_formula.h"
#include "initV_formula.h"
#include "values_operation2.h"
#include "zero_bathy.h"

#include "op_seq.h"

void InitEta(op_set cells, op_dat cellCenters, op_dat values, op_dat initValues, int fromFile) {
#ifdef DEBUG
  op_printf("InitEta...");
#endif
  if (fromFile) {
    //overwrite values.H with values stored in initValues
    int variable = 1; //bitmask 1 - H, 2 - U, 4 - V, 8 - Zb
    //TODO: we are only overwriting H, moving the whole thing
    op_par_loop(incConst, "incConst", cells,
                op_arg_dat(initValues, -1, OP_ID, 1, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_gbl(&variable, 1, "int", OP_READ));
  } else {
    //TODO: document the fact that this actually adds to the value of V
    // i.e. user should only access values[2]
    op_par_loop(initEta_formula, "initEta_formula", cells,
                op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_INC),
                op_arg_gbl(&timestamp, 1, "double", OP_READ));
  }
#ifdef DEBUG
  op_printf("done\n");
#endif
}

void InitU(op_set cells, op_dat cellCenters, op_dat values, op_dat initValues, int fromFile) {
  //TODO: document the fact that this actually adds to the value of U
  // i.e. user should only access values[1]
  if(fromFile) {
  int variable = 2;
  op_par_loop(incConst, "incConst", cells,
                op_arg_dat(initValues, -1, OP_ID, 1, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_gbl(&variable, 1, "int", OP_READ));
  } else {

#ifdef DEBUG
  op_printf("InitU...");
#endif
  op_par_loop(initU_formula, "initU_formula", cells,
              op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
              op_arg_dat(values, -1, OP_ID, 4, "float", OP_INC),
              op_arg_gbl(&timestamp, 1, "double", OP_READ));
  }
#ifdef DEBUG
  op_printf("done\n");
#endif
}

  
void InitV(op_set cells, op_dat cellCenters, op_dat values, op_dat initValues, int fromFile) {
  //TODO: document the fact that this actually adds to the value of V
  // i.e. user should only access values[2]
#ifdef DEBUG
  op_printf("InitV...");
#endif
 if(fromFile) {
  int variable = 4;
  op_par_loop(incConst, "incConst", cells,
                op_arg_dat(initValues, -1, OP_ID, 1, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_gbl(&variable, 1, "int", OP_READ));

 } else { 
  op_par_loop(initV_formula, "initV_formula", cells,
              op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
              op_arg_dat(values, -1, OP_ID, 4, "float", OP_INC),
              op_arg_gbl(&timestamp, 1, "double", OP_READ));
#ifdef DEBUG
  op_printf("done\n");
#endif
}
}
//void OutputSimulation(op_set points, op_set cells, op_dat p_x, op_dat values) {
//
//
//}

void InitBathymetry(op_set cells, op_dat cellCenters, op_dat values, op_dat initValues, int fromFile, int firstTime, op_set bathy_nodes, op_set lifted_cells, op_map liftedcellsToBathyNodes, op_map liftedcellsToCells, op_dat bathy_xy, op_dat initial_zb, float *zmin) {
  if (firstTime) {
    int result = 0;
    int leftOperand = 0;
    int rightOperand = 3;
    int operation = 0; //0 +, 1 -, 2 *, 3 /
    op_par_loop(values_operation2, "values_operation2", cells,
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_gbl(&result, 1, "int", OP_READ),
                op_arg_gbl(&leftOperand, 1, "int", OP_READ),
                op_arg_gbl(&rightOperand, 1, "int", OP_READ),
                op_arg_gbl(&operation, 1, "int", OP_READ));
  }
  if (fromFile) {
    //overwrite values.Zb with values stored in initValues
    if (new_format) {
      int variable = 8; //bitmask 1 - H, 2 - U, 4 - V, 8 - Zb
      op_par_loop(applyConst, "applyConst", cells,
                  op_arg_dat(initial_zb, -1, OP_ID, 1, "float", OP_READ),
                  op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                  op_arg_gbl(&variable, 1, "int", OP_READ));
      op_par_loop(initBathymetry_large, "initBathymetry_large", lifted_cells,
                  op_arg_dat(values, 0, liftedcellsToCells, 4, "float", OP_INC),
                  op_arg_dat(cellCenters, 0, liftedcellsToCells, 2, "float", OP_READ),
                  op_arg_dat(bathy_xy, 0, liftedcellsToBathyNodes, 2, "float", OP_READ),
                  op_arg_dat(bathy_xy, 1, liftedcellsToBathyNodes, 2, "float", OP_READ),
                  op_arg_dat(bathy_xy, 2, liftedcellsToBathyNodes, 2, "float", OP_READ),
                  op_arg_dat(initValues, 0, liftedcellsToBathyNodes, 1, "float", OP_READ),
                  op_arg_dat(initValues, 1, liftedcellsToBathyNodes, 1, "float", OP_READ),
                  op_arg_dat(initValues, 2, liftedcellsToBathyNodes, 1, "float", OP_READ));
    } else {
      int variable = 8; //bitmask 1 - H, 2 - U, 4 - V, 8 - Zb
      //TODO: we are only overwriting H, moving the whole thing
      op_par_loop(applyConst, "applyConst", cells,
                  op_arg_dat(initValues, -1, OP_ID, 1, "float", OP_READ),
                  op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                  op_arg_gbl(&variable, 1, "int", OP_READ));
    }
  } else {
    //TODO: document the fact that this actually sets to the value of Zb
    // i.e. user should only access values[3]
    if (initial_zb != NULL) {
      op_par_loop(initBathyRelative_formula, "initBathyRelative_formula", cells,
                op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_dat(initial_zb, -1, OP_ID, 1, "float", OP_READ),
                op_arg_gbl(&timestamp, 1, "double", OP_READ));
    } else {
      op_par_loop(initBathymetry_formula, "initBathymetry_formula", cells,
                op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_INC),
                op_arg_gbl(&timestamp, 1, "double", OP_READ));
    }
  }
  if (firstTime) {
    *zmin = 0.0f;
    op_par_loop(zero_bathy, "zero_bathy", cells,
             op_arg_dat(values, -1, OP_ID, 4, "float", OP_READ),
             op_arg_gbl(zmin, 1, "float", OP_MIN));
    printf("zmin %f \n", *zmin);
  }
  op_par_loop(initBathymetry_update, "initBathymetry_update", cells,
              op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
              op_arg_gbl(zmin, 1, "float", OP_READ),
              op_arg_gbl(&firstTime, 1, "int", OP_READ));
#ifdef DEBUG
  printf("InitBathymetry executing H: %g Zb: %g\n", normcomp(values, 0), normcomp(values, 3));
#endif
}

void InitBore(op_set cells, op_dat cellCenters, op_dat values, BoreParams params) {
#ifdef DEBUG
  op_printf("InitBore...");
#endif
  float g = 9.81;
  float Fl = params.ul / sqrt( g * params.Hl );
  float Fs = params.S / sqrt( g * params.Hl );

  float r = .5 * ( sqrt( 1.0 + 8.0*( Fl - Fs )*(Fl - Fs ) ) - 1.0 );

  float Hr = r * params.Hl;
  float ur = params.S + ( params.ul - params.S ) / r;
  float vr = params.vl;
  ur *= -1.0;

  op_par_loop(initBore_select, "initBore_select", cells,
              op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
              op_arg_dat(cellCenters, -1, OP_ID, 2, "float", OP_READ),
              op_arg_gbl(&params.x0, 1, "float", OP_READ),
              op_arg_gbl(&params.Hl, 1, "float", OP_READ),
              op_arg_gbl(&params.ul, 1, "float", OP_READ),
              op_arg_gbl(&params.vl, 1, "float", OP_READ),
              op_arg_gbl(&Hr, 1, "float", OP_READ),
              op_arg_gbl(&ur, 1, "float", OP_READ),
              op_arg_gbl(&vr, 1, "float", OP_READ));
#ifdef DEBUG
  op_printf("done\n");
#endif
}

void InitGaussianLandslide(op_set cells, op_dat cellCenters, op_dat values, GaussianLandslideParams params, int firstTime) {
  //again, we only need Zb
  op_par_loop(initGaussianLandslide, "initGaussianLandslide", cells,
              op_arg_dat(cellCenters, -1, OP_ID, 2, "float",OP_READ),
              op_arg_dat(values, -1, OP_ID, 4, "float",OP_RW),
              op_arg_gbl(&params.mesh_xmin, 1, "float", OP_READ),
              op_arg_gbl(&params.A, 1, "float", OP_READ),
              op_arg_gbl(&timestamp, 1, "double", OP_READ),
              op_arg_gbl(&params.lx, 1, "float", OP_READ),
              op_arg_gbl(&params.ly, 1, "float", OP_READ),
              op_arg_gbl(&params.v, 1, "float", OP_READ));

  if (firstTime) {
    int result = 0;
    int leftOperand = 0;
    int rightOperand = 3;
    int operation = 1; //0 +, 1 -, 2 *, 3 /
    op_par_loop(values_operation2, "values_operation2", cells,
                op_arg_dat(values, -1, OP_ID, 4, "float", OP_RW),
                op_arg_gbl(&result, 1, "int", OP_READ),
                op_arg_gbl(&leftOperand, 1, "int", OP_READ),
                op_arg_gbl(&rightOperand, 1, "int", OP_READ),
                op_arg_gbl(&operation, 1, "int", OP_READ));
  }
}
