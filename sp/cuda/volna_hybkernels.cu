//
// auto-generated by op2.py
//

//header
#ifdef GPUPASS
#define op_par_loop_EvolveValuesRK3_1 op_par_loop_EvolveValuesRK3_1_gpu
#define op_par_loop_EvolveValuesRK3_2 op_par_loop_EvolveValuesRK3_2_gpu
#define op_par_loop_EvolveValuesRK3_3 op_par_loop_EvolveValuesRK3_3_gpu
#define op_par_loop_EvolveValuesRK3_4 op_par_loop_EvolveValuesRK3_4_gpu
#define op_par_loop_simulation_1 op_par_loop_simulation_1_gpu
#define op_par_loop_computeGradient op_par_loop_computeGradient_gpu
#define op_par_loop_limiter op_par_loop_limiter_gpu
#define op_par_loop_computeFluxes op_par_loop_computeFluxes_gpu
#define op_par_loop_NumericalFluxes op_par_loop_NumericalFluxes_gpu
#define op_par_loop_SpaceDiscretization op_par_loop_SpaceDiscretization_gpu
#define op_par_loop_incConst op_par_loop_incConst_gpu
#define op_par_loop_initEta_formula op_par_loop_initEta_formula_gpu
#define op_par_loop_initU_formula op_par_loop_initU_formula_gpu
#define op_par_loop_initV_formula op_par_loop_initV_formula_gpu
#define op_par_loop_values_operation2 op_par_loop_values_operation2_gpu
#define op_par_loop_applyConst op_par_loop_applyConst_gpu
#define op_par_loop_initBathymetry_large op_par_loop_initBathymetry_large_gpu
#define op_par_loop_initBathyRelative_formula op_par_loop_initBathyRelative_formula_gpu
#define op_par_loop_initBathymetry_formula op_par_loop_initBathymetry_formula_gpu
#define op_par_loop_initBathymetry_update op_par_loop_initBathymetry_update_gpu
#define op_par_loop_initBore_select op_par_loop_initBore_select_gpu
#define op_par_loop_initGaussianLandslide op_par_loop_initGaussianLandslide_gpu
#define op_par_loop_getTotalVol op_par_loop_getTotalVol_gpu
#define op_par_loop_getMaxElevation op_par_loop_getMaxElevation_gpu
#define op_par_loop_getMaxSpeed op_par_loop_getMaxSpeed_gpu
#define op_par_loop_gatherLocations op_par_loop_gatherLocations_gpu
#include "volna_kernels.cu"
#undef op_par_loop_EvolveValuesRK3_1
#undef op_par_loop_EvolveValuesRK3_2
#undef op_par_loop_EvolveValuesRK3_3
#undef op_par_loop_EvolveValuesRK3_4
#undef op_par_loop_simulation_1
#undef op_par_loop_computeGradient
#undef op_par_loop_limiter
#undef op_par_loop_computeFluxes
#undef op_par_loop_NumericalFluxes
#undef op_par_loop_SpaceDiscretization
#undef op_par_loop_incConst
#undef op_par_loop_initEta_formula
#undef op_par_loop_initU_formula
#undef op_par_loop_initV_formula
#undef op_par_loop_values_operation2
#undef op_par_loop_applyConst
#undef op_par_loop_initBathymetry_large
#undef op_par_loop_initBathyRelative_formula
#undef op_par_loop_initBathymetry_formula
#undef op_par_loop_initBathymetry_update
#undef op_par_loop_initBore_select
#undef op_par_loop_initGaussianLandslide
#undef op_par_loop_getTotalVol
#undef op_par_loop_getMaxElevation
#undef op_par_loop_getMaxSpeed
#undef op_par_loop_gatherLocations
#else
#define op_par_loop_EvolveValuesRK3_1 op_par_loop_EvolveValuesRK3_1_cpu
#define op_par_loop_EvolveValuesRK3_2 op_par_loop_EvolveValuesRK3_2_cpu
#define op_par_loop_EvolveValuesRK3_3 op_par_loop_EvolveValuesRK3_3_cpu
#define op_par_loop_EvolveValuesRK3_4 op_par_loop_EvolveValuesRK3_4_cpu
#define op_par_loop_simulation_1 op_par_loop_simulation_1_cpu
#define op_par_loop_computeGradient op_par_loop_computeGradient_cpu
#define op_par_loop_limiter op_par_loop_limiter_cpu
#define op_par_loop_computeFluxes op_par_loop_computeFluxes_cpu
#define op_par_loop_NumericalFluxes op_par_loop_NumericalFluxes_cpu
#define op_par_loop_SpaceDiscretization op_par_loop_SpaceDiscretization_cpu
#define op_par_loop_incConst op_par_loop_incConst_cpu
#define op_par_loop_initEta_formula op_par_loop_initEta_formula_cpu
#define op_par_loop_initU_formula op_par_loop_initU_formula_cpu
#define op_par_loop_initV_formula op_par_loop_initV_formula_cpu
#define op_par_loop_values_operation2 op_par_loop_values_operation2_cpu
#define op_par_loop_applyConst op_par_loop_applyConst_cpu
#define op_par_loop_initBathymetry_large op_par_loop_initBathymetry_large_cpu
#define op_par_loop_initBathyRelative_formula op_par_loop_initBathyRelative_formula_cpu
#define op_par_loop_initBathymetry_formula op_par_loop_initBathymetry_formula_cpu
#define op_par_loop_initBathymetry_update op_par_loop_initBathymetry_update_cpu
#define op_par_loop_initBore_select op_par_loop_initBore_select_cpu
#define op_par_loop_initGaussianLandslide op_par_loop_initGaussianLandslide_cpu
#define op_par_loop_getTotalVol op_par_loop_getTotalVol_cpu
#define op_par_loop_getMaxElevation op_par_loop_getMaxElevation_cpu
#define op_par_loop_getMaxSpeed op_par_loop_getMaxSpeed_cpu
#define op_par_loop_gatherLocations op_par_loop_gatherLocations_cpu
#include "../openmp/volna_kernels.cpp"
#undef op_par_loop_EvolveValuesRK3_1
#undef op_par_loop_EvolveValuesRK3_2
#undef op_par_loop_EvolveValuesRK3_3
#undef op_par_loop_EvolveValuesRK3_4
#undef op_par_loop_simulation_1
#undef op_par_loop_computeGradient
#undef op_par_loop_limiter
#undef op_par_loop_computeFluxes
#undef op_par_loop_NumericalFluxes
#undef op_par_loop_SpaceDiscretization
#undef op_par_loop_incConst
#undef op_par_loop_initEta_formula
#undef op_par_loop_initU_formula
#undef op_par_loop_initV_formula
#undef op_par_loop_values_operation2
#undef op_par_loop_applyConst
#undef op_par_loop_initBathymetry_large
#undef op_par_loop_initBathyRelative_formula
#undef op_par_loop_initBathymetry_formula
#undef op_par_loop_initBathymetry_update
#undef op_par_loop_initBore_select
#undef op_par_loop_initGaussianLandslide
#undef op_par_loop_getTotalVol
#undef op_par_loop_getMaxElevation
#undef op_par_loop_getMaxSpeed
#undef op_par_loop_gatherLocations

//user kernel files

void op_par_loop_EvolveValuesRK3_1_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_EvolveValuesRK3_1(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  if (OP_hybrid_gpu) {
    op_par_loop_EvolveValuesRK3_1_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

    }else{
    op_par_loop_EvolveValuesRK3_1_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

  }
}
#else
void op_par_loop_EvolveValuesRK3_1(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  op_par_loop_EvolveValuesRK3_1_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_EvolveValuesRK3_2_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_EvolveValuesRK3_2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_EvolveValuesRK3_2_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_EvolveValuesRK3_2_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_EvolveValuesRK3_2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_EvolveValuesRK3_2_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_EvolveValuesRK3_3_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_EvolveValuesRK3_3(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  if (OP_hybrid_gpu) {
    op_par_loop_EvolveValuesRK3_3_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

    }else{
    op_par_loop_EvolveValuesRK3_3_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

  }
}
#else
void op_par_loop_EvolveValuesRK3_3(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  op_par_loop_EvolveValuesRK3_3_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_EvolveValuesRK3_4_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_EvolveValuesRK3_4(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  if (OP_hybrid_gpu) {
    op_par_loop_EvolveValuesRK3_4_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3);

    }else{
    op_par_loop_EvolveValuesRK3_4_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3);

  }
}
#else
void op_par_loop_EvolveValuesRK3_4(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  op_par_loop_EvolveValuesRK3_4_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_simulation_1_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_simulation_1(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  if (OP_hybrid_gpu) {
    op_par_loop_simulation_1_gpu(name, set,
      arg0,
      arg1);

    }else{
    op_par_loop_simulation_1_cpu(name, set,
      arg0,
      arg1);

  }
}
#else
void op_par_loop_simulation_1(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  op_par_loop_simulation_1_gpu(name, set,
    arg0,
    arg1);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_computeGradient_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_computeGradient(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  if (OP_hybrid_gpu) {
    op_par_loop_computeGradient_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9);

    }else{
    op_par_loop_computeGradient_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9);

  }
}
#else
void op_par_loop_computeGradient(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  op_par_loop_computeGradient_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_limiter_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_limiter(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  if (OP_hybrid_gpu) {
    op_par_loop_limiter_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6);

    }else{
    op_par_loop_limiter_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6);

  }
}
#else
void op_par_loop_limiter(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  op_par_loop_limiter_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_computeFluxes_gpu(char const *name, op_set set,
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
  op_arg arg14);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_computeFluxes(char const *name, op_set set,
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
  op_arg arg14){

  if (OP_hybrid_gpu) {
    op_par_loop_computeFluxes_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9,
      arg10,
      arg11,
      arg12,
      arg13,
      arg14);

    }else{
    op_par_loop_computeFluxes_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9,
      arg10,
      arg11,
      arg12,
      arg13,
      arg14);

  }
}
#else
void op_par_loop_computeFluxes(char const *name, op_set set,
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
  op_arg arg14){

  op_par_loop_computeFluxes_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9,
    arg10,
    arg11,
    arg12,
    arg13,
    arg14);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_NumericalFluxes_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_NumericalFluxes(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  if (OP_hybrid_gpu) {
    op_par_loop_NumericalFluxes_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8);

    }else{
    op_par_loop_NumericalFluxes_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8);

  }
}
#else
void op_par_loop_NumericalFluxes(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  op_par_loop_NumericalFluxes_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_SpaceDiscretization_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_SpaceDiscretization(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  if (OP_hybrid_gpu) {
    op_par_loop_SpaceDiscretization_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9);

    }else{
    op_par_loop_SpaceDiscretization_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8,
      arg9);

  }
}
#else
void op_par_loop_SpaceDiscretization(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  op_par_loop_SpaceDiscretization_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8,
    arg9);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_incConst_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_incConst(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_incConst_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_incConst_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_incConst(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_incConst_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initEta_formula_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initEta_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_initEta_formula_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_initEta_formula_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_initEta_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_initEta_formula_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initU_formula_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initU_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_initU_formula_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_initU_formula_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_initU_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_initU_formula_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initV_formula_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initV_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_initV_formula_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_initV_formula_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_initV_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_initV_formula_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_values_operation2_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_values_operation2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  if (OP_hybrid_gpu) {
    op_par_loop_values_operation2_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

    }else{
    op_par_loop_values_operation2_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4);

  }
}
#else
void op_par_loop_values_operation2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  op_par_loop_values_operation2_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_applyConst_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_applyConst(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_applyConst_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_applyConst_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_applyConst(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_applyConst_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initBathymetry_large_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initBathymetry_large(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  if (OP_hybrid_gpu) {
    op_par_loop_initBathymetry_large_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7);

    }else{
    op_par_loop_initBathymetry_large_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7);

  }
}
#else
void op_par_loop_initBathymetry_large(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  op_par_loop_initBathymetry_large_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initBathyRelative_formula_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initBathyRelative_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  if (OP_hybrid_gpu) {
    op_par_loop_initBathyRelative_formula_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3);

    }else{
    op_par_loop_initBathyRelative_formula_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3);

  }
}
#else
void op_par_loop_initBathyRelative_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  op_par_loop_initBathyRelative_formula_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initBathymetry_formula_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initBathymetry_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_initBathymetry_formula_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_initBathymetry_formula_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_initBathymetry_formula(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_initBathymetry_formula_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initBathymetry_update_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initBathymetry_update(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  if (OP_hybrid_gpu) {
    op_par_loop_initBathymetry_update_gpu(name, set,
      arg0,
      arg1);

    }else{
    op_par_loop_initBathymetry_update_cpu(name, set,
      arg0,
      arg1);

  }
}
#else
void op_par_loop_initBathymetry_update(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  op_par_loop_initBathymetry_update_gpu(name, set,
    arg0,
    arg1);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initBore_select_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initBore_select(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  if (OP_hybrid_gpu) {
    op_par_loop_initBore_select_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8);

    }else{
    op_par_loop_initBore_select_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7,
      arg8);

  }
}
#else
void op_par_loop_initBore_select(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  op_par_loop_initBore_select_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7,
    arg8);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_initGaussianLandslide_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_initGaussianLandslide(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  if (OP_hybrid_gpu) {
    op_par_loop_initGaussianLandslide_gpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7);

    }else{
    op_par_loop_initGaussianLandslide_cpu(name, set,
      arg0,
      arg1,
      arg2,
      arg3,
      arg4,
      arg5,
      arg6,
      arg7);

  }
}
#else
void op_par_loop_initGaussianLandslide(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  op_par_loop_initGaussianLandslide_gpu(name, set,
    arg0,
    arg1,
    arg2,
    arg3,
    arg4,
    arg5,
    arg6,
    arg7);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_getTotalVol_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_getTotalVol(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  if (OP_hybrid_gpu) {
    op_par_loop_getTotalVol_gpu(name, set,
      arg0,
      arg1,
      arg2);

    }else{
    op_par_loop_getTotalVol_cpu(name, set,
      arg0,
      arg1,
      arg2);

  }
}
#else
void op_par_loop_getTotalVol(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  op_par_loop_getTotalVol_gpu(name, set,
    arg0,
    arg1,
    arg2);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_getMaxElevation_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_getMaxElevation(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  if (OP_hybrid_gpu) {
    op_par_loop_getMaxElevation_gpu(name, set,
      arg0,
      arg1);

    }else{
    op_par_loop_getMaxElevation_cpu(name, set,
      arg0,
      arg1);

  }
}
#else
void op_par_loop_getMaxElevation(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  op_par_loop_getMaxElevation_gpu(name, set,
    arg0,
    arg1);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_getMaxSpeed_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_getMaxSpeed(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  if (OP_hybrid_gpu) {
    op_par_loop_getMaxSpeed_gpu(name, set,
      arg0,
      arg1);

    }else{
    op_par_loop_getMaxSpeed_cpu(name, set,
      arg0,
      arg1);

  }
}
#else
void op_par_loop_getMaxSpeed(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  op_par_loop_getMaxSpeed_gpu(name, set,
    arg0,
    arg1);

  }
#endif //OP_HYBRID_GPU

void op_par_loop_gatherLocations_gpu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1);

//GPU host stub function
#if OP_HYBRID_GPU
void op_par_loop_gatherLocations(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  if (OP_hybrid_gpu) {
    op_par_loop_gatherLocations_gpu(name, set,
      arg0,
      arg1);

    }else{
    op_par_loop_gatherLocations_cpu(name, set,
      arg0,
      arg1);

  }
}
#else
void op_par_loop_gatherLocations(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  op_par_loop_gatherLocations_gpu(name, set,
    arg0,
    arg1);

  }
#endif //OP_HYBRID_GPU
#endif
