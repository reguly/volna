//
// auto-generated by op2.py
//

// header
#include "op_lib_cpp.h"
#define double_ALIGN 128
#define float_ALIGN 64
#define int_ALIGN 64
#define VECTORIZE
#ifdef VECTORIZE
#define SIMD_VEC 4
#define ALIGNED_double __attribute__((aligned(double_ALIGN)))
#define ALIGNED_float __attribute__((aligned(float_ALIGN)))
#define ALIGNED_int __attribute__((aligned(int_ALIGN)))
#else
#define ALIGNED_double
#define ALIGNED_float
#define ALIGNED_int
#endif
#undef VECTORIZE

// global constants
extern float CFL;
extern float EPS;
extern float g;
// user kernel files
#include "EvolveValuesRK2_1_ompveckernel.cpp"
#include "EvolveValuesRK2_2_ompveckernel.cpp"
#include "simulation_1_ompveckernel.cpp"
#include "computeFluxes_ompveckernel.cpp"
#include "NumericalFluxes_ompveckernel.cpp"
#include "SpaceDiscretization_ompveckernel.cpp"
#include "getTotalVol_ompveckernel.cpp"
#include "getMaxElevation_ompveckernel.cpp"
#include "getMaxSpeed_ompveckernel.cpp"
#include "gatherLocations_ompveckernel.cpp"
#include "incConst_ompveckernel.cpp"
#include "initEta_formula_ompveckernel.cpp"
#include "initU_formula_ompveckernel.cpp"
#include "initV_formula_ompveckernel.cpp"
#include "values_operation2_ompveckernel.cpp"
#include "applyConst_ompveckernel.cpp"
#include "initBathymetry_large_ompveckernel.cpp"
#include "initBathyRelative_formula_ompveckernel.cpp"
#include "initBathymetry_formula_ompveckernel.cpp"
#include "initBathymetry_update_ompveckernel.cpp"
#include "initBore_select_ompveckernel.cpp"
#include "initGaussianLandslide_ompveckernel.cpp"
