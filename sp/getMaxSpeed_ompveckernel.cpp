//
// auto-generated by op2.py
//

//user function
inline void getMaxSpeed(const float* values, float* currentMaxSpeed) {
  /*float tmp = values[0]+values[3];
  *currentMaxSpeed = *currentMaxSpeed > tmp ? *currentMaxSpeed : tmp;*/
  if (sqrt(values[1]*values[1]+values[2]*values[2]) > sqrt(currentMaxSpeed[1]*currentMaxSpeed[1]+currentMaxSpeed[2]*currentMaxSpeed[2])) {
    currentMaxSpeed[0] = values[0];
    currentMaxSpeed[1] = values[1];
    currentMaxSpeed[2] = values[2];
    currentMaxSpeed[3] = values[3];
  }
}

// host stub function
void op_par_loop_getMaxSpeed(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_float const float * __restrict__ ptr0 = (float *) arg0.data;
  __assume_aligned(ptr0,float_ALIGN);
  ALIGNED_float       float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  getMaxSpeed");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    #pragma omp parallel for
    for ( int n=0; n<(set_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      #pragma omp simd aligned(ptr0,ptr1)
      for ( int i=0; i<SIMD_VEC; i++ ){
        getMaxSpeed(
          &(ptr0)[4 * (n+i)],
          &(ptr1)[4 * (n+i)]);
      }
    }
    //remainder
    for ( int n=(set_size/SIMD_VEC)*SIMD_VEC; n<set_size; n++ ){
    #else
    #pragma omp parallel for simd aligned(ptr0,ptr1)
    for ( int n=0; n<set_size; n++ ){
    #endif
      getMaxSpeed(
        &(ptr0)[4*n],
        &(ptr1)[4*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg0.size;
  OP_kernels[12].transfer += (float)set->size * arg1.size * 2.0f;
}
#undef VECTORIZE
