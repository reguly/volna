//
// auto-generated by op2.py
//

//user function
inline void getTotalVol(const float* cellVolume, const float* value, float* totalVol) {
  (*totalVol) += (*cellVolume) * value[0];
}

// host stub function
void op_par_loop_getTotalVol(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  //create aligned pointers for dats
  ALIGNED_float const float * __restrict__ ptr0 = (float *) arg0.data;
  __assume_aligned(ptr0,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);
  float arg2h = *(float *)arg2.data;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  getTotalVol");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    #pragma omp parallel for reduction(+:arg2h)
    for ( int n=0; n<(set_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      float dat2[SIMD_VEC] = {0.0};
      #pragma omp simd aligned(ptr0,ptr1)
      for ( int i=0; i<SIMD_VEC; i++ ){
        getTotalVol(
          &(ptr0)[1 * (n+i)],
          &(ptr1)[4 * (n+i)],
          &dat2[i]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        arg2h += dat2[i];
      }
    }
    //remainder
    for ( int n=(set_size/SIMD_VEC)*SIMD_VEC; n<set_size; n++ ){
    #else
    #pragma omp parallel for simd aligned(ptr0,ptr1) reduction(+:arg2h)
    for ( int n=0; n<set_size; n++ ){
    #endif
      getTotalVol(
        &(ptr0)[1*n],
        &(ptr1)[4*n],
        &arg2h);
    }
  }

  // combine reduction data
  *(float*)arg2.data = arg2h;
  op_mpi_reduce(&arg2,(float*)arg2.data);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;
  OP_kernels[10].time     += wall_t2 - wall_t1;
  OP_kernels[10].transfer += (float)set->size * arg0.size;
  OP_kernels[10].transfer += (float)set->size * arg1.size;
}
#undef VECTORIZE
