//
// auto-generated by op2.py
//

//user function
__device__ void getMaxSpeed_gpu( const float* values, float* currentMaxSpeed) {

  if (values[0] > 1e-3f){
    float TruncatedH = values[0];
    float u = values[1]/TruncatedH;
    float v = values[2]/TruncatedH;
    float umax = currentMaxSpeed[1]/currentMaxSpeed[0];
    float vmax = currentMaxSpeed[2]/currentMaxSpeed[0];
    if (sqrt(u*u + v*v) > sqrt(umax*umax+vmax*vmax)) {
      currentMaxSpeed[0] = values[0];
      currentMaxSpeed[1] = values[1];
      currentMaxSpeed[2] = values[2];
      currentMaxSpeed[3] = values[3];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_getMaxSpeed(
  const float *__restrict arg0,
  float *arg1,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    getMaxSpeed_gpu(arg0+n*4,
                arg1+n*4);
  }
}


//host stub function
void op_par_loop_getMaxSpeed(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(19);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[19].name      = name;
  OP_kernels[19].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  getMaxSpeed");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_19
      int nthread = OP_BLOCK_SIZE_19;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_getMaxSpeed<<<nblocks,nthread>>>(
      (float *) arg0.data_d,
      (float *) arg1.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[19].time     += wall_t2 - wall_t1;
  OP_kernels[19].transfer += (float)set->size * arg0.size;
  OP_kernels[19].transfer += (float)set->size * arg1.size * 2.0f;
}
