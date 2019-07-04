//
// auto-generated by op2.py
//

//user function
__device__ void getTotalVol_gpu( const float* cellVolume, const float* value, float* totalVol) {
  (*totalVol) += (*cellVolume) * value[0];
}

// CUDA kernel function
__global__ void op_cuda_getTotalVol(
  const float *__restrict arg0,
  const float *__restrict arg1,
  float *arg2,
  int   set_size ) {

  float arg2_l[1];
  for ( int d=0; d<1; d++ ){
    arg2_l[d]=ZERO_float;
  }

  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    getTotalVol_gpu(arg0+n*1,
                arg1+n*4,
                arg2_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg2[d+blockIdx.x*1],arg2_l[d]);
  }
}


//host stub function
void op_par_loop_getTotalVol(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  float*arg2h = (float *)arg2.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(17);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[17].name      = name;
  OP_kernels[17].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  getTotalVol");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_17
      int nthread = OP_BLOCK_SIZE_17;
    #else
      int nthread = OP_block_size;
    //  int nthread = 128;
    #endif

    int nblocks = 200;

    //transfer global reduction data to GPU
    int maxblocks = nblocks;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(float));
    reduct_size   = MAX(reduct_size,sizeof(float));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg2.data   = OP_reduct_h + reduct_bytes;
    arg2.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((float *)arg2.data)[d+b*1] = ZERO_float;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(float));
    mvReductArraysToDevice(reduct_bytes);

    int nshared = reduct_size*nthread;
    op_cuda_getTotalVol<<<nblocks,nthread,nshared>>>(
      (float *) arg0.data_d,
      (float *) arg1.data_d,
      (float *) arg2.data_d,
      set->size );
    //transfer global reduction data back to CPU
    mvReductArraysToHost(reduct_bytes);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg2h[d] = arg2h[d] + ((float *)arg2.data)[d+b*1];
      }
    }
    arg2.data = (char *)arg2h;
    op_mpi_reduce(&arg2,arg2h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[17].time     += wall_t2 - wall_t1;
  OP_kernels[17].transfer += (float)set->size * arg0.size;
  OP_kernels[17].transfer += (float)set->size * arg1.size;
}