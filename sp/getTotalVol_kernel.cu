//
// auto-generated by op2.py on 2014-10-07 13:54
//

//user function
__device__
#include "getTotalVol.h"

// CUDA kernel function
__global__ void op_cuda_getTotalVol(
  float *arg0,
  float *arg1,
  float *arg2,
  int   offset_s,    
  int   set_size ) {
  
  float arg1_l[4];
  float arg2_l[1];
  for ( int d=0; d<1; d++ ){
    arg2_l[d]=ZERO_float;
  }
  int   tid = threadIdx.x%OP_WARPSIZE;
  
  extern __shared__ char shared[];
  char *arg_s = shared + offset_s*(threadIdx.x/OP_WARPSIZE);
  
  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n = n+=blockDim.x*gridDim.x ){
    int offset = n - tid;
    int nelems = MIN(OP_WARPSIZE,set_size-offset);
    //copy data into shared memory, then into local
    for ( int m=0; m<4; m++ ){
      ((float *)arg_s)[tid+m*nelems] = arg1[tid+m*nelems+offset*4];
    }
    
    for ( int m=0; m<4; m++ ){
      arg1_l[m] = ((float *)arg_s)[m+tid*4];
    }
    
    
    //user-supplied kernel call
    getTotalVol(arg0+n,
                arg1_l,
                arg2_l);
    //copy back into shared memory, then to device
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
  
  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  getTotalVol");
  }
  
  op_mpi_halo_exchanges_cuda(set, nargs, args);
  
  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timers_core(&cpu_t1, &wall_t1);
  
  if (set->size > 0) {
    
    op_timing_realloc(17);
    OP_kernels[17].name      = name;
    OP_kernels[17].count    += 1;
    
    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_17
      int nthread = OP_BLOCK_SIZE_17;
    #else
    //  int nthread = OP_block_size;
      int nthread = 128;
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
    
    //work out shared memory requirements per element
    
    int nshared = 0;
    nshared = MAX(nshared,sizeof(float)*4);
    
    //execute plan
    int offset_s = nshared*OP_WARPSIZE;
    
    nshared = MAX(nshared*nthread,reduct_size*nthread);
    op_cuda_getTotalVol<<<nblocks,nthread,nshared>>>(
      (float *) arg0.data_d,
      (float *) arg1.data_d,
      (float *) arg2.data_d,
      offset_s,
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
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[17].time     += wall_t2 - wall_t1;
  OP_kernels[17].transfer += (float)set->size * arg0.size;
  OP_kernels[17].transfer += (float)set->size * arg1.size;
}
