//
// auto-generated by op2.py
//

//user function
#include "../Timestep.h"

// host stub function
void op_par_loop_Timestep(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  float*arg7h = (float *)arg7.data;
  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 2;
  int  inds[8] = {0,0,0,1,1,1,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: Timestep\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_25
    int part_size = OP_PART_SIZE_25;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  // allocate and initialise arrays for global reduction
  float arg7_l[nthreads*64];
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg7_l[d+thr*64]=arg7h[d];
    }
  }

  if (set_size >0) {

    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
          int map0idx;
          int map1idx;
          int map2idx;
          map0idx = arg0.map_data[n * arg0.map->dim + 0];
          map1idx = arg0.map_data[n * arg0.map->dim + 1];
          map2idx = arg0.map_data[n * arg0.map->dim + 2];


          Timestep(
            &((float*)arg0.data)[1 * map0idx],
            &((float*)arg0.data)[1 * map1idx],
            &((float*)arg0.data)[1 * map2idx],
            &((float*)arg3.data)[1 * map0idx],
            &((float*)arg3.data)[1 * map1idx],
            &((float*)arg3.data)[1 * map2idx],
            &((float*)arg6.data)[1 * n],
            &arg7_l[64*omp_get_thread_num()]);
        }
      }

      // combine reduction data
      if (col == Plan->ncolors_owned-1) {
        for ( int thr=0; thr<nthreads; thr++ ){
          for ( int d=0; d<1; d++ ){
            arg7h[d]  = MIN(arg7h[d],arg7_l[d+thr*64]);
          }
        }
      }
      block_offset += nblocks;
    }
    OP_kernels[25].transfer  += Plan->transfer;
    OP_kernels[25].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_reduce(&arg7,arg7h);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].time     += wall_t2 - wall_t1;
}
