//
// auto-generated by op2.py on 2014-10-07 13:54
//

//user function
#include "initBathymetry_large.h"

// host stub function
void op_par_loop_initBathymetry_large(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

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
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 2;
  int  inds[8] = {-1,-1,0,0,0,1,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: initBathymetry_large\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_9
    int part_size = OP_PART_SIZE_9;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set->size >0) {

    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);

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
          int map2idx = arg2.map_data[n * arg2.map->dim + 0];
          int map3idx = arg2.map_data[n * arg2.map->dim + 1];
          int map4idx = arg2.map_data[n * arg2.map->dim + 2];

          initBathymetry_large(
            &((float*)arg0.data)[4 * n],
            &((float*)arg1.data)[2 * n],
            &((float*)arg2.data)[2 * map2idx],
            &((float*)arg2.data)[2 * map3idx],
            &((float*)arg2.data)[2 * map4idx],
            &((float*)arg5.data)[1 * map2idx],
            &((float*)arg5.data)[1 * map3idx],
            &((float*)arg5.data)[1 * map4idx]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[9].transfer  += Plan->transfer;
    OP_kernels[9].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;
  OP_kernels[9].time     += wall_t2 - wall_t1;
}
