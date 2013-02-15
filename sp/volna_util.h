#include "op_lib_cpp.h"

//
// Prints the max. distance of adjacent cells
//
void printBandwidth(op_map map);

/*
 * Prints the adjacency list to a file
 */
void printAdjList(const char* filename, op_map map);

/*
 * Get permutation list based on specified reordering
 */
void op_get_permutation(int **perm, int **iperm, op_map map, op_set set);

void op_reorder_map(op_map map, int *perm, int *iperm, op_set set);

void op_reorder_dat(op_dat dat, int *iperm, op_set set);

void check_map(op_map map);

void is_all_referenced_once(int set_size, int* iperm);

void op_remove_nonref(int *cnode, int ncell, int *isRefNode, int *nnode, float **x);
