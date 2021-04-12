#ifndef __VOLNA_COMMON_H
#define __VOLNA_COMMON_H

#include "op_seq.h"

void WriteMeshToVTKBinary(const char* filename, op_dat nodeCoords, int nnode, op_map cellsToNodes, int ncell, op_dat values, const float *zmin);
void WriteMeshToVTKAscii(const char* filename, op_dat nodeCoords, int nnode, op_map cellsToNodes, int ncell, op_dat values, const float *zmin);

#endif
