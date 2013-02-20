#include "volna_writeVTK.h"
#include <stdio.h>
#include "volna_common.h"

/*
 * Utility function for binary output: swaps byte endianneses
 */
inline float swapEndiannesDouble(double d) {
  union {
    double d;
    char b[8];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.d;
}

/*
 * Utility function for binary output: swaps byte endianneses
 */
inline float swapEndiannesFloat(float f) {
  union {
    float f;
    char b[4];
  } dat1, dat2;
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

/*
 * Utility function for binary output: swaps byte endianneses
 */
inline int swapEndiannesInt(int d) {
  union {
    int d;
    char b[4];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.d;
}

/*
 * Write simulation output to binary file
 */
void WriteMeshToVTKBinary(const char* filename, op_dat nodeCoords, int nnode, op_map cellsToNodes, int ncell, op_dat values) {
  op_printf("Writing OutputSimulation to binary file: %s \n",filename);
  FILE* fp;
  fp = fopen(filename, "w");
  if(fp == NULL) {
    op_printf("can't open file for write %s\n",filename);
    exit(-1);
  }

  // write header
  char s[256];
  strcpy(s, "# vtk DataFile Version 2.0\n Output from OP2 Volna.\n"); fwrite(s, sizeof(char), strlen(s), fp);
  strcpy(s, "BINARY \nDATASET UNSTRUCTURED_GRID\n\n"); fwrite(s, sizeof(char), strlen(s), fp);

  // write vertices
  sprintf(s,"POINTS %d float\n", nnode); fwrite(s, sizeof(char), strlen(s), fp);

 float* nodeCoords_data;
  nodeCoords_data = (float*)nodeCoords->data;
  float tmp_float;
  int i = 0;
  for (i = 0; i < nnode; ++i) {
    tmp_float = swapEndiannesFloat(nodeCoords_data[i*MESH_DIM  ]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
    tmp_float = swapEndiannesFloat(nodeCoords_data[i*MESH_DIM+1]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
    tmp_float = swapEndiannesFloat(0.0);
    fwrite(&tmp_float, sizeof(float), 1, fp);
  }
  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

  // write cells
  sprintf(s, "CELLS %d %d\n", ncell, 4*ncell); fwrite(s, sizeof(char), strlen(s), fp);

  int three = 3;
  int tmp_int;
  for ( i = 0; i < ncell; ++i ) {
    tmp_int = swapEndiannesInt(three);
    fwrite(&tmp_int, sizeof(int), 1, fp);
    tmp_int = swapEndiannesInt(cellsToNodes->map[i*N_NODESPERCELL  ]);
    fwrite(&tmp_int, sizeof(int), 1, fp);
    tmp_int = swapEndiannesInt(cellsToNodes->map[i*N_NODESPERCELL+1]);
    fwrite(&tmp_int, sizeof(int), 1, fp);
    tmp_int = swapEndiannesInt(cellsToNodes->map[i*N_NODESPERCELL+2]);
    fwrite(&tmp_int, sizeof(int), 1, fp);
  }
  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

  // write cell types (5 for triangles)
  sprintf(s, "CELL_TYPES %d\n", ncell); fwrite(s, sizeof(char), strlen(s), fp);

  int five=5;
  for ( i=0; i<ncell; ++i ) {
    tmp_int = swapEndiannesInt(five);
    fwrite(&tmp_int, sizeof(int), 1, fp);
  }

  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);
  float* values_data;
  values_data = (float*) values->data;
    sprintf(s, "CELL_DATA %d\n"
                  "SCALARS Eta float 1\n"
                  "LOOKUP_TABLE default\n",
                  ncell); fwrite(s, sizeof(char), strlen(s), fp);

  for ( i=0; i<ncell; ++i ) {
    tmp_float = swapEndiannesFloat(values_data[i*N_STATEVAR] + values_data[i*N_STATEVAR+3]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
  }

    strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

    strcpy(s, "SCALARS U float 1\nLOOKUP_TABLE default\n"); fwrite(s, sizeof(char), strlen(s), fp);
  for ( i=0; i<ncell; ++i ){
    tmp_float = swapEndiannesFloat(values_data[i*N_STATEVAR+1]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
  }
  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);


  strcpy(s, "SCALARS V float 1\nLOOKUP_TABLE default\n"); fwrite(s, sizeof(char), strlen(s), fp);
  for ( i=0; i<ncell; ++i ) {
    tmp_float = swapEndiannesFloat(values_data[i*N_STATEVAR+2]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
  }

  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

  strcpy(s, "SCALARS Bathymetry float 1\nLOOKUP_TABLE default\n"); fwrite(s, sizeof(char), strlen(s), fp);
  for ( i=0; i<ncell; ++i ) {
    tmp_float = swapEndiannesFloat(values_data[i*N_STATEVAR+3]);
    fwrite(&tmp_float, sizeof(float), 1, fp);
  }

  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

  strcpy(s, "SCALARS Visual float 1\nLOOKUP_TABLE default\n"); fwrite(s, sizeof(char), strlen(s), fp);
  float hundred = 100.0;
  for ( i=0; i<ncell; ++i ) {
    if(values_data[i*N_STATEVAR] < 1e-3){
      tmp_float = swapEndiannesFloat(hundred);
      fwrite(&tmp_float, sizeof(float), 1, fp);
    }
    else {
      tmp_float = swapEndiannesFloat(values_data[i*N_STATEVAR] + values_data[i*N_STATEVAR+3]);
      fwrite(&tmp_float, sizeof(float), 1, fp);
    }
  }
  strcpy(s, "\n"); fwrite(s, sizeof(char), strlen(s), fp);

  if(fclose(fp) != 0) {
    op_printf("can't close file %s\n",filename);
    exit(-1);
  }
}

/*
 * Write simulation output to ASCII file
 */
void WriteMeshToVTKAscii(const char* filename, op_dat nodeCoords, int nnode, op_map cellsToNodes, int ncell, op_dat values) {
  op_printf("Writing OutputSimulation to ASCII file: %s \n",filename);
  FILE* fp;
  fp = fopen(filename, "w");
  if(fp == NULL) {
    op_printf("can't open file for write %s\n",filename);
    exit(-1);
  }

  // write header
  fprintf(fp,"# vtk DataFile Version 2.0\n Output from OP2 Volna.\n");
  fprintf(fp,"ASCII \nDATASET UNSTRUCTURED_GRID\n\n");
  // write vertices
  fprintf(fp,"POINTS %d float\n", nnode);
  float* nodeCoords_data;
  nodeCoords_data = (float*)nodeCoords->data;
  int i = 0;
  for (i = 0; i < nnode; ++i) {
    fprintf(fp, "%g %g %g \n",
        (float)nodeCoords_data[i*MESH_DIM  ],
        (float)nodeCoords_data[i*MESH_DIM+1],
        0.0);
  }
  fprintf(fp, "\n");
  fprintf(fp, "CELLS %d %d\n", ncell, 4*ncell);
  for ( i = 0; i < ncell; ++i ) {
    fprintf(fp, "3 %d %d %d\n",
        cellsToNodes->map[i*N_NODESPERCELL  ],
        cellsToNodes->map[i*N_NODESPERCELL+1],
        cellsToNodes->map[i*N_NODESPERCELL+2]);
  }
  fprintf(fp, "\n");
  // write cell types (5 for triangles)
  fprintf(fp, "CELL_TYPES %d\n", ncell);
  for ( i=0; i<ncell; ++i )
    fprintf(fp, "5\n");
  fprintf(fp, "\n");

  float* values_data;
  values_data = (float*) values->data;
  fprintf(fp, "CELL_DATA %d\n"
              "SCALARS Eta float 1\n"
              "LOOKUP_TABLE default\n",
              ncell);
  float tmp = 0.0;
  for ( i=0; i<ncell; ++i ) {
    tmp = values_data[i*N_STATEVAR] + values_data[i*N_STATEVAR+3];
    fprintf(fp, "%10.20g\n", values_data[i*N_STATEVAR] + values_data[i*N_STATEVAR+3]);
  }

  fprintf(fp, "\n");
  fprintf(fp, "SCALARS U float 1\n"
              "LOOKUP_TABLE default\n");
  for ( i=0; i<ncell; ++i )
    fprintf(fp, "%10.20g\n", values_data[i*N_STATEVAR+1]);
  fprintf(fp, "\n");
  fprintf(fp, "SCALARS V float 1\n"
              "LOOKUP_TABLE default\n");
  for ( i=0; i<ncell; ++i )
    fprintf(fp, "%10.20g\n", values_data[i*N_STATEVAR+2]);
  fprintf(fp, "\n");
  fprintf(fp, "SCALARS Bathymetry float 1\n"
              "LOOKUP_TABLE default\n");
  for ( i=0; i<ncell; ++i )
    fprintf(fp, "%10.20g\n", values_data[i*N_STATEVAR+3]);
  fprintf(fp, "\n");

  fprintf(fp, "SCALARS Visual float 1\n"
              "LOOKUP_TABLE default\n");

  for ( i=0; i<ncell; ++i ) {
    if(values_data[i*N_STATEVAR] < 1e-3)
      fprintf(fp, "%10.20g\n", 100.0);
    else
      fprintf(fp, "%10.20g\n", values_data[i*N_STATEVAR] + values_data[i*N_STATEVAR+3]);
  }
  fprintf(fp, "\n");

  if(fclose(fp) != 0) {
    op_printf("can't close file %s\n",filename);
    exit(-1);
  }
}
