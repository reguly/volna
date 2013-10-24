inline void SpaceDiscretization(float *left, //OP_INC
              float *right, //OP_INC
              float *edgeFluxes, //OP_READ
              float *bathySource, //OP_READ
              float *edgeNormals, int *isRightBoundary,
	      float *cellVolumes0, //OP_READ
	      float *cellVolumes1 //OP_READ
)
{
  left[0] -= (edgeFluxes[0])/cellVolumes0[0];
  left[1] -= (edgeFluxes[1] + bathySource[0] * edgeNormals[0])/cellVolumes0[0];
  left[2] -= (edgeFluxes[2] + bathySource[0] * edgeNormals[1])/cellVolumes0[0];

  right[0] += select(*isRightBoundary==0,edgeFluxes[0]/cellVolumes1[0],0.0f);
  right[1] += select(*isRightBoundary==0,(edgeFluxes[1] + bathySource[1] * edgeNormals[0])/cellVolumes1[0],0.0f);
  right[2] += select(*isRightBoundary==0,(edgeFluxes[2] + bathySource[1] * edgeNormals[1])/cellVolumes1[0],0.0f);
}