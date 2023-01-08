inline void NumericalFluxes_sph(float *left, //OP_INC
              float *right, //OP_INC
              const float *leftcellCenters, const float *rightcellCenters,
              const float *edgeFluxes, //OP_READ
              const float *bathySource, //OP_READ
              const float *edgeNormals, const int *isRightBoundary,
              const float *cellVolumes0, //OP_READ
              const float *cellVolumes1) //OP_READ)
{
  left[0] -= (edgeFluxes[0])/(Radius*cellVolumes0[0]);
  left[1] -= (edgeFluxes[1] + bathySource[0] * edgeNormals[0])/(Radius*cellVolumes0[0]*cos(3.14159265358979323846f*leftcellCenters[1]/180.0f));
  left[2] -= (edgeFluxes[2] + bathySource[0] * edgeNormals[1])/(Radius*cellVolumes0[0]);
  // Added centered Source Term
  left[1] += (bathySource[2] *edgeNormals[0])/(Radius*cellVolumes0[0]*cos(3.14159265358979323846f*leftcellCenters[1]/180.0f));
  left[2] += (bathySource[2] *edgeNormals[1])/(Radius*cellVolumes0[0]);
  if (!*isRightBoundary) {
    right[0] += edgeFluxes[0]/cellVolumes1[0];
    right[1] += (edgeFluxes[1] + bathySource[1] * edgeNormals[0])/(Radius*cellVolumes1[0]*cos(3.14159265358979323846f*rightcellCenters[1]/180.0f));
    right[2] += (edgeFluxes[2] + bathySource[1] * edgeNormals[1])/(Radius*cellVolumes1[0]);
    // Added centered Source Term
    right[1] -= (bathySource[3] *edgeNormals[0])/(Radius*cellVolumes1[0]*cos(3.14159265358979323846f*rightcellCenters[1]/180.0f));
    right[2] -= (bathySource[3] *edgeNormals[1])/(Radius*cellVolumes1[0]);
   } else {
   right[0] += 0.0f;
   right[1] += 0.0f;
   right[2] += 0.0f;
  }

}
