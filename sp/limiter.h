inline void limiter(const float *q, float *lim,
                    const float *value, const float *gradient,
                    const float *edgecenter1, const float *edgecenter2,
                    const float *edgecenter3,
                    float *zeroInit,
                    const float *cellcenter)
{

  float facevalue[3], dx[3], dy[3];
  int i, j;
  float max[3], edgealpha[3];
  dx[0] = (edgecenter1[0] - cellcenter[0]);
  dy[0] = (edgecenter1[1] - cellcenter[1]);
  dx[1] = (edgecenter2[0] - cellcenter[0]);
  dy[1] = (edgecenter2[1] - cellcenter[1]);
  dx[2] = (edgecenter3[0] - cellcenter[0]);
  dy[2] = (edgecenter3[1] - cellcenter[1]);
  // If the cell is not on the wet/dry boundary
  if(q[0] > EPS){
  // The limiter is calculated for all physical variables using the
  // Barth-Jesperson formula and then the minimum limiter is used.
  // q[0] - Hmin , q[1] - Hmax
  // q[2] - Umin , q[3] - Umax
  // q[4] - Vmin , q[5] - Vmax
  // q[6] - Zmin , q[7] - Zmax
  for(j=0;j<4;j++){
   for(i =0 ; i<3; i++){
    facevalue[i] = value[j] + (((gradient[2*j]*dx[i]) + (gradient[2*j + 1]*dy[i])));
     if(facevalue[i] > q[2*j + 1]) {
      edgealpha[i] = (q[2*j + 1] - value[j]) / (facevalue[i] - value[j]);
     } else if (facevalue[i] < q[2*j]){
      edgealpha[i] = (q[2*j] - value[j]) / (facevalue[i] - value[j]);
     } else{
      edgealpha[i] = 1.0f;
     }
    max[i] = edgealpha[i] < 1.0f ? edgealpha[i] : 1.0f;
   }
   lim[j] = max[0] < max[1] ? max[0] : max[1];
   lim[j] = lim[j] < max[2] ? lim[j]: max[2];
  }
  lim[0] = lim[0] < lim[1] ? lim[0]: lim[1];
  lim[0] = lim[0] < lim[2] ? lim[0]: lim[2];
  lim[0] = lim[0] < lim[3] ? lim[0]: lim[3];
  } else {
    lim[0] = 0.0f;
    lim[1] = 0.0f;
    lim[2] = 0.0f;
    lim[3] = 0.0f;
  }
  zeroInit[0] = 0.0f;
  zeroInit[1] = 0.0f;
  zeroInit[2] = 0.0f;
  zeroInit[3] = 0.0f;
}
