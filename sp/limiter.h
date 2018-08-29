inline void limiter(const float *q, float *q2,
                    const float *value, const float *gradient,
                    const float *edgecenter1, const float *edgecenter2,
                    const float *edgecenter3, const float *cellcenter)
{

  float facevalue[3], dx[3], dy[3], y;
  int i, j;
  float max[3], edgealpha[3], left[3],right[3];

  dx[0] = (edgecenter1[0] - cellcenter[0]);
  dy[0] = (edgecenter1[1] - cellcenter[1]);
  dx[1] = (edgecenter2[0] - cellcenter[0]);
  dy[1] = (edgecenter2[1] - cellcenter[1]);
  dx[2] = (edgecenter3[0] - cellcenter[0]);
  dy[2] = (edgecenter3[1] - cellcenter[1]);
  // If the cell is not on the wet/dry boundary
  if((value[0] > 2.0f*EPS) && (q[0]> 2.0f*EPS)){
  // The limiter is calculated for all physical variables using the
  // Barth-Jesperson formula and then the minimum limiter is used.
  // q[0] - Hmin , q[1] - Hmax , q[8] - Halpha
  // q[2] - Umin , q[3] - Umax , q[9] - Ualpha
  // q[4] - Vmin , q[5] - Vmax , q[10] - Valpha
  // q[6] - Zmin , q[7] - Zmax , q[11] - Zalpha
  for(j=0;j<4;j++){
   for(i =0 ; i<3; i++){
    facevalue[i] = value[j] + ((gradient[2*j]*dx[i]) + (gradient[2*j + 1]*dy[i]));
     if(facevalue[i] > value[j]) {
      edgealpha[i] = (q[2*j + 1] - value[j]) / (facevalue[i] - value[j]);
     } else if (facevalue[i] < value[j]){
      edgealpha[i] = (q[2*j] - value[j]) / (facevalue[i] - value[j]);
     } else{
      edgealpha[i] = 1.0f;
     }
    // Alter the multiplication factor to obtain other limiters. 
    left[i] = 1.0f*edgealpha[i] < 1.0f ? 1.0f*edgealpha[i] : 1.0f;
    right[i] = edgealpha[i] < 1.0f ? edgealpha[i] : 1.0f;
    max[i] = left[i] > right[i] ? left[i] : right[i] ;
   } 
   q2[j] = max[0] < max[1] ? max[0] : max[1];
   q2[j] = q2[j] < max[2] ? q2[j]: max[2];
  }
  q2[0] = q2[0] < q2[1] ? q2[0]: q2[1];
  q2[0] = q2[0] < q2[2] ? q2[0]: q2[2];
  q2[0] = q2[0] < q2[3] ? q2[0]: q2[3];
  
  } else {
    q2[0] = 0.0f;
    q2[1] = 0.0f;
    q2[2] = 0.0f;
    q2[3] = 0.0f;
  }
}
