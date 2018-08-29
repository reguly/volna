inline void computeGradient(const float *center,
                            const float *neighbour1,
                            const float *neighbour2,
                            const float *neighbour3,
                            const float *centerVol,
                            const float *nb1Vol,
                            const float *nb2Vol,
                            const float *nb3Vol,
                            const float *edgeLength1,
                            const float *edgeLength2,
                            const float *edgeLength3,
                            const float *edgeNormal1,
                            const float *edgeNormal2,
                            const float *edgeNormal3,
                            float *q, float *out)
{
  float theta[3];

  theta[0] = (nb1Vol[0])/(centerVol[0] + nb1Vol[0]);
  theta[1] = (nb2Vol[0])/(centerVol[0] + nb2Vol[0]);
  theta[2] = (nb3Vol[0])/(centerVol[0] + nb3Vol[0]);

  out[0] = edgeLength1[0]*((1.0f - theta[0])*center[0] + (theta[0] *neighbour1[0]))*edgeNormal1[0];
  out[0] += edgeLength2[0]*((1.0f - theta[1])*center[0] + (theta[1] *neighbour2[0]))*edgeNormal2[0];
  out[0] += edgeLength3[0]*((1.0f - theta[2])*center[0] + (theta[2] *neighbour3[0]))*edgeNormal3[0];
  out[0] /= centerVol[0];
 
  out[1] = edgeLength1[0]*((1.0f - theta[0])*center[0] + (theta[0] *neighbour1[0]))*edgeNormal1[1];
  out[1] += edgeLength2[0]*((1.0f - theta[1])*center[0] + (theta[1] *neighbour2[0]))*edgeNormal2[1];
  out[1] += edgeLength3[0]*((1.0f - theta[2])*center[0] + (theta[2] *neighbour3[0]))*edgeNormal3[1];
  out[1] /= centerVol[0];

  out[2] = edgeLength1[0]*((1.0f - theta[0])*center[1] + (theta[0] *neighbour1[1]))*edgeNormal1[0];
  out[2] += edgeLength2[0]*((1.0f - theta[1])*center[1] + (theta[1] *neighbour2[1]))*edgeNormal2[0];
  out[2] += edgeLength3[0]*((1.0f - theta[2])*center[1] + (theta[2] *neighbour3[1]))*edgeNormal3[0];
  out[2] /= centerVol[0];

  out[3] = edgeLength1[0]*((1.0f - theta[0])*center[1] + (theta[0] *neighbour1[1]))*edgeNormal1[1];
  out[3] += edgeLength2[0]*((1.0f - theta[1])*center[1] + (theta[1] *neighbour2[1]))*edgeNormal2[1];
  out[3] += edgeLength3[0]*((1.0f - theta[2])*center[1] + (theta[2] *neighbour3[1]))*edgeNormal3[1];
  out[3] /= centerVol[0];

  out[4] = edgeLength1[0]*((1.0f - theta[0])*center[2] + (theta[0] *neighbour1[2]))*edgeNormal1[0];
  out[4] += edgeLength2[0]*((1.0f - theta[1])*center[2] + (theta[1] *neighbour2[2]))*edgeNormal2[0];
  out[4] += edgeLength3[0]*((1.0f - theta[2])*center[2] + (theta[2] *neighbour3[2]))*edgeNormal3[0];
  out[4] /= centerVol[0];

  out[5] = edgeLength1[0]*((1.0f - theta[0])*center[2] + (theta[0] *neighbour1[2]))*edgeNormal1[1];
  out[5] += edgeLength2[0]*((1.0f - theta[1])*center[2] + (theta[1] *neighbour2[2]))*edgeNormal2[1];
  out[5] += edgeLength3[0]*((1.0f - theta[2])*center[2] + (theta[2] *neighbour3[2]))*edgeNormal3[1];
  out[5] /= centerVol[0];

  out[6] = edgeLength1[0]*((1.0f - theta[0])*center[3] + (theta[0] *neighbour1[3]))*edgeNormal1[0];
  out[6] += edgeLength2[0]*((1.0f - theta[1])*center[3] + (theta[1] *neighbour2[3]))*edgeNormal2[0];
  out[6] += edgeLength3[0]*((1.0f - theta[2])*center[3] + (theta[2] *neighbour3[3]))*edgeNormal3[0];
  out[6] /= centerVol[0];

  out[7] = edgeLength1[0]*((1.0f - theta[0])*center[3] + (theta[0] *neighbour1[3]))*edgeNormal1[1];
  out[7] += edgeLength2[0]*((1.0f - theta[1])*center[3] + (theta[1] *neighbour2[3]))*edgeNormal2[1];
  out[7] += edgeLength3[0]*((1.0f - theta[2])*center[3] + (theta[2] *neighbour3[3]))*edgeNormal3[1];
  out[7] /= centerVol[0];
  
  q[0] = center[0] < neighbour1[0] ? center[0] : neighbour1[0];
  q[0] = q[0] < neighbour2[0] ? q[0] : neighbour2[0];
  q[0] = q[0] < neighbour3[0] ? q[0] : neighbour3[0];
  q[1] = center[0] > neighbour1[0] ? center[0] : neighbour1[0];
  q[1] = q[1] > neighbour2[0] ? q[1] : neighbour2[0];
  q[1] = q[1] > neighbour3[0] ? q[1] : neighbour3[0];

  q[2] = center[1] < neighbour1[1] ? center[1] : neighbour1[1];
  q[2] = q[2] < neighbour2[1] ? q[2] : neighbour2[1];
  q[2] = q[2] < neighbour3[1] ? q[2] : neighbour3[1];
  q[3] = center[1] > neighbour1[1] ? center[1] : neighbour1[1];
  q[3] = q[3] > neighbour2[1] ? q[3] : neighbour2[1];
  q[3] = q[3] > neighbour3[1] ? q[3] : neighbour3[1];

  q[4] = center[2] < neighbour1[2] ? center[2] : neighbour1[2];
  q[4] = q[4] < neighbour2[2] ? q[4] : neighbour2[2];
  q[4] = q[4] < neighbour3[2] ? q[4] : neighbour3[2];
  q[5] = center[2] > neighbour1[2] ? center[2] : neighbour1[2];
  q[5] = q[5] > neighbour2[2] ? q[5] : neighbour2[2];
  q[5] = q[5] > neighbour3[2] ? q[5] : neighbour3[2];

  q[6] = center[3] < neighbour1[3] ? center[3] : neighbour1[3];
  q[6] = q[6] < neighbour2[3] ? q[6] : neighbour2[3];
  q[6] = q[6] < neighbour3[3] ? q[6] : neighbour3[3];
  q[7] = center[3] > neighbour1[3] ? center[3] : neighbour1[3];
  q[7] = q[7] > neighbour2[3] ? q[7] : neighbour2[3];
  q[7] = q[7] > neighbour3[3] ? q[7] : neighbour3[3];

}
