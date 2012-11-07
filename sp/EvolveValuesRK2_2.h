inline void EvolveValuesRK2_2(const float *dT, float *outConservative, //OP_RW, discard
            float *inConservative, //OP_READ, discard
            float *midPointConservative, //OP_RW, null
            float *out) //OP_WRITE

{
  float temp[4];
  temp[0] = 0.5*(outConservative[0] * *dT + midPointConservative[0] + inConservative[0]);
  temp[1] = 0.5*(outConservative[1] * *dT + midPointConservative[1] + inConservative[1]);
  temp[2] = 0.5*(outConservative[2] * *dT + midPointConservative[2] + inConservative[2]);

  temp[0] = temp[0] <= EPS ? EPS : temp[0];
  temp[3] = inConservative[3];

  //call to ToPhysicalVariables inlined
  float TruncatedH = temp[0] < EPS ? EPS : temp[0];
  out[0] = temp[0];
  out[1] = temp[1] / TruncatedH;
  out[2] = temp[2] / TruncatedH;
  out[3] = temp[3];
  
  midPointConservative[0] = 0.0f;
  midPointConservative[1] = 0.0f;
  midPointConservative[2] = 0.0f;
  midPointConservative[3] = 0.0f;
  outConservative[0] = 0.0f;
  outConservative[1] = 0.0f;
  outConservative[2] = 0.0f;
  outConservative[3] = 0.0f;
}
