inline void EvolveValuesRK2_2(const float *dT, float *outConservative, //OP_RW, discard
            const float *midPointConservative, //OP_READ, discard
            const float *in, //OP_READ, discard
            float *out) //OP_WRITE

{
  out[0] = 0.5*(outConservative[0] * *dT + midPointConservative[0] + in[0]);
  out[1] = 0.5*(outConservative[1] * *dT + midPointConservative[1] + in[1]);
  out[2] = 0.5*(outConservative[2] * *dT + midPointConservative[2] + in[2]);

  out[0] = out[0] <= EPS ? EPS : out[0];
  out[3] = in[3];

  //call to ToPhysicalVariables inlined
  float TruncatedH = outConservative[0] < EPS ? EPS : outConservative[0];
  //out[0] = outConservative[0];
  //out[1] = outConservative[1];
  //out[2] = outConservative[2];
  //out[3] = outConservative[3];
}
