inline void EvolveValuesRK2_1(const float *dT, float *inConservative, //OP_RW //temp
            const float *in, //OP_READ
            float *midPointConservative, //OP_WRITE //temp
            float *midPoint) //OP_WRITE
{
  inConservative[0] *= *dT;
  inConservative[1] *= *dT;
  inConservative[2] *= *dT;
  //inConservative[3] *= *dT;
  //call to ToConservativeVariables inlined
  //inConservative[0] = in[0];
  //inConservative[1] = in[1];
  //inConservative[2] = in[2];
  //inConservative[3] = in[3];

  midPointConservative[0] = inConservative[0] + in[0];
  midPointConservative[1] = inConservative[1] + in[1];
  midPointConservative[2] = inConservative[2] + in[2];
  midPointConservative[3] = in[3];

  //call to ToPhysicalVariables inlined
  float TruncatedH = midPointConservative[0] < EPS ? EPS : midPointConservative[0];
  midPoint[0] = midPointConservative[0];
  midPoint[1] = midPointConservative[1];
  midPoint[2] = midPointConservative[2];
  midPoint[3] = midPointConservative[3];
}
