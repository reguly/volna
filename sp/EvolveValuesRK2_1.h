inline void EvolveValuesRK2_1(const float *dT, const float *inConservative, //OP_RW //temp
            const float *in, //OP_READ
            float *Conservative) //OP_WRITE
{
  Conservative[0] = inConservative[0] * *dT + in[0];
  Conservative[1] = inConservative[1] * *dT + in[1];
  Conservative[2] = inConservative[2] * *dT + in[2];

  //Conservative[0] = inConservative[0] + in[0];
  //Conservative[1] = inConservative[1] + in[1];
  //Conservative[2] = inConservative[2] + in[2];
  Conservative[3] = in[3];

  Conservative[0] = Conservative[0] <= EPS ? EPS : Conservative[0];
}
