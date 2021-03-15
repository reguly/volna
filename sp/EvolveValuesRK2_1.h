inline void EvolveValuesRK2_1(const float *dT, const float *inConservative, //OP_RW //temp
            const float *in, //OP_READ
            float *out) //OP_WRITE
{
  out[0] = inConservative[0] * *dT + in[0];
  out[1] = inConservative[1] * *dT + in[1]*in[0];
  out[2] = inConservative[2] * *dT + in[2]*in[0];
  out[3] = in[3];

  float TruncatedH = out[0] < EPS ? EPS : out[0];
  out[1] = out[1] / TruncatedH;
  out[2] = out[2] / TruncatedH;

}
