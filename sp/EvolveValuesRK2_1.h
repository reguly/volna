inline void EvolveValuesRK2_1(const float *dT, const float *Lw_n, //OP_RW //temp
            const float *in, //OP_READ
            float *out) //OP_WRITE
{
  out[0] = Lw_n[0] * *dT + in[0];
  out[1] = Lw_n[1] * *dT + in[1];
  out[2] = Lw_n[2] * *dT + in[2];
  out[3] = in[3]-in[0];
  
  float TruncatedH = out[0] < EPS ? EPS : out[0];
  out[0] = TruncatedH;
  out[3] += TruncatedH;
}
