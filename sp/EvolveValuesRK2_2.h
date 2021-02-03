inline void EvolveValuesRK2_2(const float *dT,const float *Lw_1, //OP_RW, discard
            const float *values, //OP_READ, discard
            const float *Conservative, //OP_READ, discard
            float *out) //OP_WRITE

{
  out[0] = 0.5*(Lw_1[0] * *dT + Conservative[0] + values[0]);
  out[1] = 0.5*(Lw_1[1] * *dT + Conservative[1] + values[1]);
  out[2] = 0.5*(Lw_1[2] * *dT + Conservative[2] + values[2]);

  out[0] = out[0] > EPS ? out[0] : EPS;
  out[3] = values[3];

}
