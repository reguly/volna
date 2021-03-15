inline void Friction_manning(const float *dT,const float *M_n, //OP_RW, discard
            float *values) //OP_WRITE

{
  float F;
  float speed = sqrt(values[1]*values[1] + values[2]*values[2]);
  F = g* (*M_n * *M_n) *speed;
  float TruncatedH = values[0] < 10.f*EPS ? EPS : values[0];
  F = F/(pow(TruncatedH,4.0f/3.0f));
  values[0] = TruncatedH;
  values[1] = values[1] / (1.0f + F * *dT);
  values[2] = values[2] / (1.0f + F * *dT);
  values[3] = values[3];
}
