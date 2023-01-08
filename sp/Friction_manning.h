inline void Friction_manning(const float *dT,const float *M_n, //OP_RW, discard
            float *values) //OP_WRITE

{
  float Fr;
  float TruncatedH = values[0] < EPS ? EPS : values[0];
  float u = values[1]/TruncatedH;
  float v = values[2]/TruncatedH;
  float speed = sqrt(u*u + v*v);
  Fr = g* (*M_n * *M_n) *speed;
  Fr = Fr/(pow(TruncatedH,4.0f/3.0f));
  //float Fx = F*TruncatedH*values[1];
  //float Fy = F*TruncatedH*values[2];
  // Update Momentum
  values[0] = TruncatedH;
  values[3] = values[3];
  if (values[0] <= 1e-3f){
     values[1] = 0.0f;
     values[2] = 0.0f;
  } else if (1e-3f < values[0] <= 50.0f) {
     values[1] = values[1] / (1.0f + Fr * *dT);
     values[2] = values[2] / (1.0f + Fr * *dT);
  } else {
     values[1] = values[1];
     values[2] = values[2];
  }
}
