inline void Friction_manning(const float *dT,const float *M_n, //OP_RW, discard
            float *values) //OP_WRITE

{
  float F;
  float TruncatedH = values[0] < EPS ? EPS : values[0];
  float u = values[1]/TruncatedH;
  float v = values[2]/TruncatedH;
  float speed = sqrt(u*u + v*v);
  //float speed = sqrt(values[1]*values[1] + values[2]*values[2]);
  F = g* (*M_n * *M_n) *speed;
  F = F/(pow(TruncatedH,4.0f/3.0f));
  //float Fx = F*TruncatedH*values[1];
  //float Fy = F*TruncatedH*values[2];
  // Update Momentum 
  if (values[0] <= EPS){
     values[0] = TruncatedH;
     values[1] = 0.0f;
     values[2] = 0.0f;
     values[3] = values[3];
  } else {  
  values[0] = TruncatedH;
  values[1] = values[1] / (1.0f + F * *dT);
  values[2] = values[2] / (1.0f + F * *dT);
  values[3] = values[3];
  }
}
