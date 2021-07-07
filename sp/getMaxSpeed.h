inline void getMaxSpeed(const float* values, float* currentMaxSpeed) {
  /*float tmp = values[0]+values[3];
  *currentMaxSpeed = *currentMaxSpeed > tmp ? *currentMaxSpeed : tmp;*/
  if (values[0] > 1e-3){
    float TruncatedH = values[0];// < EPS ? EPS : values[0];
    float u = values[1]/TruncatedH;
    float v = values[2]/TruncatedH;
    float umax = currentMaxSpeed[1]/currentMaxSpeed[0];
    float vmax = currentMaxSpeed[2]/currentMaxSpeed[0];
    if (sqrt(u*u + v*v) > sqrt(umax*umax+vmax*vmax)) {
      currentMaxSpeed[0] = values[0];
      currentMaxSpeed[1] = values[1];
      currentMaxSpeed[2] = values[2];
      currentMaxSpeed[3] = values[3];
    }
  }
}
