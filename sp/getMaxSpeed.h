inline void getMaxSpeed(const float* values, float* currentMaxSpeed) {
  /*float tmp = values[0]+values[3];
  *currentMaxSpeed = *currentMaxSpeed > tmp ? *currentMaxSpeed : tmp;*/
  float TruncatedH = values[0] < EPS ? EPS : values[0];
  float u = values[1]/TruncatedH;
  float v = values[2]/TruncatedH;
  if (sqrt(u*u + v*v) > sqrt(currentMaxSpeed[1]*currentMaxSpeed[1]+currentMaxSpeed[2]*currentMaxSpeed[2])) {
    currentMaxSpeed[0] = values[0];
    currentMaxSpeed[1] = u;
    currentMaxSpeed[2] = v;
    currentMaxSpeed[3] = values[3];
  }
}
