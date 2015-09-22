inline void getMaxSpeed(const float* values, float* currentMaxSpeed) {
  /*float tmp = values[0]+values[3];
  *currentMaxSpeed = *currentMaxSpeed > tmp ? *currentMaxSpeed : tmp;*/
  if (sqrt(values[1]*values[1]+values[2]*values[2]) > sqrt(currentMaxSpeed[1]*currentMaxSpeed[1]+currentMaxSpeed[2]*currentMaxSpeed[3])) {
    currentMaxSpeed[0] = values[0];
    currentMaxSpeed[1] = values[1];
    currentMaxSpeed[2] = values[2];
    currentMaxSpeed[3] = values[3];
  }
}
