inline void getMaxElevation(float* values, float* currentMaxElevation) {
  /*float tmp = values[0]+values[3];
  *currentMaxElevation = *currentMaxElevation > tmp ? *currentMaxElevation : tmp;*/
  if (values[0] > currentMaxElevation[0]) {
    currentMaxElevation[0] = values[0];
    currentMaxElevation[1] = values[1];
    currentMaxElevation[2] = values[2];
    currentMaxElevation[3] = values[3];
  }
}
