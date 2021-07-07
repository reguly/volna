inline void getMaxElevation(const float* values, float* currentMaxElevation) {
  float tmp = values[0]+values[3];
  float tmpmax = currentMaxElevation[0]+currentMaxElevation[3];
  
  if (tmp > tmpmax ) {
  //currentMaxElevation[0] = currentMaxElevation[0] > tmp ? currentMaxElevation[0] : tmp;
  currentMaxElevation[0] = values[0];
  currentMaxElevation[1] = values[1];
  currentMaxElevation[2] = values[2];
  currentMaxElevation[3] = values[3];
  } //if (values[0] > currentMaxElevation[0]) {
  //if (values[3] > 0.0f){
    //currentMaxElevation[0] = EPS;
    //currentMaxElevation[1] = values[1];
    //currentMaxElevation[2] = values[2];
    //currentMaxElevation[3] = values[3];
  //}
}
