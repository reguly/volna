inline void initU_formula(const float *coords, float *values, const double *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val = 0.0f;

  // if (values[0] > EPS) {
  //    val = -1.0f*sqrt(9.81f) * (values[0] + values[3]);
  //    val *= values[0];
  // } else {
  //    val = 0.0f;
  // }
	//
	//
  values[1] += val;
}
