inline void initEta_formula(float *coords, float *values, const double *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val = 1.5f*exp(-0.000001f*(x+9990.0f)*(x+9990.0f));;
  values[0] += val;
}
