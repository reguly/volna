inline void initBathymetry_formula(float *coords, float *values, const double *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val = .2f*(-5.0f-x)*(x<0.0f)-(x>=0.0f)+.2f*(t<1.0f)*exp(-(x+3.0f-2.0f*t)*(x+3.0f-2.0f*t)-y*y)+.2f*(t>=1.0f)*exp(-(x+1.0f)*(x+1.0f)-y*y);;
  values[3] = val;
}