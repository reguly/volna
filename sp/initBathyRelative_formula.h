inline void initBathyRelative_formula(const float *coords, float *values, const float *bathy0, const double *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val = 0.0f;;
  values[3] += *bathy0 + val;
}