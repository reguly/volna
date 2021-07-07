inline void gatherLocations(const float *values, const float *zmin, float *dest) {
  dest[0] = values[0] + (values[3] - values[0] + *zmin);
  dest[1] = values[0];
  dest[2] = values[1]/values[0];
  dest[3] = values[2]/values[0];
  dest[4] = values[3];
}
