inline void gatherLocations(const float *values, float *dest) {
	dest[0] = values[0] + values[3];
  dest[1] = values[0];
  dest[2] = values[1];
  dest[3] = values[2];
  dest[4] = values[3];
}
