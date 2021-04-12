inline void initBathymetry_update(float *values, const float *zmin, const int *firstTime) {
  if (*firstTime)
    values[0] -= values[3];
    values[0] = values[0] < EPS ? EPS : values[0];
    //values[1] *= values[0];
    //values[2] *= values[0];
    values[3] += fabs(*zmin) + values[0];

  values[0] = values[0] < EPS ? EPS : values[0];
}
