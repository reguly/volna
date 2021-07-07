inline void initBathymetry_update(float *values, const float *z_zero, const float *zmin, const int *firstTime) {
    if (*firstTime){
      if (*z_zero > 0.0f){
        values[0] = EPS;
        values[3] = -1.0f* *zmin + *z_zero;
      } else {
        values[0] = -1.0f* *z_zero;
        values[3] = -1.0f* *zmin;
      }
    } else {
      if (*z_zero > 0.0f){
        values[3] = -1.0f* *zmin + *z_zero + values[0];
      } else {
        values[3] = values[0] + *z_zero - *zmin;
      }
    }
}
