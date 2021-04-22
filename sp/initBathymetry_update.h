inline void initBathymetry_update(float *values, const float *zmin, const int *firstTime) {
    if (*firstTime){
      if (values[3] >0.0f){
        values[0] = EPS;
        values[3] = -1.0f* *zmin + values[3];
      } else {
        values[0] = -1.0f*values[3];
        values[3] = -1.0f* *zmin;
      }
    } else {
      if (values[3] > 0.0f){
        values[3] = -1.0f* *zmin + values[3] + values[0];
      } else {
        values[3] = (values[0]+values[3]) -1.0f* *zmin;
      }

    }
}
