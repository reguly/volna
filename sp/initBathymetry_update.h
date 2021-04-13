inline void initBathymetry_update(float *values, const float *zmin, const int *firstTime) {
    if (*firstTime){
      //op_printf("values %g %g \n", values[0], values[3]);
      values[0] -= values[3];
      values[0] = values[0] < EPS ? EPS : values[0];
      values[3] += fabs(*zmin) + values[0];
      //op_printf("After values %g %g \n", values[0], values[3]);
    } else {
    values[0] = values[0] < EPS ? EPS : values[0];
    values[3] += fabs(*zmin) + values[0];
    }
   //op_printf("Second Time values %g %g \n", values[0], values[3]);
}
