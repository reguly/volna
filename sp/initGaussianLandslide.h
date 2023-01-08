inline void initGaussianLandslide(const float *center, float *values, const float *mesh_xmin, const float *A, const float *t, const float *lx, const float *ly, const float *v) {
  float x = center[0];
  float y = center[1];
  values[3] = (*mesh_xmin-x)*(x<0.0f)-5.0f*(x>=0.0f)+
      *A*(*t<1.0f/(*v))*exp(-1.0f* *lx* *lx*(x+3.0f-*v**t)*(x+3.0f-*v**t)-*ly**ly*y*y)
      +*A*(*t>=1.0f/(*v))*exp(-*lx*(x+3.0f-1.0f)**lx*(x+3.0f-1.0f)-*ly**ly*y*y);
}
