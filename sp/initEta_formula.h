inline void initEta_formula(const float *coords, float *values, const double *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val =   0.1f*(0.9756f/(1.f-0.2195f) -1.0f - ((x-2.f)*(x-2.f) +(y-2.f)*(y-2.f))/(1.f) *
   (0.9518f/((1.0f - 0.2195f)*(1.0f - 0.2195f))   -1.0f) )    ;;
  values[0] += val;
}