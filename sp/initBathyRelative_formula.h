inline void initBathyRelative_formula(const float *coords, float *values, const float *bathy0, const float *time) {
  float x = coords[0];
  float y = coords[1];
  float t = *time;
  float val = exp(-(2.f*sqrt(x*0.01f*0.01f/(tan((5.7f*2.f*3.14159265358979323846f)/360.f)))-sqrt(g)*0.01f*t)*(2.f*sqrt(x*0.01f*0.01f/(tan((5.7f*2.f*3.14159265358979323846f)/360.f)))-sqrt(g)*0.01f*t));;
  *values = *bathy0 + val;
}
