inline void incConst(const float *in, float *out, const int *variables) {
  float H;
  if (*variables & 1) {
    out[0] += *in;
    out[3] += *in;
  }
  if (*variables & 2) {
    H = out[0] > EPS ? out[0] : EPS;
    out[1] += *in * H;
  }
  if (*variables & 4) {
    H = out[0] > EPS ? out[0] : EPS;
    out[2] += *in * H;
  }
  if (*variables & 8) {
    out[3] += *in;
  }
}
