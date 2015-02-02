inline void getTotalVol(const float* cellVolume, const float* value, float* totalVol) {
  (*totalVol) += (*cellVolume) * value[0];
}
