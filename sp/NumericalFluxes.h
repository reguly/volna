inline void NumericalFluxes(float **maxEdgeEigenvalues, float **EdgeVolumes, float *cellVolumes, //OP_READ
                            float *minTimeStep ) //OP_MIN
{
  float local = 0.0f;
  for (int j = 0; j < 3; j++) {
    local += *maxEdgeEigenvalues[j] * *(EdgeVolumes[j]);
  }

  *minTimeStep = MIN(*minTimeStep, 2.0f * *cellVolumes / local);
}
