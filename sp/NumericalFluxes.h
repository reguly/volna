inline void NumericalFluxes(float *maxEdgeEigenvalues0,
			    float *maxEdgeEigenvalues1, 
			    float *maxEdgeEigenvalues2, 
			    float *EdgeVolumes0,
			    float *EdgeVolumes1,
			    float *EdgeVolumes2,
			    float *cellVolumes, //OP_READ
            float *zeroInit, float *minTimeStep ) //OP_MIN
{
  float local = 0.0f;
  local += *maxEdgeEigenvalues0 * *(EdgeVolumes0);
  local += *maxEdgeEigenvalues1 * *(EdgeVolumes1);
  local += *maxEdgeEigenvalues2 * *(EdgeVolumes2);
  zeroInit[0] = 0.0f;
  zeroInit[1] = 0.0f;
  zeroInit[2] = 0.0f;
  zeroInit[3] = 0.0f;

  *minTimeStep = min(*minTimeStep, 2.0f * *cellVolumes / local);
}