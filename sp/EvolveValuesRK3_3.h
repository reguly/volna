inline void EvolveValuesRK3_3(const float *dT, const float *values, //OP_READ
            const float *midPointConservative, //OP_READ
            float *Conservative, //OP_RW //temp
            float *midPoint) //OP_WRITE
{
  Conservative[0] = Conservative[0]/3.0f;
  Conservative[1] = Conservative[1]/3.0f;
  Conservative[2] = Conservative[2]/3.0f;
  //call to ToConservativeVariables inlined
  Conservative[0] += ((2.0f * values[0])/3.0f) + ((*dT * midPointConservative[0])/6.0f);
  Conservative[1] += ((2.0f * (values[0]* values[1]))/3.0f) + (( *dT * midPointConservative[1])/6.0f);
  Conservative[2] += ((2.0f * (values[0]* values[2]))/3.0f) + (( *dT * midPointConservative[2])/6.0f);;
  Conservative[3] = Conservative[3];

  Conservative[0] = Conservative[0] <= EPS ? EPS : Conservative[0];
  //call to ToPhysicalVariables inlined
  float TruncatedH = Conservative[0] < EPS ? EPS : Conservative[0];
  midPoint[0] = Conservative[0];
  midPoint[1] = Conservative[1] / TruncatedH;
  midPoint[2] = Conservative[2] / TruncatedH;
  midPoint[3] = Conservative[3];
}
