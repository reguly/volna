inline void EvolveValuesRK3_4(const float *dT, const float *midPointConservative3, //OP_READ
            const float *Conservative, //OP_READ //temp
            float *values_new) //OP_WRITE
{
  //call to ToConservativeVariables inlined
  values_new[0] = Conservative[0] + (0.5f)* *dT * midPointConservative3[0];
  values_new[1] = Conservative[1] + (0.5f)* *dT * midPointConservative3[1];
  values_new[2] = Conservative[2] + (0.5f)* *dT * midPointConservative3[2];
  values_new[3] = Conservative[3];

  values_new[0] = values_new[0] < EPS ? EPS : values_new[0];
  //call to ToPhysicalVariables inlined
  float TruncatedH = values_new[0] < EPS ? EPS : values_new[0];
  values_new[0] = values_new[0];
  values_new[1] = values_new[1] / TruncatedH;
  values_new[2] = values_new[2] / TruncatedH;
  values_new[3] = values_new[3];
}
