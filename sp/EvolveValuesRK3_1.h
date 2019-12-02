inline void EvolveValuesRK3_1(const float *dT, const float *midPointConservative, //OP_RW //temp
            const float *in, //OP_READ
            float *Conservative, //OP_WRITE //temp
            float *midPoint) //OP_WRITE
{
  //call to ToConservativeVariables inlined
 
  Conservative[0] = in[0] + (0.5f)* *dT * midPointConservative[0];
  Conservative[1] = in[1] + (0.5f)* *dT * midPointConservative[1];
  Conservative[2] = in[2] + (0.5f)* *dT * midPointConservative[2];
  Conservative[3] = in[3]; 
  
  Conservative[0] = Conservative[0] <= EPS ? EPS : Conservative[0];
  //call to ToPhysicalVariables inlined
  float TruncatedH = Conservative[0] < EPS ? EPS : Conservative[0];
  midPoint[0] = Conservative[0];
  midPoint[1] = Conservative[1] / TruncatedH;
  midPoint[2] = Conservative[2] / TruncatedH;
  midPoint[3] = Conservative[3];
}
