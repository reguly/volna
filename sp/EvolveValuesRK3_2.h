inline void EvolveValuesRK3_2(const float *dT, const float *midPointConservative1, //OP_READ
            float *Conservative) //OP_RW //temp)
{

  //call to ToConservativeVariables inlined
  Conservative[0] +=  (0.5f)* *dT * midPointConservative1[0];
  Conservative[1] +=  (0.5f)* *dT * midPointConservative1[1];
  Conservative[2] +=  (0.5f)* *dT * midPointConservative1[2];
  Conservative[3] = Conservative[3];

  Conservative[0] = Conservative[0] <= EPS ? EPS : Conservative[0];
}
