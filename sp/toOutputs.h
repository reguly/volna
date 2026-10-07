inline void toOutputs(const float *values,
            const float *zmin,
            float *output )
{
  float truncatedH = values[0] > EPS ? values[0] : EPS;
  output[0] = values[3] + *zmin;
  output[1] = values[1]/truncatedH;
  output[2] = values[2]/truncatedH;
  output[3] = values[3] - values[0] + *zmin;
  // Output inundated cells flow depth
  // Inundated threshold depth set at 1e-3
  if ((values[3] - values[0] + *zmin) > 0.0f){
     if (values[0] > 1e-3){
         output[4] = values[0];
     } else {
         output[4] = 0.0f;
     }
  } else {
    output[4] = 0.0f;
  } 
}
