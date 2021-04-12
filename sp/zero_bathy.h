inline void zero_bathy(const float *values,
            float *zmin ) //OP_MIN
{
  *zmin = MIN(*zmin, values[3]);
}
