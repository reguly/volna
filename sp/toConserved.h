inline void toConserved(float *center)
{

 if (center[0] < EPS){
   center[1] = 0.0f;
   center[2] = 0.0f;
 } else {
   center[1] *= center[0];
   center[2] *= center[0];
 }
}
