inline void computeFluxes1ComputeCellValues(const float *cellLeft, const float *cellRight,
                                const float *alphaleft, const float *alpharight,
                                const float *edgeLength, const float *edgeNormals,
                                const float *leftcellCenters, const float *rightcellCenters,
                                const float *edgeCenters,
                                const float *leftGradient, const float *rightGradient,
                                const int *isRightBoundary, //OP_READ
                                float *bathySource, //OP_WRITE
                                float* leftCellValues,  //OP_WRITE
                                float* rightCellValues) //OP_WRITE
{
  //begin EdgesValuesFromCellValues
  float InterfaceBathy;
  leftCellValues[0] = cellLeft[0];
  leftCellValues[1] = cellLeft[1];
  leftCellValues[2] = cellLeft[2];
  leftCellValues[3] = cellLeft[3];
  float dxl, dyl, dxr, dyr;
  dxl = (edgeCenters[0] - leftcellCenters[0]);
  dyl = (edgeCenters[1] - leftcellCenters[1]);

  dxr = (edgeCenters[0] - rightcellCenters[0]);
  dyr = (edgeCenters[1] - rightcellCenters[1]);
  
  if (!*isRightBoundary) {
    rightCellValues[0] = cellRight[0];
    rightCellValues[1] = cellRight[1];
    rightCellValues[2] = cellRight[2];
    rightCellValues[3] = cellRight[3];
    
  } else {
    rightCellValues[3] = cellLeft[3];
    float nx = edgeNormals[0];
    float ny = edgeNormals[1];
    float inNormalVelocity = cellLeft[1] * nx + cellLeft[2] * ny;
    float inTangentVelocity = -1.0f *  cellLeft[1] * ny + cellLeft[2] * nx;
    float outNormalVelocity = 0.0f;
    float outTangentVelocity = 0.0f;

    //WALL
    rightCellValues[0] = cellLeft[0];
    outNormalVelocity =  -1.0f*inNormalVelocity;
    outTangentVelocity = inTangentVelocity;

    /* //HEIGHTSUBC
     rightCellValues[0] = -1.0 * rightCellValues[3];
     rightCellValues[0] += 0.1 * sin(10.0*t);
     outNormalVelocity = inNormalVelocity;
     outNormalVelocity -=
     2.0 * sqrt( g * cellLeft[0] );
     outNormalVelocity +=
     2.0 * sqrt( g * rightCellValues[0] );

     outTangentVelocity = inTangentVelocity;
     */ //end HEIGHTSUBC

    /* //FLOWSUBC
     outNormalVelocity = 1.0f;

     //rightCellValues[0] = - rightCellValues[3];

     rightCellValues[0] = (inNormalVelocity - outNormalVelocity);
     rightCellValues[0] *= .5 / sqrt( g );

     rightCellValues[0] += sqrt( cellLeft[0] );

     outTangentVelocity = inTangentVelocity;
     */ 
    rightCellValues[1] = outNormalVelocity * nx - outTangentVelocity * ny;
    rightCellValues[2] = outNormalVelocity * ny + outTangentVelocity * nx;
  }

  // ------------------------------------------------------------------------------------
  // Second order Reconstruction
  if (!*isRightBoundary) {
    leftCellValues[0] +=  alphaleft[0] * ((dxl * leftGradient[0])+(dyl * leftGradient[1]));
    leftCellValues[0] = leftCellValues[0] > 0.0f ? leftCellValues[0] : 0.0f;

    leftCellValues[3] += alphaleft[0] * ((dxl * leftGradient[6])+(dyl * leftGradient[7]));
    leftCellValues[1] += alphaleft[0] * ((dxl * leftGradient[2])+(dyl * leftGradient[3]));
    leftCellValues[2] += alphaleft[0] * ((dxl * leftGradient[4])+(dyl * leftGradient[5]));
  
    rightCellValues[0] +=  alpharight[0] * ((dxr * rightGradient[0])+(dyr * rightGradient[1]));
    rightCellValues[0] = rightCellValues[0] > 0.0f ? rightCellValues[0] : 0.0f;
    rightCellValues[3] += alpharight[0] * ((dxr * rightGradient[6])+(dyr * rightGradient[7]));
    rightCellValues[1] += alpharight[0] * ((dxr * rightGradient[2])+(dyr * rightGradient[3]));
    rightCellValues[2] += alpharight[0] * ((dxr * rightGradient[4])+(dyr * rightGradient[5]));
  }
  // Audusse Reconstruction(2004) 1st order Source Discretization
  InterfaceBathy = leftCellValues[3] > rightCellValues[3] ? leftCellValues[3] : rightCellValues[3];
  bathySource[0] =.5f * g * (leftCellValues[0]*leftCellValues[0]);
  bathySource[1] =.5f * g * (rightCellValues[0]*rightCellValues[0]);
  float hL = (leftCellValues[0] + leftCellValues[3] - InterfaceBathy);
  hL = hL > 0.0f? hL : 0.0f;
  float hR = (rightCellValues[0] + rightCellValues[3] - InterfaceBathy);
  hR = hR > 0.0f ? hR : 0.0f;
  bathySource[0] -= .5f * g * (hL * hL);
  bathySource[1] -= .5f * g * (hR * hR);
  // Audusse Reconstruction(2005) 2nd order Centered term
  bathySource[2] = -.5f * g *(leftCellValues[0] + cellLeft[0])*(leftCellValues[3] - cellLeft[3]);
  bathySource[3] = -.5f * g *(rightCellValues[0] + cellRight[0])*(rightCellValues[3] - cellRight[3]);
  
  leftCellValues[0] = hL;
  rightCellValues[0] = hR;
  
  bathySource[0] *= *edgeLength;
  bathySource[1] *= *edgeLength;
  bathySource[2] *= *edgeLength;
  bathySource[3] *= *edgeLength;
}