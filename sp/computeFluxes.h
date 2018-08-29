inline void computeFluxes(const float *cellLeft, const float *cellRight,
                                const float *alphaleft, const float *alpharight,
                                const float *edgeLength, const float *edgeNormals,
                                const float *leftcellCenters, const float *rightcellCenters,
                                const float *edgeCenters,
                                const float *leftGradient, const float *rightGradient,
                                const int *isRightBoundary, //OP_READ
                                float *bathySource, float *out, //OP_WRITE
                                float *maxEdgeEigenvalues) //OP_WRITE
{
  //begin EdgesValuesFromCellValues
  float leftCellValues[4];
  float rightCellValues[4];
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
 
  // ------------------------------------------------------------------------------------
  // HLL Riemann Solver
  // Estimation of the wave speeds at the interface.
  float cL = sqrt(g * leftCellValues[0]);
  cL = cL > 0.0f ? cL : 0.0f;
  float cR = sqrt(g * rightCellValues[0]);
  cR = cR > 0.0f ? cR : 0.0f;
  float uLn = leftCellValues[1] * edgeNormals[0] + leftCellValues[2] * edgeNormals[1];
  float uRn = rightCellValues[1] * edgeNormals[0] + rightCellValues[2] * edgeNormals[1];
  
  float unStar = 0.5f * (uLn + uRn) + (cL-cR);
  float cStar = 0.5f * (cL + cR) - 0.25f* (uRn-uLn);
  float sL = (uLn - cL) < (unStar - cStar) ? (uLn - cL) : (unStar - cStar);
  float sR = (uRn + cR) > (unStar + cStar) ? (uRn + cR) : (unStar + cStar);

  // sStar is needed for the HLLC extension.
  float sStar;
  sStar = (sL*rightCellValues[0]*(uRn - sR) - sR*leftCellValues[0]*(uLn - sL))/
          (rightCellValues[0]*(uRn - sR) - leftCellValues[0]*(uLn - sL));

  if ((leftCellValues[0] <= EPS) && (rightCellValues[0] > EPS)) {
      sL = uRn - 2.0f*cR;
      sR = uRn + cR;
      sStar = sL;
  }
  if ((rightCellValues[0] <= EPS) && (leftCellValues[0] > EPS)) {
      sR = uLn + 2.0f*cL;
      sL =  uLn - cL;
      sStar = sR;
  }

  float sLMinus = sL < 0.0f ? sL : 0.0f;
  float sRPlus = sR > 0.0f ? sR : 0.0f;
  float sRMinussL = sRPlus - sLMinus;
  sRMinussL = sRMinussL < EPS ?  EPS : sRMinussL;
  //-------------------------------------------------
  float t1 = sRPlus / sRMinussL;
  //assert( ( 0 <= t1 ) && ( t1 <= 1 ) );

  float t2 = ( -1.0 * sLMinus ) / sRMinussL;
  //assert( ( 0 <= t2 ) && ( t2 <= 1 ) );

  float t3 = ( sRPlus * sLMinus ) / sRMinussL;
  // ------------------------------------------------------------------------------------
  // Velocities parallel to the interface.
  float uLp = leftCellValues[2]*edgeNormals[0] - leftCellValues[1]*edgeNormals[1];
  float uRp = rightCellValues[2]*edgeNormals[0] - rightCellValues[1]*edgeNormals[1];

  float LeftFluxes_H, LeftFluxes_U, LeftFluxes_V, LeftFluxes_N;
  //inlined ProjectedPhysicalFluxes(leftCellValues, Normals, params, LeftFluxes);
  float HuDotN = (leftCellValues[0] * leftCellValues[1]) * edgeNormals[0] +
  (leftCellValues[0] * leftCellValues[2]) * edgeNormals[1];

  LeftFluxes_H = HuDotN;
  LeftFluxes_U = HuDotN * leftCellValues[1];
  LeftFluxes_V = HuDotN * leftCellValues[2];
  // Normal Momentum flux term
  LeftFluxes_N = HuDotN * uLn;

  LeftFluxes_U += (.5f * g * edgeNormals[0] ) * ( leftCellValues[0] * leftCellValues[0] );
  LeftFluxes_V += (.5f * g * edgeNormals[1] ) * ( leftCellValues[0] * leftCellValues[0] );
  LeftFluxes_N += (.5f * g ) * ( leftCellValues[0] * leftCellValues[0] );
  //end of inlined

  float RightFluxes_H,RightFluxes_U, RightFluxes_V, RightFluxes_N;
  //inlined ProjectedPhysicalFluxes(rightCellValues, Normals, params, RightFluxes);
  HuDotN = (rightCellValues[0] * rightCellValues[1] * edgeNormals[0]) +
  (rightCellValues[0] * rightCellValues[2] * edgeNormals[1]);

  RightFluxes_H =   HuDotN;
  RightFluxes_U =   HuDotN * rightCellValues[1];
  RightFluxes_V =   HuDotN * rightCellValues[2];
  // Normal Momentum flux term
  RightFluxes_N =   HuDotN * uRn;

  RightFluxes_U += (.5f * g * edgeNormals[0] ) * ( rightCellValues[0] * rightCellValues[0] );
  RightFluxes_V += (.5f * g * edgeNormals[1] ) * ( rightCellValues[0] * rightCellValues[0] );
  RightFluxes_N += (.5f * g) * ( rightCellValues[0] * rightCellValues[0] );

  float FStar[3];
  FStar[0] =
  ( t1 * LeftFluxes_H ) +
  ( t2 * RightFluxes_H ) +
  ( t3 * ( rightCellValues[0] - leftCellValues[0] ) );


  FStar[1] =
  ( t1 * LeftFluxes_N ) +
  ( t2 * RightFluxes_N ) +
  ( t3 * ( (rightCellValues[0] * uRn) -
          (leftCellValues[0] * uLn) ) );

  // ------------------------------------------------------------------------
  // HLLC Flux Solver
  if( sL >= 0.0f) {
    out[0] = LeftFluxes_H;
    out[1] = LeftFluxes_U;
    out[2] = LeftFluxes_V;
  } else if ((sL < 0.0f) && (sStar >= 0.0f)){
    out[0] = FStar[0];
    FStar[2] = FStar[0] * uLp;
    out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
    out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  } else if((sStar < 0.0f) && (sR >= 0.0f)){
    out[0] = FStar[0];
    FStar[2] = FStar[0] * uRp;
    out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
    out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  } else {
    out[0] = RightFluxes_H;
    out[1] = RightFluxes_U;
    out[2] = RightFluxes_V;
  }
  out[0] *= *edgeLength;
  out[1] *= *edgeLength;
  out[2] *= *edgeLength;
  float maximum = fabs(uLn + cL);
  maximum = maximum > fabs(uLn - cL) ? maximum : fabs(uLn - cL);
  maximum = maximum > fabs(uRn + cR) ? maximum : fabs(uRn + cR);
  maximum = maximum > fabs(uRn - cR) ? maximum : fabs(uRn - cR);
  *maxEdgeEigenvalues = maximum;
}
