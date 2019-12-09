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
  float uL, vL, uR, vR;
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
    float HL = leftCellValues[0] > 1.0f ? leftCellValues[0] : 1.0f;
    uL = leftCellValues[1]/ HL;
    vL = leftCellValues[2]/ HL;
    float inNormalVelocity = uL * nx + vL * ny;
    float inTangentVelocity = -1.0f *  uL * ny + vL * nx;
    float outNormalVelocity = 0.0f;
    float outTangentVelocity = 0.0f;


    //Outflow
    rightCellValues[0] = leftCellValues[0];
    if(leftCellValues[3] < 0.0f){
      rightCellValues[1] = uL;
      rightCellValues[2] = vL;
    } else {
    // Wall
      outNormalVelocity =  -1.0f*inNormalVelocity;
      outTangentVelocity = inTangentVelocity;
      rightCellValues[1] = outNormalVelocity * nx - outTangentVelocity * ny;
      rightCellValues[2] = outNormalVelocity * ny + outTangentVelocity * nx;
    }

    rightCellValues[3] = leftCellValues[3];

    rightCellValues[1] *= rightCellValues[0];
    rightCellValues[2] *= rightCellValues[0];

  }

  // ------------------------------------------------------------------------------------
  // Second order Reconstruction
  if (!*isRightBoundary){
    leftCellValues[0] +=  alphaleft[0] * ((dxl * leftGradient[0])+(dyl * leftGradient[1]));
    //leftCellValues[0] = leftCellValues[0] > EPS ? leftCellValues[0] : EPS;

    leftCellValues[3] += alphaleft[3] * ((dxl * leftGradient[6])+(dyl * leftGradient[7]));
    leftCellValues[1] += alphaleft[1] * ((dxl * leftGradient[2])+(dyl * leftGradient[3]));
    leftCellValues[2] += alphaleft[2] * ((dxl * leftGradient[4])+(dyl * leftGradient[5]));

    rightCellValues[0] +=  alpharight[0] * ((dxr * rightGradient[0])+(dyr * rightGradient[1]));
    //rightCellValues[0] = rightCellValues[0] > EPS ? rightCellValues[0] : EPS;
    rightCellValues[3] += alpharight[3] * ((dxr * rightGradient[6])+(dyr * rightGradient[7]));
    rightCellValues[1] += alpharight[1] * ((dxr * rightGradient[2])+(dyr * rightGradient[3]));
    rightCellValues[2] += alpharight[2] * ((dxr * rightGradient[4])+(dyr * rightGradient[5]));
  }

  //float TruncatedHL = leftCellValues[0] > EPS ? leftCellValues[0] : EPS;
  //float TruncatedHR = rightCellValues[0] > EPS ? rightCellValues[0] : EPS;
  float TruncatedHL = leftCellValues[0] > 1.0f ? leftCellValues[0] : 1.0f;
  float TruncatedHR = rightCellValues[0] > 1.0f ? rightCellValues[0] : 1.0f;
  uL = leftCellValues[1]/TruncatedHL;
  vL = leftCellValues[2]/TruncatedHL;

  uR = rightCellValues[1]/TruncatedHR;
  vR = rightCellValues[2]/TruncatedHR;

  // Audusse Reconstruction(2004) 1st order Source Discretization
  InterfaceBathy = leftCellValues[3] > rightCellValues[3] ? leftCellValues[3] : rightCellValues[3];
  bathySource[0] =.5f * g * (leftCellValues[0]*leftCellValues[0]);
  bathySource[1] =.5f * g * (rightCellValues[0]*rightCellValues[0]);
  float hL = (leftCellValues[0] + leftCellValues[3] - InterfaceBathy);
  hL = hL > EPS ? hL : EPS;
  float hR = (rightCellValues[0] + rightCellValues[3] - InterfaceBathy);
  hR = hR > EPS ? hR : EPS;
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

  float uLn = uL * edgeNormals[0] + vL * edgeNormals[1];
  float uRn = uR * edgeNormals[0] + vR * edgeNormals[1];


  float unStar = 0.5f * (uLn + uRn) + (cL-cR);
  float cStar = 0.5f * (cL + cR) - 0.25f* (uRn-uLn);
  float sL = (uLn - cL) < (unStar - cStar) ? (uLn - cL) : (unStar - cStar);
  float sR = (uRn + cR) > (unStar + cStar) ? (uRn + cR) : (unStar + cStar);

  // sStar is needed for the HLLC extension.
  /*float sStar;
  sStar = (sL*rightCellValues[0]*(uRn - sR) - sR*leftCellValues[0]*(uLn - sL))/
          (rightCellValues[0]*(uRn - sR) - leftCellValues[0]*(uLn - sL));

  */
  if ((cellLeft[0] <= EPS) && (cellRight[0] > EPS)) {
      sL = uRn - 2.0f*cR;
      sR = uRn + cR;
  }
  if ((cellRight[0] <= EPS) && (cellLeft[0] > EPS)) {
      sR = uLn + 2.0f*cL;
      sL =  uLn - cL;
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
  //float uLp = leftCellValues[2]*edgeNormals[0] - leftCellValues[1]*edgeNormals[1];
  //float uRp = rightCellValues[2]*edgeNormals[0] - rightCellValues[1]*edgeNormals[1];

  float LeftFluxes_H, LeftFluxes_U, LeftFluxes_V, LeftFluxes_N;
  //inlined ProjectedPhysicalFluxes(leftCellValues, Normals, params, LeftFluxes);
  float HuDotN = (leftCellValues[1]) * edgeNormals[0] +
  (leftCellValues[2]) * edgeNormals[1];

  LeftFluxes_H = HuDotN;
  LeftFluxes_U = HuDotN * uL;
  LeftFluxes_V = HuDotN * vL;
  // Normal Momentum flux term
  LeftFluxes_N = HuDotN * uLn;

  LeftFluxes_U += (.5f * g * edgeNormals[0] ) * ( leftCellValues[0] * leftCellValues[0] );
  LeftFluxes_V += (.5f * g * edgeNormals[1] ) * ( leftCellValues[0] * leftCellValues[0] );
  LeftFluxes_N += (.5f * g ) * ( leftCellValues[0] * leftCellValues[0] );
  //end of inlined

  float RightFluxes_H,RightFluxes_U, RightFluxes_V, RightFluxes_N;
  //inlined ProjectedPhysicalFluxes(rightCellValues, Normals, params, RightFluxes);
  HuDotN = (rightCellValues[1]) * edgeNormals[0] +
  (rightCellValues[2]) * edgeNormals[1];

  RightFluxes_H =   HuDotN;
  RightFluxes_U =   HuDotN * uR;
  RightFluxes_V =   HuDotN * vR;
  // Normal Momentum flux term
  RightFluxes_N =   HuDotN * uRn;

  RightFluxes_U += (.5f * g * edgeNormals[0] ) * ( rightCellValues[0] * rightCellValues[0] );
  RightFluxes_V += (.5f * g * edgeNormals[1] ) * ( rightCellValues[0] * rightCellValues[0] );
  RightFluxes_N += (.5f * g) * ( rightCellValues[0] * rightCellValues[0] );

  /*float FStar[3];
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
  }*/

  out[0] =
  ( t1 * LeftFluxes_H ) +
  ( t2 * RightFluxes_H ) +
  ( t3 * ( rightCellValues[0] - leftCellValues[0] ) );

  out[1] =
  ( t1 * LeftFluxes_U ) +
  ( t2 * RightFluxes_U ) +
  ( t3 * ( (rightCellValues[1]) -
          (leftCellValues[1]) ) );


  out[2] =
  ( t1 * LeftFluxes_V ) +
  ( t2 * RightFluxes_V ) +
  ( t3 * ( (rightCellValues[2]) -
          (leftCellValues[2]) ) );


  out[0] *= *edgeLength;
  out[1] *= *edgeLength;
  out[2] *= *edgeLength;
  float maximum = fabs(uLn + cL);
  maximum = maximum > fabs(uLn - cL) ? maximum : fabs(uLn - cL);
  maximum = maximum > fabs(uRn + cR) ? maximum : fabs(uRn + cR);
  maximum = maximum > fabs(uRn - cR) ? maximum : fabs(uRn - cR);
  *maxEdgeEigenvalues = maximum;
}
