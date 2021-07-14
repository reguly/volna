inline void computeFluxes_sph(const float *cellLeft, const float *cellRight,
                                const float *alphaleft, const float *alpharight,
                                const float *edgeLength, const float *edgeNormals,
                                const float *leftcellCenters, const float *rightcellCenters,
                                const float *edgeCenters,
                                const float *leftGradient, const float *rightGradient,
                                const int *isRightBoundary, //OP_READ
                                float *bathySource, float *out, //OP_WRITE
                                float *maxEdgeEigenvalues, const float *zmin) //OP_WRITE
{
  //begin EdgesValuesFromCellValues
  float leftCellValues[4];
  float rightCellValues[4];
  float InterfaceBathy;
  float zL, zR;
  float uR,vR,uL,vL;
  leftCellValues[0] = cellLeft[0];
  leftCellValues[1] = cellLeft[1];
  leftCellValues[2] = cellLeft[2];
  leftCellValues[3] = cellLeft[3];
  float dxl, dyl, dxr, dyr;
  // Convert to arc length
  dxl = (edgeCenters[0] - leftcellCenters[0])*Radius;
  dyl = (edgeCenters[1] - leftcellCenters[1])*Radius;
  dxr = (edgeCenters[0] - rightcellCenters[0])*Radius;
  dyr = (edgeCenters[1] - rightcellCenters[1])*Radius;

  // Second order Reconstruction
  // leftCellValues[0] += alphaleft[0] * ((dxl * leftGradient[0])+(dyl * leftGradient[1]));
  // leftCellValues[1] += alphaleft[0] * ((dxl * leftGradient[2])+(dyl * leftGradient[3]));
  // leftCellValues[2] += alphaleft[0] * ((dxl * leftGradient[4])+(dyl * leftGradient[5]));
  // leftCellValues[3] += alphaleft[0] * ((dxl * leftGradient[6])+(dyl * leftGradient[7]));
  if (leftCellValues[0] >= EPS){
     uL = leftCellValues[1]/leftCellValues[0];
     vL = leftCellValues[2]/leftCellValues[0];
  } else {
     uL = 0.0f;
     vL = 0.0f;
  }
  zL = cellLeft[3] - cellLeft[0];
  zR = cellRight[3] - cellRight[0];

  if (!*isRightBoundary) {
    rightCellValues[0] = cellRight[0];
    rightCellValues[1] = cellRight[1];
    rightCellValues[2] = cellRight[2];
    rightCellValues[3] = cellRight[3];

    // Second order Reconstruction
    // rightCellValues[0] += alpharight[0] * ((dxr * rightGradient[0])+(dyr * rightGradient[1]));
    // rightCellValues[1] += alpharight[0] * ((dxr * rightGradient[2])+(dyr * rightGradient[3]));
    // rightCellValues[2] += alpharight[0] * ((dxr * rightGradient[4])+(dyr * rightGradient[5]));
    // rightCellValues[3] += alpharight[0] * ((dxr * rightGradient[6])+(dyr * rightGradient[7]));
    if (rightCellValues[0] >= EPS){
       uR = rightCellValues[1]/rightCellValues[0];
       vR = rightCellValues[2]/rightCellValues[0];
     } else {
       uR = 0.0f;
       vR = 0.0f;
     }
   } else {
     float nx = edgeNormals[0];
     float ny = edgeNormals[1];
     float inNormalVelocity = uL * nx + vL * ny;
     float inTangentVelocity = -1.0f * uL * ny + vL * nx;
     float outNormalVelocity = 0.0f;
     float outTangentVelocity = 0.0f;

     // Outflow
     float wet = (fabs(*zmin) - zL) > 0.0f ? (fabs(*zmin) - zL) : EPS;
     float critical = sqrt(uL*uL + vL*vL);
     outTangentVelocity = inTangentVelocity;
     if (critical < sqrt(g*leftCellValues[0])){
       rightCellValues[0] = wet;
       rightCellValues[3] = wet + zL;
       outNormalVelocity = 0.0f;
     } else {
       rightCellValues[0] = leftCellValues[0];
       rightCellValues[3] = leftCellValues[3];
       outNormalVelocity = inNormalVelocity;
     }
     /*
     // Wall
     rightCellValues[3] = leftCellValues[3];
     rightCellValues[0] = leftCellValues[0];
     outNormalVelocity =  -1.0f*inNormalVelocity;
     outTangentVelocity = inTangentVelocity;
     */
     uR = outNormalVelocity * nx - outTangentVelocity * ny;
     vR = outNormalVelocity * ny + outTangentVelocity * nx;
     rightCellValues[1] = uR*rightCellValues[0];
     rightCellValues[2] = vR*rightCellValues[0];
     zR = zL;
  }

  rightCellValues[3] -= rightCellValues[0];
  leftCellValues[3] -= leftCellValues[0];
  // ------------------------------------------------------------------------------------
  // Audusse Reconstruction(2004) 1st order Source Discretization
  // ------------------------------------------------------------------------------------
  InterfaceBathy = leftCellValues[3] > rightCellValues[3] ? leftCellValues[3] : rightCellValues[3];
  bathySource[0] =0.5f * g * (leftCellValues[0]*leftCellValues[0]);
  bathySource[1] =0.5f * g * (rightCellValues[0]*rightCellValues[0]);
  float hL = (leftCellValues[0] + leftCellValues[3] - InterfaceBathy);
  //hL = hL > EPS ? hL : EPS;
  hL = hL > 0.0f ? hL : 0.0f;
  float hR = (rightCellValues[0] + rightCellValues[3] - InterfaceBathy);
  //hR = hR > EPS ? hR : EPS;
  hR = hR > 0.0f ? hR : 0.0f;
  bathySource[0] -= .5f * g * (hL * hL);
  bathySource[1] -= .5f * g * (hR * hR);
  // Audusse Reconstruction(2005) 2nd order centered term
  // bathySource[2] = -.5f * g *(leftCellValues[0] + cellLeft[0])*(leftCellValues[3] - zL);
  // bathySource[3] = -.5f * g *(rightCellValues[0] + cellRight[0])*(rightCellValues[3] - zR);
  bathySource[2] = .5f * g * (hL * hL)*sin(M_PI*leftcellCenters[1]/180.0);
  bathySource[3] = .5f * g * (hR * hR)*sin(M_PI*rightcellCenters[1]/180.0);
  bathySource[0] *= *edgeLength;
  bathySource[1] *= *edgeLength;
  bathySource[2] *= *edgeLength;
  bathySource[3] *= *edgeLength;

  // ------------------------------------------------------------------------------------
  // HLL Riemann Solver
  // Estimation of the wave speeds at the interface.
  float cL = sqrt(g * hL);
  cL = cL > 0.0f ? cL : 0.0f;
  float cR = sqrt(g * hR);
  cR = cR > 0.0f ? cR : 0.0f;

  float uLn = uL * edgeNormals[0] + vL * edgeNormals[1];
  float uRn = uR * edgeNormals[0] + vR * edgeNormals[1];

  float unStar = 0.5f * (uLn + uRn) + (cL-cR);
  float cStar = 0.5f * (cL + cR) - 0.25f* (uRn-uLn);
  float sL = (uLn - cL) < (unStar - cStar) ? (uLn - cL) : (unStar - cStar);
  float sR = (uRn + cR) > (unStar + cStar) ? (uRn + cR) : (unStar + cStar);

  // sStar is needed for the HLLC extension.
  float sStar;
  sStar = (sL*hR*(uRn - sR) - sR*hL*(uLn - sL))/
          (hR*(uRn - sR) - hL*(uLn - sL));

  if ((hL <= EPS) && (hR > EPS)) {
      sL = uRn - 2.0f*cR;
      sR = uRn + cR;
      sStar = sL;
  }
  if ((hR <= EPS) && (hL > EPS)) {
      sR = uLn + 2.0f*cL;
      sL =  uLn - cL;
      sStar = sR;
  }
  // ------------------------------------------------------------------------------------
  // Velocities parallel to the interface.
  float uLp = vL*edgeNormals[0] - uL*edgeNormals[1];
  float uRp = vR*edgeNormals[0] - uR*edgeNormals[1];

  float LeftFluxes_H, LeftFluxes_U, LeftFluxes_V, LeftFluxes_N;
  // float HuDotN = (hL*uL) * edgeNormals[0] + (hL*vL) * edgeNormals[1];
  float HuDotN = (hL*uL) * edgeNormals[0] + (hL*vL*cos(M_PI*edgeCenters[1]/180.0)) * edgeNormals[1];
  LeftFluxes_H = HuDotN;
  LeftFluxes_U = HuDotN * uL;
  LeftFluxes_V = HuDotN * vL;
  // Normal Momentum flux term
  LeftFluxes_N = HuDotN * uLn;

  LeftFluxes_U += (.5f * g * edgeNormals[0] ) * ( hL * hL);
  LeftFluxes_V += (.5f * g * edgeNormals[1] ) * ( hL * hL* cos(M_PI*edgeCenters[1]/180.0));
  LeftFluxes_N += (.5f * g ) * ( hL * hL );

  float RightFluxes_H, RightFluxes_U, RightFluxes_V, RightFluxes_N;
  HuDotN = (hR*uR) * edgeNormals[0] + (hR*vR*cos(M_PI*edgeCenters[1]/180.0)) * edgeNormals[1];
  // HuDotN = (hR*uR) * edgeNormals[0] + (hR*vR) * edgeNormals[1];
  RightFluxes_H =   HuDotN;
  RightFluxes_U =   HuDotN * uR;
  RightFluxes_V =   HuDotN * vR;
  // Normal Momentum flux term
  RightFluxes_N =   HuDotN * uRn;
  RightFluxes_U += (.5f * g * edgeNormals[0] ) * ( hR * hR);
  RightFluxes_V += (.5f * g * edgeNormals[1] ) * ( hR * hR* cos(M_PI*edgeCenters[1]/180.0));
  RightFluxes_N += (.5f * g ) * ( hR * hR );
  // ------------------------------------------------------------------------
  //    HLLC Flux Solver (Batten et al. 1997)
  //    "On the choice of wavespeeds for the HLLC Reimann solver"
  // ------------------------------------------------------------------------
  /*if( sL >= 0.0f) {
    out[0] = LeftFluxes_H;
    out[1] = LeftFluxes_U;
    out[2] = LeftFluxes_V;
  } else if ((sL < 0.0f) && (sStar >= 0.0f)){
    float UStarL;
    float tStar = (sL - uLn)/(sL - sStar);
    UStarL = hL*tStar;
    out[0] = LeftFluxes_H + sL*(UStarL-hL);
    out[1] = LeftFluxes_U + sL*(UStarL*(uL+ edgeNormals[0]*(sStar-uLn))-leftCellValues[1]);
    out[2] = LeftFluxes_V + sL*(UStarL*(vL+ edgeNormals[1]*(sStar-uLn))-leftCellValues[2]);
  } else if((sStar < 0.0f) && (sR >= 0.0f)){
    float UStarR;
    float tStar = (sR - uRn)/(sR - sStar);
    UStarR = hR*tStar;
    out[0] = RightFluxes_H + sR*(UStarR-hR);
    out[1] = RightFluxes_U + sR*(UStarR*(uR+ edgeNormals[0]*(sStar-uRn))-rightCellValues[1]);
    out[2] = RightFluxes_V + sR*(UStarR*(vR+ edgeNormals[1]*(sStar-uRn))-rightCellValues[2]);
  } else {
    out[0] = RightFluxes_H;
    out[1] = RightFluxes_U;
    out[2] = RightFluxes_V;
  }*/

  //------------------------------------------------------------------------
  // HLL Flux Solver
  // ------------------------------------------------------------------------
  float sLMinus = sL < 0.0f ? sL : 0.0f;
  float sRPlus = sR > 0.0f ? sR : 0.0f;
  float sRMinussL = sRPlus - sLMinus;
  sRMinussL = sRMinussL < EPS ?  EPS : sRMinussL;
  float t1 = sRPlus / sRMinussL;
  float t2 = ( -1.0 * sLMinus ) / sRMinussL;
  float t3 = ( sRPlus * sLMinus ) / sRMinussL;
  out[0] =
  ( t1 * LeftFluxes_H ) +
  ( t2 * RightFluxes_H ) +
  ( t3 * ( hR - hL ) );
  // op_printf("out %f \n", out[0]);
  out[1] =
  ( t1 * LeftFluxes_U ) +
  ( t2 * RightFluxes_U ) +
  ( t3 * ( (hR*uR) - (hL*uL) ) );

  out[2] =
  ( t1 * LeftFluxes_V ) +
  ( t2 * RightFluxes_V ) +
  ( t3 * ( (hR*vR) - (hL*vL) ) );

  // ------------------------------------------------------------------------
  //     HLLC Flux Solver (Huang et al. 2013)
  //     "Well-balanced finite volume scheme for shallow water flooding and drying
  //     over arbitrary topography"
  // ------------------------------------------------------------------------
  // float sLMinus = sL < 0.0f ? sL : 0.0f;
  // float sRPlus = sR > 0.0f ? sR : 0.0f;
  // float sRMinussL = sRPlus - sLMinus;
  // sRMinussL = sRMinussL < EPS ?  EPS : sRMinussL;
  // float t1 = sRPlus / sRMinussL;
  // float t2 = ( -1.0 * sLMinus ) / sRMinussL;
  // float t3 = ( sRPlus * sLMinus ) / sRMinussL;
  // float FStar[3];
  // FStar[0] =
  // ( t1 * LeftFluxes_H ) +
  // ( t2 * RightFluxes_H ) +
  // ( t3 * ( hR - hL ) );
  //
  //
  // FStar[1] =
  // ( t1 * LeftFluxes_N ) +
  // ( t2 * RightFluxes_N ) +
  // ( t3 * ( (hR * uRn) -
  //         (hL * uLn) ) );
  //
  // if( sL >= 0.0f) {
  //   out[0] = t1*LeftFluxes_H;
  //   out[1] = t1*LeftFluxes_U;
  //   out[2] = t1*LeftFluxes_V;
  // } else if ((sL < 0.0f) && (sStar >= 0.0f)){
  //   out[0] = FStar[0];
  //   FStar[2] = FStar[0] * uLp;
  //   out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
  //   out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  // } else if((sStar < 0.0f) && (sR >= 0.0f)){
  //   out[0] = FStar[0];
  //   FStar[2] = FStar[0] * uRp;
  //   out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
  //   out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  // } else {
  //   out[0] = t2*RightFluxes_H;
  //   out[1] = t2*RightFluxes_U;
  //   out[2] = t2*RightFluxes_V;
  // }
  out[0] *= *edgeLength;
  out[1] *= *edgeLength;
  out[2] *= *edgeLength;

  float maximum = fabs(uLn + cL);
  maximum = maximum > fabs(uLn - cL) ? maximum : fabs(uLn - cL);
  maximum = maximum > fabs(uRn + cR) ? maximum : fabs(uRn + cR);
  maximum = maximum > fabs(uRn - cR) ? maximum : fabs(uRn - cR);
  *maxEdgeEigenvalues = maximum;
}
