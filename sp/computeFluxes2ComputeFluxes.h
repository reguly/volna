inline void computeFluxes2ComputeFluxes(const float *edgeLength, const float *edgeNormals,
                                const float* leftCellValues,
                                const float* rightCellValues,
                                float *out, //OP_WRITE
                                float *maxEdgeEigenvalues) //OP_WRITE
{
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