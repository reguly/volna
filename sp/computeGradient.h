inline void computeGradient(const float *center,
                            const float *neighbour1,
                            const float *neighbour2,
                            const float *neighbour3,
                            const float *cellCenter,
                            const float *nb1Center,
                            const float *nb2Center,
                            const float *nb3Center,
                            float *q, float *out) //OP_WRITE
{
  // Only reconstruct if the cell is not a touching the edge
  // Least-Squares Gradient Reconstruction
  if( (cellCenter[0] != nb3Center[0]) && (cellCenter[1] != nb3Center[1])){
    if(center[0]> 100.0f*EPS){
    float total, Rhs[8];
    float dh[3], dz[3],du[3], dv[3], weights[3];
    float Gram[2][2], inverse[2][2], delta[3][2];
    float x = cellCenter[0];
    float y = cellCenter[1];
    // Finding delta_x,delta_y for the neighbouring cells.
    delta[0][0] =  (nb1Center[0] - x);
    delta[0][1] =  (nb1Center[1] - y);

    delta[1][0] =  (nb2Center[0] - x);
    delta[1][1] =  (nb2Center[1] - y);

    delta[2][0] =  (nb3Center[0] - x);
    delta[2][1] =  (nb3Center[1] - y);
    // Calculating the weights coefficients based on the distance between
    // neighbouring cells and center cell.
    weights[0] = sqrt(delta[0][0] * delta[0][0] + delta[0][1] * delta[0][1]);
    weights[1] = sqrt(delta[1][0] * delta[1][0] + delta[1][1] * delta[1][1]);
    weights[2] = sqrt(delta[2][0] * delta[2][0] + delta[2][1] * delta[2][1]);
    total = weights[0] + weights[1] + weights[2];
    weights[0] = total/weights[0];
    weights[1] = total/weights[1];
    weights[2] = total/ weights[2];
    delta[0][0] *= weights[0];
    delta[0][1] *= weights[0];

    delta[1][0] *= weights[1];
    delta[1][1] *= weights[1];

    delta[2][0] *= weights[2];
    delta[2][1] *= weights[2];
    // Setting up the Gram matrix
    Gram[0][0] = ((delta[0][0]*delta[0][0]) + (delta[1][0] *delta[1][0]) + (delta[2][0] *delta[2][0]));
    Gram[0][1] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][0] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][1] = ((delta[0][1]*delta[0][1]) + (delta[1][1] *delta[1][1]) + (delta[2][1] *delta[2][1]));
    // Finding the inverse of the determinant
    float det = 1.0 / (Gram[0][0]*Gram[1][1] - Gram[0][1]*Gram[1][0]);
    inverse[0][0] = det * Gram[1][1];
    inverse[0][1] = det * (- Gram[0][1]);
    inverse[1][0] = det * (-Gram[1][0]);
    inverse[1][1] = det * Gram[0][0];
    // Setting up the RHS
    dh[0] = neighbour1[0] - center[0];
    dh[1] = neighbour2[0] - center[0];
    dh[2] = neighbour3[0] - center[0];
    dh[0] *= weights[0];
    dh[1] *= weights[1];
    dh[2] *= weights[2];

    dz[0] = neighbour1[3] - center[3];
    dz[1] = neighbour2[3] - center[3];
    dz[2] = neighbour3[3] - center[3];
    dz[0] *= weights[0];
    dz[1] *= weights[1];
    dz[2] *= weights[2];

    du[0] = neighbour1[1] - center[1];
    du[1] = neighbour2[1] - center[1];
    du[2] = neighbour3[1] - center[1];
    du[0] *= weights[0];
    du[1] *= weights[1];
    du[2] *= weights[2];

    dv[0] = neighbour1[2] - center[2];
    dv[1] = neighbour2[2] - center[2];
    dv[2] = neighbour3[2] - center[2];
    dv[0] *= weights[0];
    dv[1] *= weights[1];
    dv[2] *= weights[2];

    Rhs[0] = (delta[0][0]*dh[0]) + (delta[1][0]*dh[1]) + (delta[2][0]*dh[2]);
    Rhs[1] = (delta[0][1]*dh[0]) + (delta[1][1]*dh[1]) + (delta[2][1]*dh[2]);
    out[0] = (inverse[0][0] * Rhs[0]) + (inverse[0][1] * Rhs[1]);
    out[1] = (inverse[1][0] * Rhs[0]) + (inverse[1][1] * Rhs[1]);

    Rhs[2] = (delta[0][0]*du[0]) + (delta[1][0]*du[1]) + (delta[2][0]*du[2]);
    Rhs[3] = (delta[0][1]*du[0]) + (delta[1][1]*du[1]) + (delta[2][1]*du[2]);
    out[2] = (inverse[0][0] * Rhs[2]) + (inverse[0][1] * Rhs[3]);
    out[3] = (inverse[1][0] * Rhs[2]) + (inverse[1][1] * Rhs[3]);

    Rhs[4] = (delta[0][0]*dv[0]) + (delta[1][0]*dv[1]) + (delta[2][0]*dv[2]);
    Rhs[5] = (delta[0][1]*dv[0]) + (delta[1][1]*dv[1]) + (delta[2][1]*dv[2]);
    out[4] = (inverse[0][0] * Rhs[4]) + (inverse[0][1] * Rhs[5]);
    out[5] = (inverse[1][0] * Rhs[4]) + (inverse[1][1] * Rhs[5]);

    Rhs[6] = (delta[0][0]*dz[0]) + (delta[1][0]*dz[1]) + (delta[2][0]*dz[2]);
    Rhs[7] = (delta[0][1]*dz[0]) + (delta[1][1]*dz[1]) + (delta[2][1]*dz[2]);
    out[6] = (inverse[0][0] * Rhs[6]) + (inverse[0][1] * Rhs[7]);
    out[7] = (inverse[1][0] * Rhs[6]) + (inverse[1][1] * Rhs[7]);
  }
  }else {
    // Gradients for the edge cells are set to zero.
    out[0] = 0.0f;
    out[1] = 0.0f;
    out[2] = 0.0f;
    out[3] = 0.0f;
    out[4] = 0.0f;
    out[5] = 0.0f;
    out[6] = 0.0f;
    out[7] = 0.0f;
 }
  // Computed the local max and min values for H,U,V,Z
  // q[0] - Hmin , q[1] - Hmax
  // q[2] - Umin , q[3] - Umax
  // q[4] - Vmin , q[5] - Vmax
  // q[6] - Zmin , q[7] - Zmax
  q[0] = center[0] < neighbour1[0] ? center[0] : neighbour1[0];
  q[0] = q[0] < neighbour2[0] ? q[0] : neighbour2[0];
  q[0] = q[0] < neighbour3[0] ? q[0] : neighbour3[0];
  q[1] = center[0] > neighbour1[0] ? center[0] : neighbour1[0];
  q[1] = q[1] > neighbour2[0] ? q[1] : neighbour2[0];
  q[1] = q[1] > neighbour3[0] ? q[1] : neighbour3[0];

  q[2] = center[1] < neighbour1[1] ? center[1] : neighbour1[1];
  q[2] = q[2] < neighbour2[1] ? q[2] : neighbour2[1];
  q[2] = q[2] < neighbour3[1] ? q[2] : neighbour3[1];
  q[3] = center[1] > neighbour1[1] ? center[1] : neighbour1[1];
  q[3] = q[3] > neighbour2[1] ? q[3] : neighbour2[1];
  q[3] = q[3] > neighbour3[1] ? q[3] : neighbour3[1];

  q[4] = center[2] < neighbour1[2] ? center[2] : neighbour1[2];
  q[4] = q[4] < neighbour2[2] ? q[4] : neighbour2[2];
  q[4] = q[4] < neighbour3[2] ? q[4] : neighbour3[2];
  q[5] = center[2] > neighbour1[2] ? center[2] : neighbour1[2];
  q[5] = q[5] > neighbour2[2] ? q[5] : neighbour2[2];
  q[5] = q[5] > neighbour3[2] ? q[5] : neighbour3[2];

  q[6] = center[3] < neighbour1[3] ? center[3] : neighbour1[3];
  q[6] = q[6] < neighbour2[3] ? q[6] : neighbour2[3];
  q[6] = q[6] < neighbour3[3] ? q[6] : neighbour3[3];
  q[7] = center[3] > neighbour1[3] ? center[3] : neighbour1[3];
  q[7] = q[7] > neighbour2[3] ? q[7] : neighbour2[3];
  q[7] = q[7] > neighbour3[3] ? q[7] : neighbour3[3];
}
