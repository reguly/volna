inline void initBathymetry_large(float *values, float *cellCenter,
	float *node0, float *node1, float *node2,
	float *bathy0, float *bathy1, float *bathy2) {

	//First, check whether the cell is in the large triangle
	bool isInside = false;

  	// First, check if the point is in the bounding box of the triangle
  	// vertices (else, the algorithm is not nearly robust enough)
	float xmin = MIN(MIN(node0[0], node1[0]), node2[0]);
	float xmax = MAX(MAX(node0[0], node1[0]), node2[0]);
	float ymin = MIN(MIN(node0[1], node1[1]), node2[1]);
	float ymax = MAX(MAX(node0[1], node1[1]), node2[1]);

	if ( ( cellCenter[0] < xmin ) || ( cellCenter[0] > xmax ) ||
		( cellCenter[1] < ymin ) || ( cellCenter[1] > ymax ) ) {
		isInside = false;
	}else{
	    // Case where the point is in the bounding box. Here, if abc is not
	    // Check if the Triangle vertices are clockwise or
	    // counter-clockwise
		float insider = 1.0f;
		float p[2] = {cellCenter[0], cellCenter[1]};

#define ORIENT2D(pA, pB, pC) (pA[0] - pC[0]) * (pB[1] - pC[1]) - (pA[1] - pC[1]) * (pB[0] - pC[0])
    if ( ORIENT2D(node0, node1, node2) > 0 ) {  // counter clockwise
    	insider =  ORIENT2D( node0, p, node2);
    	insider *= ORIENT2D( node0, node1, p);
    	insider *= ORIENT2D( node1, node2, p);
    }
    else {      // clockwise
    	insider =  ORIENT2D( node0, p, node1);
    	insider *= ORIENT2D( node0, node2, p);
    	insider *= ORIENT2D( node2, node1, p);
    }
    isInside = insider >= 0.0f;
	}

  //If it's inside, then calculate the bathymetry value, with a simple linear interpolation
  if (isInside) {
    //get the normal vector
    float a =  (node1[1]-node0[1])*(*bathy2-*bathy0)-(node2[1]-node0[1])*(*bathy1-*bathy0);
    float b = -(node1[0]-node0[0])*(*bathy2-*bathy0)+(node2[0]-node0[0])*(*bathy1-*bathy0);
    float c =  (node1[0]-node0[0])*(node2[1]-node0[1])-(node2[0]-node0[0])*(node1[1]-node0[1]);
    //get z so that the vector is orthogonal to the normal
    values[3] = *bathy0 - (a*(cellCenter[0]-node0[0]) + b*(cellCenter[1]-node0[1]))/c;
  }
}