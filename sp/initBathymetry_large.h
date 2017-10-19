inline void initBathymetry_large(float *values, const float *cellCenter,
 const float *node0, const float *node1, const float *node2,
 const float *bathy0, const float *bathy1, const float *bathy2) {


 bool isInside = false;



 float xmin = MIN(MIN(node0[0], node1[0]), node2[0]);
 float xmax = MAX(MAX(node0[0], node1[0]), node2[0]);
 float ymin = MIN(MIN(node0[1], node1[1]), node2[1]);
 float ymax = MAX(MAX(node0[1], node1[1]), node2[1]);

 if ( ( cellCenter[0] < xmin ) || ( cellCenter[0] > xmax ) ||
  ( cellCenter[1] < ymin ) || ( cellCenter[1] > ymax ) ) {
  isInside = false;
 }else{



  float insider = 1.0f;
  float p[2] = {cellCenter[0], cellCenter[1]};


    if ( (node0[0] - node2[0]) * (node1[1] - node2[1]) - (node0[1] - node2[1]) * (node1[0] - node2[0]) > 0 ) {
     insider = (node0[0] - node2[0]) * (p[1] - node2[1]) - (node0[1] - node2[1]) * (p[0] - node2[0]);
     insider *= (node0[0] - p[0]) * (node1[1] - p[1]) - (node0[1] - p[1]) * (node1[0] - p[0]);
     insider *= (node1[0] - p[0]) * (node2[1] - p[1]) - (node1[1] - p[1]) * (node2[0] - p[0]);
    }
    else {
     insider = (node0[0] - node1[0]) * (p[1] - node1[1]) - (node0[1] - node1[1]) * (p[0] - node1[0]);
     insider *= (node0[0] - p[0]) * (node2[1] - p[1]) - (node0[1] - p[1]) * (node2[0] - p[0]);
     insider *= (node2[0] - p[0]) * (node1[1] - p[1]) - (node2[1] - p[1]) * (node1[0] - p[0]);
    }
    isInside = insider >= 0.0f;
 }


  if (isInside) {

    float a = (node1[1]-node0[1])*(*bathy2-*bathy0)-(node2[1]-node0[1])*(*bathy1-*bathy0);
    float b = -(node1[0]-node0[0])*(*bathy2-*bathy0)+(node2[0]-node0[0])*(*bathy1-*bathy0);
    float c = (node1[0]-node0[0])*(node2[1]-node0[1])-(node2[0]-node0[0])*(node1[1]-node0[1]);

    values[3] += *bathy0 - (a*(cellCenter[0]-node0[0]) + b*(cellCenter[1]-node0[1]))/c;
  }
}
