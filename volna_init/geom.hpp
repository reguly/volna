#ifndef GEOM_HPP
#define GEOM_HPP

#include "external/eigen2/Eigen/Core"
#include "config.hpp"

typedef Eigen::Matrix<RealType, 3, 1> Vector;
typedef Eigen::Matrix<RealType, 3, 1> Point;
typedef Eigen::Matrix<RealType, 1, Eigen::Dynamic> ScalarValue;

class GeomValues {
public:
  GeomValues( size_t N ):
	  x( ScalarValue::Zero(N) ), 
	  y( ScalarValue::Zero(N) ),
	  z( ScalarValue::Zero(N) ) {};
  ScalarValue x;
  ScalarValue y;
  ScalarValue z;
};

RealType orient2d( const Point &pa, const Point &pb, const Point &pc ) {
  return (pa - pc).x() * (pb - pc).y() - (pa - pc).y() * (pb - pc).x();
}


#endif // GEOM_HPP
