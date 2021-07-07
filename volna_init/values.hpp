#ifndef VALUES_HPP
#define VALUES_HPP

#include <fstream>
#include <iostream>
#include "external/eigen2/Eigen/Core"
#include "external/eigen2/Eigen/Array"

typedef Eigen::Matrix<RealType, 1, Eigen::Dynamic> ScalarValue;

struct PhysicalParams;

enum ValueType { PHYSICAL, CONSERVATIVE };

class Values {
public:
	ValueType type;
  ScalarValue H;
  ScalarValue U;
  ScalarValue V;
  ScalarValue Zb;
  Values( size_t N ):
	  type( PHYSICAL ), H( ScalarValue::Zero(N) ), U( ScalarValue::Zero(N) ), 
	  V( ScalarValue::Zero(N) ), Zb( ScalarValue::Zero(N) ) {};
  void ReadFromFile( std::ifstream &, int );
};


#endif // VALUES_HPP
