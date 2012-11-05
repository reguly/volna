#ifndef CONFIG_H
#define CONFIG_H

#include <limits>
#include <cassert>
// #include <omp.h> // open MP header

// #define NDEBUG
// #define EIGEN_NO_DEBUG
// #define BOOST_DISABLE_ASSERTS
// #define EIGEN_VECTORIZE

#define ORDER2

// Problem parameters
unsigned int const neq = 3;  // number of scalar equations
unsigned int const ndim = 2; // ambient dimension

// CURRENTLY NOT IN USE
unsigned int const ncell = 3; // number of nodes in a mesh element
unsigned int const nfacet = 2; // number of nodes in a facet

#define RealType float

// template defining a small value (used for testing for dry cells)
// at compile time
template<typename T>
struct EPSILON {
  static T value() { return static_cast<T>(0); }
};
// EPSILON specialization for floats
template<>
struct EPSILON<float> {
  static float value() { return 1e-6; }
};
// EPSILON specialization for doubles
template<>
struct EPSILON<double> {
  static double value() { return 1e-11; }
};


RealType const EPS = EPSILON<RealType>::value();

RealType const INFTY = std::numeric_limits<RealType>::max();
RealType const MINUS_INFTY = std::numeric_limits<RealType>::min();

int const MAXINT = std::numeric_limits<int>::max();


#endif // CONFIG_H
