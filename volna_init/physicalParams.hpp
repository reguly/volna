#ifndef PHYSICALPARAMS_HPP
#define PHYSICALPARAMS_HPP

#include "config.hpp"

struct PhysicalParams {
  RealType g;
  RealType Cf;
  RealType coriolis;
  RealType Mn;
  PhysicalParams():
    g(9.81), Cf(0.), coriolis(0.),Mn(0.) {}
};

#endif // PHYSICALPARAMS_HPP
