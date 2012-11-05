#ifndef EVENT_HPP
#define EVENT_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "config.hpp"
//#include "mesh.hpp"
//#include "geom.hpp"

using std::string;
using std::istream;

class Simulation;
class Values;
class Mesh;


struct TimerParams {
  RealType start, end, step;
  unsigned int istart, iend, istep;
  TimerParams():
    start(0.), end(0.), step(0.), istart(0), iend(0), istep(0) {}
};

struct EventParams {
  RealType location_x, location_y;
  bool post_update;
  std::string className;
  std::string formula;
  std::string streamName;
  EventParams():
    location_x(0.), location_y(0.), post_update(0), className(""), formula(""), streamName("") {}
};

class Timer {
public:
  RealType start, end, step;
  unsigned int istart, iend, istep;
  RealType t; // current time
  int iter; // current iteration nuber
  int LocalIter; // # of iteration since last time event happened
  int LocalTime; // time since last time event happened
  Timer();
  Timer( int, int, int, RealType, RealType, RealType );
  void update( RealType );
  void LocalReset();
  bool happens();
  bool IsFinished();
  void dump(TimerParams &p);
};



class Event {
public:
  Event();
  Event( const Timer& );
  virtual ~Event() {}
  Timer timer;
  bool post_update; // does event occur before or after time updating ?
  virtual void init( Simulation *s );
  virtual void finalize();
  virtual void dump(EventParams &p) = 0;
};


#endif // EVENT_HPP
