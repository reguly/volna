#include "event.hpp"

Timer::Timer():
  start(0.), end(INFTY), step(INFTY),
  istart(0), iend(MAXINT), istep(1),
  t(0.), iter(0), LocalIter(0), LocalTime(0.) {};

Timer::Timer( int istart_, int istep_ ,
       int iend_ , RealType start_,
       RealType step_, RealType end_ ):
  istart(istart_), istep(istep_), iend(iend_),
  start(start_), step(step_), end(end_) {};

void Timer::update( RealType dt ) {
  t+= dt;
  iter += 1;
  LocalIter += 1;
  LocalTime += dt;
}

void Timer::LocalReset() {
  LocalTime = 0.;
  LocalIter = 0;
}

void Timer::dump(TimerParams &p) {
  p.start = start;
  p.end = end;
  p.step = step;
  p.istart = istart;
  p.iend = iend;
  p.istep = istep;
}

bool Timer::happens() {
  bool result;
  result = ( t <= end && t >= start &&
	     (unsigned int)iter <= iend && (unsigned int)iter >= istart );
  result =  result &&
    ( ( iter == 0) ||
      ( (unsigned int)LocalIter == istep ) ||
      ( LocalTime >= step ) );
  return result;
}


bool Timer::IsFinished() {
  return ((unsigned int)iter >= iend) || (t >= end);
}


Event::Event():
  post_update(false) {};

Event::Event( const Timer &timer_ ):
  post_update(false), timer(timer_) {};

void Event::init( Simulation *s) {
};

void Event::finalize() {
};

