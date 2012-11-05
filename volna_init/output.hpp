#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "event.hpp"

class Output : public Event {
public:
  Output( std::string &, const Timer & );
  virtual ~Output() {}
  std::string streamName;
  virtual void dump(EventParams &p)=0;
};

///////////////////////////////////////////////////////////////////////////////
Output::Output( std::string &streamName_, const Timer &timer ):
  Event( timer ), streamName( streamName_ )
{
  post_update = true;
};

///////////////////////////////////////////////////////////////////////////////

class OutputSimulation : public Output {
public:
  OutputSimulation( std::string &, const Timer & );

  void dump(EventParams &p);
};

//////////////////////////////////////////////////////////////////////////////////

OutputSimulation::OutputSimulation( std::string & streamName_, const Timer & timer_ ):
  Output( streamName_, timer_ ) {}

void OutputSimulation::dump(EventParams &p) {
  p.className = "OutputSimulation";
  p.formula = "";
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class OutputLocation : public Output {
public:
  RealType x, y;
  OutputLocation( std::string &, const Timer &, Mesh &,
		  const RealType, const RealType );

  void dump(EventParams &p);
};

//////////////////////////////////////////////////////////////////////////////////

OutputLocation::OutputLocation( std::string &streamName_,
				const Timer &timer_,  Mesh & mesh,
				const RealType x_, const RealType y_):
  Output( streamName_, timer_ ), x(x_), y(y_) {
};

void OutputLocation::dump(EventParams &p) {
  p.className = "OutputLocation";
  p.formula = "";
  p.location_x = x;
  p.location_y = y;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////
class OutputTime : public Output {
public:
  OutputTime( std::string &, const Timer & );

  void dump(EventParams &p);
};

///////////////////////////////////////////////////////////////////////////////
OutputTime::OutputTime( std::string &streamName_,
			const Timer &timer_ ):
  Output( streamName_, timer_ ) {};

void OutputTime::dump(EventParams &p) {
  p.className = "OutputTime";
  p.formula = "";
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class OutputConservedQuantities : public Output {
public:
  OutputConservedQuantities( std::string &, const Timer & );

  void dump(EventParams &p);
};

///////////////////////////////////////////////////////////////////////////////

OutputConservedQuantities::OutputConservedQuantities( std::string &streamName_,
						      const Timer &timer_ ):
  Output( streamName_, timer_ ) {};

void OutputConservedQuantities::dump(EventParams &p) {
  p.className = "OutputConservedQuantities";
  p.formula = "";
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class OutputMaxElevation : public Output {
public:
  OutputMaxElevation( std::string &, const Timer & );
  void execute (Mesh &, Values & );
  void dump(EventParams &p);
private:
  ScalarValue CurrentMaxElevation;
};

// ///////////////////////////////////////////////////////////////////////////////

OutputMaxElevation::OutputMaxElevation( std::string &streamName_,
					const Timer &timer_ ):
  Output( streamName_, timer_) {};

void OutputMaxElevation::dump(EventParams &p) {
  p.className = "OutputMaxElevation";
  p.formula = "";
  p.streamName = streamName;
  p.post_update = post_update;
}


///////////////////////////////////////////////////////////////////////////////

class OutputLine : public Output {
public:
  RealType x_start, y_start, x_end, y_end;
  std::vector<int> ids;
  OutputLine( std::string &, const Timer &,
	      RealType, RealType,
	      RealType, RealType );

  void dump(EventParams &p);
};


///////////////////////////////////////////////////////////////////////////////

OutputLine::OutputLine( std::string &streamName_, const Timer &timer_,
			RealType x_start_, RealType y_start_,
			RealType x_end_, RealType y_end_ ):
  Output( streamName_, timer_ ), x_start(x_start_), y_start(y_start_),
  x_end(x_end_), y_end(y_end_) {
};

void OutputLine::dump(EventParams &p) {
  p.className = "OutputLine";
  p.formula = "";
  p.streamName = streamName;
  p.post_update = post_update;
}


#endif
