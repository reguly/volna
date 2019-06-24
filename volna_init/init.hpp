#ifndef INIT_HPP
#define INIT_HPP

#include <iostream>

#include "external/eigen2/Eigen/Array"

// #include "config.hpp"
// #include "physicalParams.hpp"
// #include "mesh.hpp"
// #include "values.hpp"
// #include "mathParser.hpp"
#include "event.hpp"

class Init : public Event {
public:
  Init( std::string &, std::string&, const Timer & );
  std::string formula;
  std::string streamName;
  virtual void dump(EventParams &p) = 0;
};

Init::Init( std::string &formula_, std::string &streamName_,
	    const Timer &timer ):
  Event( timer ), formula( formula_ ), streamName( streamName_ )
{
  post_update = false;
};

///////////////////////////////////////////////////////////////////////////////

class InitEta : public Init {
public:
  InitEta( std::string &, std::string&, const Timer & );

  virtual void dump(EventParams &p);
};

InitEta::InitEta( std::string &formula_, std::string &streamName_, const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitEta::dump(EventParams &p) {
  p.className = "InitEta";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class InitU : public Init {
public:
  InitU( std::string &, std::string&, const Timer & );

  virtual void dump(EventParams &p);
};

InitU::InitU( std::string &formula_, std::string &streamName_, const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitU::dump(EventParams &p) {
  p.className = "InitU";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class InitV : public Init {
public:
  InitV( std::string &, std::string&, const Timer & );

  virtual void dump(EventParams &p);
};

InitV::InitV( std::string &formula_, std::string &streamName_,
	      const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitV::dump(EventParams &p) {
  p.className = "InitV";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class InitBathymetry : public Init {
public:
  InitBathymetry( std::string &, std::string &, const Timer & );
  virtual void init( Simulation * );

  virtual void dump(EventParams &p);
};

InitBathymetry::InitBathymetry( std::string &formula_, std::string
				&streamName_, const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitBathymetry::init( Simulation *sim ) {
}

void InitBathymetry::dump(EventParams &p) {
  p.className = "InitBathymetry";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

class InitBathymetry_HDF : public Init {
public:
  InitBathymetry_HDF( std::string &, std::string &, const Timer & );
  virtual void init( Simulation * );

  virtual void dump(EventParams &p);
};

InitBathymetry_HDF::InitBathymetry_HDF( std::string &formula_, std::string
        &streamName_, const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitBathymetry_HDF::init( Simulation *sim ) {
}

void InitBathymetry_HDF::dump(EventParams &p) {
  p.className = "InitBathymetryHDF";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

class InitBathyRelative : public Init {
public:
  InitBathyRelative( std::string &, std::string &, const Timer & );
  virtual void init( Simulation * );

  virtual void dump(EventParams &p);
};

InitBathyRelative::InitBathyRelative( std::string &formula_, std::string
        &streamName_, const Timer &timer ):
  Init( formula_, streamName_, timer ) {};

void InitBathyRelative::init( Simulation *sim ) {
}

void InitBathyRelative::dump(EventParams &p) {
  p.className = "InitBathyRelative";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class InitBore : public Init {
private:
  RealType x0; // Location of the bore
  RealType Hl, ul, vl; // left state
  RealType S; // speed of the bore
public:
  InitBore( std::string&, std::string &, const Timer &, RealType, RealType, RealType, RealType, RealType );

  virtual void dump(EventParams &p);
};

InitBore::InitBore( std::string &formula_, std::string &streamName_, const Timer &timer,
		    RealType x0_, RealType Hl_,
		    RealType ul_, RealType
		    vl_, RealType S_ ):
  Init( formula_, streamName_, timer ),  x0(x0_), Hl(Hl_), ul(ul_), vl(vl_), S(S_) {};


void InitBore::dump(EventParams &p) {
  p.className = "InitBore";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

///////////////////////////////////////////////////////////////////////////////

class InitGaussianLandslide : public Init {
private:
  RealType A; // Amplitude of the gaussian bump
  RealType v; // velocity of the gaussian bump
  RealType lx, ly;
public:
  InitGaussianLandslide( std::string&, std::string&, const Timer &,
			 const RealType&, const RealType&, const
			 RealType&, const RealType& );
  virtual void dump(EventParams &p);
};

InitGaussianLandslide::InitGaussianLandslide( std::string &formula_, std::string &streamName_, const Timer &timer,
		    const RealType &A_, const RealType &v_,
		    const RealType &lx_, const RealType &ly_):
  Init( formula_, streamName_, timer ),  A(A_), v(v_), lx(lx_), ly(ly_) {};

void InitGaussianLandslide::dump(EventParams &p) {
  p.className = "InitGaussianLandslide";
  p.formula = formula;
  p.streamName = streamName;
  p.post_update = post_update;
}

#endif // INIT_HPP
