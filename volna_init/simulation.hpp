
#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <boost/shared_ptr.hpp>

#include "config.hpp"
#include "values.hpp"
#include "mesh.hpp"
#include "physicalParams.hpp"
#include "init.hpp"
#include "output.hpp"

typedef boost::shared_ptr<Event> EventPtr;
typedef std::vector<EventPtr> EventCollection;
typedef EventCollection::iterator EventIterator;

struct BoreParams {
  RealType x0, Hl, ul, vl;
  RealType S;
  BoreParams():
    x0(0.), Hl(1.), ul(0.), vl(0.), S(0.) {}
};

struct GaussianLandslideParams {
  RealType A, v, lx, ly;
  GaussianLandslideParams():
    A(0.), v(1.), lx(0.), ly(0.) {}
};

struct InitExpression {
  std::string eta;
  std::string U;
  std::string V;
  std::string bathymetry;
  std::string bathyrelative;
  InitExpression():
    eta(""), U(""), V(""), bathymetry(""), bathyrelative("") {}
};

class Simulation {
public:
  BoreParams bore_params;
  GaussianLandslideParams gaussian_landslide_params;
  PhysicalParams Params;
  Timer EventTimer;
  std::string EventName, EventFilename;
  unsigned int BoundaryRegionNumber;
  Timer timer;
  RealType FinalTime, Dtmax;
  RealType CFL;
  std::string MeshFileName;
  Point LocationPoint;
  Mesh mesh;
  // mathematical formulas for Initial value and bathymetry
  InitExpression InitFormulas;
  EventCollection events;
  Values CellValues;
  std::string InitVar;
  std::string InitFilename;
  std::string InitFormula;
  ScalarValue Bathymetry;
  ScalarValue BathyRelative;
  //yac::controller mathParser;
  Simulation();
  Simulation( istream &is );
  void init();
private:

};

Simulation::Simulation():
  FinalTime(INFTY), Dtmax(INFTY), CFL(.9), CellValues(10000),
  InitFilename(""), InitFormula("") {};

void Simulation::init() {

  std::ifstream ifs( MeshFileName.c_str() );

  if (ifs)
    mesh.readGmsh( ifs );

  else if ( mesh.nx!= 0 && mesh.ny!=0 )
    mesh.InitRectangle();

  else {
        std::cerr << "Could not initialize the mesh, aborting\n";
        std::abort();
      }

  // mesh.RCMRenumbering();
  mesh.ComputeConnectivity();

  // std::ofstream stream3( "bandwith.ppm" );
  // mesh.WriteMeshBandwith( stream3 );
  // stream3.close();

  //mesh.RCMRenumbering();
  //mesh.ComputeConnectivity();
  mesh.LegacyInterface();
  mesh.ComputeGeometricQuantities();
  mesh.ComputeGradientInterpolator();

  // std::ofstream stream2( "bandwith_rcm.ppm" );
  // mesh.WriteMeshBandwith( stream2 );
  // stream2.close();

  std::vector<RealType> X, Y;
  for ( int i =0; i < mesh.NVolumes; ++i ) {
    X.push_back( mesh.CellCenters.x(i) );
    Y.push_back( mesh.CellCenters.y(i) );
  }

  mesh.mathParser.updateReservedVariables( X, Y );

	CellValues.H.resize( mesh.NVolumes );
	CellValues.U.resize( mesh.NVolumes );
	CellValues.V.resize( mesh.NVolumes );
	CellValues.Zb.resize( mesh.NVolumes );

	for ( int i = 0; i<mesh.NVolumes; ++i ) {
		CellValues.H(i) = 0.0f;
		CellValues.U(i) = 0.0f;
		CellValues.V(i) = 0.0f;
		CellValues.Zb(i) = 0.0f;
	}

	if (InitFormulas.eta.begin() != InitFormulas.eta.end() ||
			InitFormulas.U.begin() != InitFormulas.U.end() ||
			InitFormulas.V.begin() != InitFormulas.V.end() ||
      InitFormulas.bathyrelative.begin() != InitFormulas.bathyrelative.end() ||
			InitFormulas.bathymetry.begin() != InitFormulas.bathymetry.end()) {
		printf("Unsupported method of specifiying initial values: please use Init events\n");
		exit(-1);
	}
}

#endif // SIMULATION_HPP
