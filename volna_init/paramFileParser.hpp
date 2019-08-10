#ifndef PARAMFILEPARSER_HPP
#define PARAMFILEPARSER_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/dynamic/switch.hpp>
#include <boost/spirit/iterator/file_iterator.hpp>

#include "simulation.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using boost::shared_ptr;

namespace spirit = boost::spirit;

///////////////////////////////////////////////////////////////////////////////

struct skip_grammar : spirit::grammar<skip_grammar>
{
  Simulation simulation;
  template <typename ScannerT>
  struct definition {
    definition( skip_grammar const & )
    {
      using spirit::blank_p;
      using spirit::anychar_p;

      skip
        = blank_p | '#' >> *(anychar_p - '\n');
    }

    typedef typename spirit::rule<ScannerT> rule_t;

    rule_t const & start() const { return skip; }

  private:
    rule_t skip;

  };
};

///////////////////////////////////////////////////////////////////////////////
// semantic actions for the parser

struct assign_string_v
{
public:
  assign_string_v( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.V = string(first, last);
  }
private:
  Simulation &simulation;
};

struct assign_string_u
{
public:
  assign_string_u( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.U = string(first, last);
  }
private:
  Simulation &simulation;
};

struct assign_string_eta
{
public:
  assign_string_eta( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.eta = string(first, last);
  }
private:
  Simulation &simulation;
};

struct assign_string_bathymetry
{
public:
  assign_string_bathymetry( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.bathymetry = string(first, last);
  }
private:
  Simulation &simulation;
};

struct assign_string_bathymetry_hdf
{
public:
  assign_string_bathymetry_hdf( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.bathymetry_hdf = string(first, last);
  }
private:
  Simulation &simulation;
};

struct assign_string_bathyrelative
{
public:
  assign_string_bathyrelative( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    simulation.InitFormulas.bathyrelative = string(first, last);
  }
private:
  Simulation &simulation;
};


struct add_event
{
public:
  add_event( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
    const Timer timer = simulation.EventTimer;

    if ( simulation.EventName == "OutputSimulation" ) {
      boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputSimulation( simulation.EventFilename, timer ) );
      simulation.events.push_back( ptr );
    }

    else if ( simulation.EventName == "OutputTime" ) {
      boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputTime( simulation.EventFilename, timer ) );
      simulation.events.push_back( ptr );
    }

    else if ( simulation.EventName == "OutputConservedQuantities" ) {
      boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputConservedQuantities( simulation.EventFilename, timer ) );
      simulation.events.push_back( ptr );
    }


    else if (simulation.EventName == "OutputLocation" ) {
      RealType x = simulation.LocationPoint.x();
      RealType y = simulation.LocationPoint.y();
	boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputLocation(simulation.EventFilename,timer, 
      			     simulation.mesh, x, y) );
      simulation.events.push_back( ptr );
    }

    else if (simulation.EventName == "OutputMaxElevation" ) {
      boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputMaxElevation(simulation.EventFilename,timer));
	simulation.events.push_back( ptr );
    }
    else if (simulation.EventName == "OutputMaxSpeed" ) {
      boost::shared_ptr<Event> ptr = boost::shared_ptr<Event>
        ( new OutputMaxSpeed(simulation.EventFilename,timer));
	simulation.events.push_back( ptr );
    }
    
    // reset the timer
    simulation.EventTimer = Timer();
  }
private:
  Simulation &simulation;
};

struct add_init_event
{
public:
  add_init_event( Simulation &s ):
    simulation( s ) {}
  template <typename ItT>
  void operator() ( ItT first, ItT last ) const
  {
  
    const Timer timer = simulation.EventTimer;
    
    // add the event to the list of events
  
    boost::shared_ptr<Event> ptr;

    if ( simulation.InitVar == "Eta" ) {
      ptr = boost::shared_ptr<Event>
	( new InitEta(simulation.InitFormula, 
		      simulation.InitFilename, timer) );
    }
    else if ( simulation.InitVar == "U" ) {
      ptr = boost::shared_ptr<Event>
	( new InitU(simulation.InitFormula,
		    simulation.InitFilename, timer) );
    }
    else if ( simulation.InitVar == "V" ) {
      ptr = boost::shared_ptr<Event>
	( new InitV(simulation.InitFormula,
		    simulation.InitFilename, timer) );
    }

    else if ( simulation.InitVar == "Bathymetry" ) {
      ptr = boost::shared_ptr<Event>
	( new InitBathymetry( simulation.InitFormula, 
			      simulation.InitFilename, timer ) );
    }

    else if ( simulation.InitVar == "BathyHDF" ) {
      ptr = boost::shared_ptr<Event>
  ( new InitBathymetry_HDF( simulation.InitFormula, 
            simulation.InitFilename, timer ) );
    }

    else if ( simulation.InitVar == "BathyRelative" ) {
      ptr = boost::shared_ptr<Event>
  ( new InitBathyRelative( simulation.InitFormula, 
            simulation.InitFilename, timer ) );
    }

    else if ( simulation.InitVar == "Bore" ) {
      RealType x0 = simulation.bore_params.x0;
      RealType Hl = simulation.bore_params.Hl;
      RealType ul = simulation.bore_params.ul;
      RealType vl = simulation.bore_params.vl;
      RealType S = simulation.bore_params.S;
      ptr = boost::shared_ptr<Event>
      	( new InitBore(simulation.InitFormula, 
		       simulation.InitFilename, timer, x0, Hl, ul, vl, S) );
    }

    else if ( simulation.InitVar == "GaussianLandslide" ) {
      std::cerr << "landslide\n";
      RealType A = simulation.gaussian_landslide_params.A;
      RealType v = simulation.gaussian_landslide_params.v;
      RealType lx = simulation.gaussian_landslide_params.lx;
      RealType ly = simulation.gaussian_landslide_params.ly;
      ptr = boost::shared_ptr<Event>
      	( new InitGaussianLandslide(simulation.InitFormula, 
		       simulation.InitFilename, timer, A, v, lx, ly) );
    }


    simulation.events.push_back( ptr );
  
    // reset the timer
    simulation.EventTimer = Timer();
  }
private:
  Simulation &simulation;
};

///////////////////////////////////////////////////////////////////////////////

struct volna_grammar : public spirit::grammar<volna_grammar>
{

  Simulation &sim;

  volna_grammar( Simulation & s )
    : sim( s ) {}

  template <typename ScannerT>
  struct definition
  {
    definition( volna_grammar const & self)
    {
      using spirit::ch_p;
      using spirit::anychar_p;
      using spirit::str_p;
      using spirit::alpha_p;
      using spirit::alnum_p;
      using spirit::space_p;
      using spirit::blank_p;
      using spirit::real_p;
      using spirit::lexeme_d;
      using spirit::uint_p;
      using spirit::assign_a;
      using spirit::nothing_p;
      using spirit::switch_p;
      using spirit::case_p;
      using spirit::if_p;

      // a sequence of simulation objects separated by (optional) blank lines)
      simulation_objects_list
        = *(blank_p) >> *(simulation_object >> +(ch_p('\n')));

      simulation_object
        = time | physical_params | numerical_params
        | mesh | event | initial_value | boundary_condition;

      numerical_params
        = str_p("NumericalParams") >> ch_p('{')
                                   >> *(num_param_options)
                                   >> ch_p('}');

      num_param_options
        = ( str_p("cfl") >> ch_p('=')
            >> real_p[ assign_a( self.sim.CFL ) ] );

      physical_params
        = str_p("PhysicalParams") >> ch_p('{')
                                  >> *(param_options)
                                  >> ch_p('}');
      param_options
        = ( str_p("g") >> ch_p('=')
            >> real_p[ assign_a( self.sim.Params.g ) ] );

      time
        = str_p("Time") >> ch_p('{')
                        >> *(time_options)
                        >> ch_p('}');

      time_options
        = ( (str_p("end") >> ch_p("=") >> real_p
             [ assign_a( self.sim.FinalTime ) ] ) |
            (str_p("dtmax") >> ch_p("=") >> real_p
             [ assign_a( self.sim.Dtmax ) ] )
            );
      
      mesh
        = str_p("Mesh")
        >> ( (ch_p('"')
              >> meshfilename [ assign_a( self.sim.MeshFileName ) ]
              >> ch_p('"'))
             |
             ( str_p("Rectangle")
               >> real_p[assign_a(self.sim.mesh.xmin)]
               >> real_p[assign_a(self.sim.mesh.xmax)]
               >> real_p[assign_a(self.sim.mesh.ymin)]
               >> real_p[assign_a(self.sim.mesh.ymax)]
               >> uint_p[assign_a(self.sim.mesh.nx)]
               >> uint_p[assign_a(self.sim.mesh.ny)] ));

      meshfilename = lexeme_d[*(alnum_p|ch_p('/')|ch_p('_')|ch_p('-'))
			      >> str_p(".msh") ];
      
      event = output_generic | output_location | init;

      output_generic
	= ( ( (str_p("OutputSimulation")|
	       str_p("OutputTime")|
	       str_p("OutputConservedQuantities")|
	       str_p("OutputMaxElevation")|
	       str_p("OutputMaxSpeed"))
	     [ assign_a( self.sim.EventName ) ] )
	    >> ch_p('{') 
	    >> *(timer_option) 
	    >> ch_p('}')
	    >> ch_p('"')
	    >> event_filename[ assign_a( self.sim.EventFilename ) ]
	    >> ch_p('"')
	    >> ( !(ch_p('{') 
		   >> nothing_p
		   >> ch_p('}')
		   )))
	[ add_event(self.sim) ];

      output_location 
	= ( (str_p("OutputLocation")
	     [ assign_a( self.sim.EventName ) ] 
	     >> ch_p('{') 
	     >> *(timer_option) 
	     >> ch_p('}')
	     >> ch_p('"')
	     >> event_filename[ assign_a( self.sim.EventFilename ) ]
	     >> ch_p('"')
	     >> ( !(ch_p('{') 
		    >> output_location_options
		    >> ch_p('}')
		    ))))
	[ add_event(self.sim) ];

      output_location_options 
	= str_p("x") 
	>> ch_p("=") 
	>> real_p[assign_a(self.sim.LocationPoint.x())]
	>> str_p("y")
	>> ch_p("=")
	>> real_p[assign_a(self.sim.LocationPoint.y())];

      timer_option =
        ( (str_p("istart") >> ch_p("=") >> uint_p[assign_a(self.sim.EventTimer.istart)]) |
          (str_p("istep") >> ch_p("=") >> uint_p[assign_a(self.sim.EventTimer.istep)])|
          (str_p("iend") >> ch_p("=") >> uint_p[assign_a(self.sim.EventTimer.iend)])|
          (str_p("start") >> ch_p("=") >> real_p[assign_a(self.sim.EventTimer.start)])|
          (str_p("step") >> ch_p("=") >> real_p[assign_a(self.sim.EventTimer.step)])|
          (str_p("end") >> ch_p("=") >> real_p[assign_a(self.sim.EventTimer.end)])
	  );
      
      init = init_var | init_bore | init_gaussian_landslide;

      init_var
= ( str_p("Init")
      >> ( ( 
        ch_p('{') 
        >> *(timer_option) 
        >> ch_p('}')
        >>
        (str_p("Bathymetry") |str_p("BathyRelative")|str_p("Eta")|str_p("U")|str_p("V")) 
         [ assign_a( self.sim.InitVar ) ] ) | ( (str_p("BathyHDF")) [ assign_a( self.sim.InitVar ) ] ) )
      >> 
      (
       ( ch_p('{') >> math_expression
         [ assign_a( self.sim.InitFormula ) ]
         >> ch_p('}') ) |
       ( ch_p('"') >> data_filename
        [ assign_a( self.sim.InitFilename ) ]
        >> ch_p('"') )
        )
      )
  [ add_init_event( self.sim ) ];


      init_bore
	= ( str_p("Init")
	    >> ch_p('{') 
	    >> *(timer_option) 
	    >> ch_p('}')
	    >>
	    (str_p("Bore"))
	    [ assign_a( self.sim.InitVar ) ]
	    >> ch_p('{') 
	    >> *(init_bore_option)
	    [ assign_a( self.sim.InitFormula ) ]
	    >> ch_p('}') )
	[ add_init_event( self.sim ) ];

      init_gaussian_landslide
	= ( str_p("Init")
	    >> ch_p('{') 
	    >> *(timer_option) 
	    >> ch_p('}')
	    >>
	    (str_p("GaussianLandslide"))
	    [ assign_a( self.sim.InitVar ) ]
	    >> ch_p('{') 
	    >> *(init_gaussian_landslide_option)
	    [ assign_a( self.sim.InitFormula ) ]
	    >> ch_p('}') )
	[ add_init_event( self.sim ) ]; 

      init_gaussian_landslide_option = 
	( (str_p("A") >> ch_p("=") >> real_p[assign_a(self.sim.gaussian_landslide_params.A)]) |
          (str_p("V") >> ch_p("=") >> real_p[assign_a(self.sim.gaussian_landslide_params.v)])|
	  (str_p("lx") >> ch_p("=") >> real_p[assign_a(self.sim.gaussian_landslide_params.lx)])|
	  (str_p("ly") >> ch_p("=") >> real_p[assign_a(self.sim.gaussian_landslide_params.ly)])
	  );

      init_bore_option = 
	( (str_p("x0") >> ch_p("=") >> real_p[assign_a(self.sim.bore_params.x0)]) |
          (str_p("Hl") >> ch_p("=") >> real_p[assign_a(self.sim.bore_params.Hl)])|
	  (str_p("ul") >> ch_p("=") >> real_p[assign_a(self.sim.bore_params.ul)])|
	  (str_p("S") >> ch_p("=") >> real_p[assign_a(self.sim.bore_params.S)])
	  );
      
      event_filename =
	lexeme_d[*(alnum_p|ch_p('/')|ch_p('_')|ch_p('.')|ch_p('%'))];

      data_filename = 
	lexeme_d[*(alnum_p|ch_p('/')|ch_p('_')|ch_p('.')|ch_p('%')|ch_p('-'))];

      initial_value
        = str_p("Init") >>
        ( ( str_p("eta") >> ch_p("{")
            >> math_expression[ assign_string_eta( self.sim ) ]
            ) |
          ( str_p("U") >> ch_p("{")
            >> math_expression[ assign_string_u( self.sim ) ]
            ) |
          ( str_p("V") >> ch_p("{")
            >> math_expression[ assign_string_v( self.sim ) ]
            ) |
          ( str_p("BathyRelative") >> ch_p("{")
            >> math_expression[ assign_string_bathyrelative( self.sim ) ]
            ) |
          ( str_p("Bathymetry") >> ch_p("{")
            >> math_expression[ assign_string_bathymetry( self.sim ) ]
            ) |
          ( str_p("BathyHDF") >> ch_p("{")
            >> math_expression[ assign_string_bathymetry_hdf( self.sim ) ]
            )) >> ch_p('}');

      math_expression = lexeme_d[*(anychar_p - str_p('}'))];

      boundary_condition = 
	(
	 str_p("BoundaryCondition") 
	 >> uint_p[assign_a(self.sim.BoundaryRegionNumber)]
	 >> (
	     str_p("SubCriticalInflowEta") |
	     str_p("SubCriticalOutflowEta") |
	     str_p("SubcriticalInFlowSpeed") |
	     str_p("SubcriticalOutflowSpeed")
	     )
	 >> math_expression
	 );


    }

    typedef typename spirit::rule<ScannerT> rule_t;

    rule_t const & start() const { return simulation_objects_list; }

  private:
    rule_t simulation_objects_list, simulation_object;
    rule_t param_options, num_param_options;
    rule_t time, time_options;
    rule_t mesh, initial_value, event, bathymetry, bathyrelative, bathymetry_hdf, physical_params;
    rule_t boundary_condition;
    rule_t numerical_params;
    rule_t meshfilename, event_filename;
    rule_t event_name, timer_option, event_options;
    rule_t math_expression;
    rule_t output, output_option_list, filename;
    rule_t output_generic, output_location;
    rule_t output_location_options;
    rule_t data_filename;
    rule_t init, init_var, init_bore, init_bore_option;
    rule_t init_gaussian_landslide, init_gaussian_landslide_option;

  };
};

///////////////////////////////////////////////////////////////////////////////

class ParamFileController {
public:
  ParamFileController() {};

  template <typename ItT>
  bool parse(ItT first, ItT last, Simulation & s )
  {
    typedef spirit::parse_info<ItT> parse_info_t;

    skip_grammar skipper;

    parse_info_t info = spirit::parse( first, last,
				       volna_grammar( s ),
				       skipper );
    std::cout << "stopped at: \n" << std::string(info.stop, last) << std::endl;

    if (info.full) {
      cerr << "Parsing succeeded." << endl;
    }

    else
      cerr << "Parsing failed." << endl;


    return (info.full);

  }


};


#endif // PARAMFILEPARSER_HPP
