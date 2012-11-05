#ifndef MESHIO_HPP
#define MESHIO_HPP


// Parsing functions for Mesh I/O. Currently, gmsh and VTK ascii formats
// are supported
#include <string>
#include <sstream>
#include <map>

#include "external/eigen2/Eigen/Core"
#include "external/eigen2/Eigen/StdVector"

#include "config.hpp"
#include "precomputed.hpp"
//#include "geom.hpp"
#include "meshObjects.hpp"
#include "utils.hpp"

typedef std::map<int, Point> Nodes_t;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void vtk_parse_nodes_ascii( std::istream &is, unsigned int &nb_nodes,
                            Nodes_t &nodes ) {

  unsigned int node_counter = 0;
  std::string line;

  while ( node_counter < nb_nodes ) {

    getline( is, line );
    std::istringstream isstream(line);
    std::string token;

    // get node coordinates
    double x, y, z;
    getline(isstream, token, ' ');
    if (!parse<double>(x, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }
    getline(isstream, token, ' ');
    if (!parse<double>(y, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }
    getline(isstream, token, ' ');
    if (!parse<double>(z, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }

    Point point = Point::Zero();
    point.x() = x;
    point.y() = y;
    point.z() = 0;
// #ifdef VOLNA_3D
//     point.z() = z;
// #endif

    nodes.insert(std::pair<int, Point>(node_counter, point));
    ++node_counter;
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<int N, int d>
void vtk_parse_elements_ascii( std::istream &is, unsigned int &nb_cells,
                               std::vector< Cell<N,d> > &cells ) {

  cells.reserve( nb_cells );

  unsigned int cell_counter = 0;
  std::string line;

  for ( int i = 0; i < nb_cells; ++i ) {

    getline( is, line );

    std::istringstream isstream(line);
    std::string token;

    // get cell type
    int cell_length = 0;
    getline(isstream, token, ' ');
    if (!parse<int>(cell_length, token)){
      std::cerr << "Bad input for cell type (expecting an integer): ";
      std::cerr << token << std::endl;
      exit(1);
    }

    if ( cell_length != N ) {
      std::cerr << i << "\n";
      std::cerr << "Got invalid cell length " << cell_length;
      std::cerr << " from file (expected length is ";
      std::cerr << N << ")\n";
    }

    // get cell vertices
    boost::array<int, N> vertices;

    for ( unsigned int j = 0; j < N; ++j ){

      int vertice_id = 0;
      getline(isstream, token, ' ');
      if (!parse<int>(vertice_id, token)){
        std::cerr << "Bad input for vertice index (expecting an integer): ";
        std::cerr << token << std::endl;
        exit(1);
      }

      // FIXME -- problem with VTK pixel class node ordering (which is
      // not the same as quadrangle class...)
      if ( j == 2 )
        vertices[j+1] = vertice_id;
      else if (j == 3 )
        vertices[j-1] = vertice_id;
      else
        vertices[j] = vertice_id;
    }

    cells.push_back( Cell<N,d>( vertices ) );
  }

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void gmsh_parse_nodes_ascii( std::istream &is, int &nb_nodes, Nodes_t &nodes ) {

  int node_counter = 0;
  std::string line;

  while ( node_counter < nb_nodes ) {

    getline( is, line );
    std::istringstream isstream(line);
    std::string token;
    // get node index
    int node_index;
    getline(isstream, token, ' ');
    if (!parse<int>(node_index, token)){
      std::cerr << "Bad input for node index (expecting an integer): " ;
      std::cerr << token << std::endl;
      exit(1);
    }

    // get node coordinates
    double x, y, z;
    getline(isstream, token, ' ');
    if (!parse<double>(x, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }
    getline(isstream, token, ' ');
    if (!parse<double>(y, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }
    getline(isstream, token, ' ');
    if (!parse<double>(z, token)){
      std::cerr << "Bad input for node coordinate (expecting a float/double): ";
      std::cerr << token << std::endl;
      exit(1);
    }

    Point point = Point::Zero();
    point.x() = x;
    point.y() = y;
    point.z() = 0;
// #ifdef VOLNA_3D
//     point.z() = z;
// #endif

    nodes.insert(std::pair<int, Point>(node_index, point));
    ++node_counter;
  }
}

///////////////////////////////////////////////////////////////////////////////

void gmsh_parse_nodes_binary( std::istream &is, int &nb_nodes,
                              Nodes_t &nodes ) {
  for (int i = 0; i < nb_nodes; ++i){
    int node_index;
    double x, y, z;
    is.read( (char*)&node_index, 4 );
    is.read( (char*)&x, 8 );
    is.read( (char*)&y, 8 );
    is.read( (char*)&z, 8 );
    
    Point point = Point::Zero();
    point.x() = x;
    point.y() = y;
    //point.z() = z;
    point.z() = 0;
    nodes.insert( std::pair<int,Point>( node_index, point ) );
  }
}

///////////////////////////////////////////////////////////////////////////////

template<int N, int d>
Cell<N,d> gmsh_parse_element_binary ( std::istream &is, int &nb_tags ){

  int tmp = 0;
  is.read( (char*)&tmp, sizeof(int) ); // element number
  // read tags
  for ( int i = 0; i < nb_tags; ++i ) {
    if ( i==0 ) {
      int elementary = 0;
      is.read( (char*)&elementary, sizeof(int) );
    }
    else if ( i==1 ) {
      int physical = 0;
      is.read( (char*)&physical, sizeof(int) );
    }
    else if ( i==2 ) {
      int partition = 0;
      is.read( (char*)&partition, sizeof(int) );
    }
  }
  // read vertices
  boost::array<int, N> vertices;
  for ( int i = 0; i < N; ++i ) {
    int vertice;
    is.read( (char*)&vertice, sizeof(int) );
    vertices[i] = vertice;
  }

  return Cell<N,d>( vertices );
}

///////////////////////////////////////////////////////////////////////////////

template<int N, int d>
Facet<N,d> gmsh_parse_facet_binary ( std::istream &is, int &nb_tags ){

  int tmp = 0;
  is.read( (char*)&tmp, sizeof(int) ); // element number
  // read tags
  for ( int i = 0; i < nb_tags; ++i ) {
    if ( i==0 ) {
      int elementary = 0;
      is.read( (char*)&elementary, sizeof(int) );
    }
    else if ( i==1 ) {
      int physical = 0;
      is.read( (char*)&physical, sizeof(int) );
    }
    else if ( i==2 ) {
      int partition = 0;
      is.read( (char*)&partition, sizeof(int) );
    }
  }
  // read vertices
  boost::array<int, N> vertices;
  for ( int i = 0; i < N; ++i ) {
    int vertice;
    is.read( (char*)&vertice, sizeof(int) );
    vertices[i] = vertice;
  }

  return Facet<N,d>( vertices );
}

///////////////////////////////////////////////////////////////////////////////

template<size_t N, int d>
void gmsh_parse_elements_ascii( std::istream &is, int &nb_cells,
                                std::vector< Facet<N-1,d-1> > &boundary_facets,
                                std::vector< Cell<N,d> > &cells ) {

  cells.reserve( nb_cells );
  std::string line;

  for ( unsigned int i = 0; i < nb_cells; ++i ) {

    getline( is, line );

    std::istringstream isstream( line );
    std::string token;

    // get cell type

    // get cell vertices
    boost::array<int, N> vertices;

    for ( unsigned int j = 0; j<N; ++j ) {

      int vertice_id = 0;
      getline( isstream, token, ' ' );
      if (!parse<int>( vertice_id, token )) {
        std::cerr << "Bad input for vertice index "
                  << "(expecting an integer):  "
                  << token << std::endl;
        exit(1);
      }

      vertices[j] = vertice_id;
    }

    cells.push_back( Cell<N,d>( vertices ) );

  }

}

///////////////////////////////////////////////////////////////////////////////

template<size_t N, int d>
void gmsh_parse_elements_binary ( std::istream &is, int &nb_elements,
                                  std::vector< Facet<N-1,d-1> > &boundary_faces,
                                  std::vector< Cell<N,d> > &elements) {

  int element_counter = 0;
  while ( element_counter < nb_elements ){
    // elements are grouped by their type
    int cell_type = 0;
    is.read( (char*)&cell_type, sizeof(int) );
    int nb_cells = 0;
    is.read( (char*)&nb_cells, sizeof(int) );
    int nb_tags = 0;
    is.read( (char*)&nb_tags, sizeof(int) );
    // There are only 3 valid values of cell type: 15 (point), and
    // the gmsh type for faces and elements
    if ( cell_type == 15 ) { // point -- do nothing
      for ( int i = 0; i < nb_cells; ++i ) {
        Cell<1,0> cell = gmsh_parse_element_binary<1,0>( is, nb_tags );
        ++element_counter;
      }
    }

    else if ( cell_type == 1 && d == 3) { // line in dim 3 -- do nothing
      for ( int i = 0; i < nb_cells; ++i ) {
        Cell<2,1> cell = gmsh_parse_element_binary<2,1>( is, nb_tags );
        ++element_counter;
      }
    }

    else if ( cell_type == GMSH_TYPE<N,d>::value() ) {
      for ( int i = 0; i < nb_cells; ++i ) {
        Cell<N,d> cell = gmsh_parse_element_binary<N,d>( is, nb_tags );
        elements.push_back( cell );
        ++element_counter;
      }
    }

    else if ( cell_type == GMSH_FACE_TYPE<N,d>::value() ) {
      for ( int i = 0; i < nb_cells; ++i ) {
        Facet<N-1,d-1> facet = gmsh_parse_facet_binary<N-1,d-1>( is, nb_tags );
        boundary_faces.push_back( facet );
        ++element_counter;
      }
    }

    else {
      std::cerr << "Error: invalid element, of type ";
      std::cerr << cell_type << ", at line "  << element_counter ;
      std::cerr << " (valid types are 15, 1, " << GMSH_TYPE<N,d>::value() << ", and ";
      std::cerr <<  GMSH_FACE_TYPE<N,d>::value() << " )\n";
      std::abort();
    }
  }
}

#endif // MESHIO_HPP
