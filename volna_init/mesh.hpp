#ifndef MESH_HPP
#define MESH_HPP

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "external/eigen2/Eigen/StdVector"
#include "external/eigen2/Eigen/LU"
#include "external/eigen2/Eigen/Core"

#include <queue>

#include "config.hpp"
#include "geom.hpp"
#include "meshObjects.hpp"
#include "meshIo.hpp"
#include "mathParser.hpp"

inline bool PairCompare( const std::pair<int,Face> &elem1,
                         const std::pair<int,Face> &elem2 )
{
  return ( elem1.second < elem2.second );
}

typedef std::vector<Triangle> Cells_t;



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class Mesh {
public:
  RealType hmin, hmax; // min/max length of faces
  int nx, ny;
  GeomValues CellCenters;
  GeomValues FacetCenters;
  GeomValues FacetNormals;
  ScalarValue CellVolumes;
  ScalarValue FacetVolumes;
  yac::controller mathParser;
  // Matrix for Least Square gradient reconstruction
  typedef Eigen::Matrix<RealType,ndim,ndim> LSQMatrix;
  Mesh();
  Nodes_t Nodes;
  Cells_t Cells;
  Cells_t::const_iterator CellsIter;
  std::vector<Face> Facets;
  std::vector<Face> BoundaryFaces;
  std::vector< std::pair<int, Face> > tracked_faces_;
  std::vector<LSQMatrix> GradientInterpolator;
  Eigen::VectorXi CellsVerticesInterpolator;
  //std::vector<WeightsAtPoint> W;
  RealType xmin, xmax, ymin, ymax;
  int InitRectangle();
  int NFaces, NVolumes, NPoints;
  int readGmsh( std::istream &);
  int readOceanMesh( std::istream &);
  void RCMRenumbering();
  void WriteMeshBandwith( std::ofstream & );
  void ComputeConnectivity();
  void LegacyInterface();
  void ComputeGeometricQuantities();
  // void ComputeGradientInterpolator();
};

int Mesh::InitRectangle() {

#ifdef VOLNA_3D
  std::cerr << "Error : Initializing a 2D mesh ";
  std::cerr << "while Volna is compiled for 3D. Aborting\n";
  std::abort();
#endif

  hmin = INFTY;
  hmax = MINUS_INFTY;
  NPoints = (nx+1)*(ny+1);
  NVolumes = 2*(nx)*(ny);

   int counter = 1;
  for ( int j = 0; j <= ny; ++j) {
    for ( int i=0; i <= nx; ++i) {

      Point point = Point::Zero();
      point.x() = xmin + i* (xmax-xmin)/(RealType) nx;
      point.y() = ymin + j* (ymax-ymin)/(RealType) ny;

      Nodes.insert( std::pair<int,Point> (counter, point) );

      ++counter;
    }
  }

  for (  int i = 0; i < nx; ++i ) {
    for (  int j = 0; j < ny; ++j ) {

      // Vertices numbers start at 1
      boost::array<int,3> vertices =
        {{ j * (nx+1) + i + 1, 
	   (j+1) * (nx+1) + i+ 1, 
	   (j+1)* (nx+1) + i + 2    }};
      Triangle t(vertices);

      Cells.push_back(t);

      boost::array<int,3> vertices2 =
        { { (j*(nx+1)) + i+1, (j*(nx+1)) + i+2, ((j+1)*(nx+1))+i+2 } };
      Triangle t2(vertices2);

      Cells.push_back(t2);

    }
  }

  return 0;
}


Mesh::Mesh():
  hmin( INFTY ), hmax( MINUS_INFTY ), nx(0), ny(0),
  CellCenters(GeomValues::GeomValues(10000)),
  FacetCenters(GeomValues::GeomValues(10000)),
  FacetNormals(GeomValues::GeomValues(10000)) {};

int Mesh::readGmsh( std::istream &is ) {

  int filetype = 0;             // 0 --> ascii, 1 --> binary
  int datasize = sizeof( double ); // the gmsh standard for now

  if (!is) {
    std::cerr << "Error while reading mesh from stream. Aborting.\n";
    std::abort();
  }

  std::string line = "";                // string holding the current line

  while ( getline(is ,line) ) { // Iterate through stream, line by line

    if ( line.substr(0, 2) == std::string("//") ) // ignore comments
      {}

    // Read file header
    else if (line.substr(0, 11) == std::string("$MeshFormat")){

      std::cerr << "Parsing header..." << std::endl;
      getline(is, line);

      std::istringstream isstream(line);
      std::string token;
      // get file version -- for now, do nothing
      getline(isstream, token, ' ');

      // get file type (ascii or binary)
      getline(isstream, token, ' ');
      if (!parse<int>(filetype, token)){
        std::cerr << "  bad input for gmsh filetype (expecting an integer): " ;
        std::cerr << token << std::endl;
        exit(1);
      }
      if ( !(filetype <= 1) ) {
        std::cerr << "  invalid value for gmsh filetype (expecting 0 or 1): ";
        std::cerr << filetype << std::endl;
        exit(1);
      }
      if ( filetype == 0){
        std::cerr << "  gmsh filetype: ascii" << std::endl;
      }
      else {
        std::cerr << "  gmsh filetype: binary" << std::endl;
      }

      std::cerr << "done.\n" << std::endl;
    }

    // Read mesh nodes
    else if ( line.substr(0, 6) == std::string("$Nodes") ) {
      std::cerr << "Parsing nodes..." << std::endl;
      getline( is, line );
      // Get node numbers
      int nb_nodes;
      if( parse<int>(nb_nodes, line) )
        {
          std::cerr << "  number of nodes: " << nb_nodes << std::endl;
        }
      else
        {
          std::cerr << "  bad input (number of nodes): ";
          std::cerr << "  expected an integer, got \'" << line <<
            "\'" << std::endl;
          exit(1);
        }

      // Get node list

      if (filetype == 0)   // ascii
        gmsh_parse_nodes_ascii( is, nb_nodes, Nodes );


      else {                    // binary
        gmsh_parse_nodes_binary( is, nb_nodes, Nodes );
        // Sanity check: is the next character a newline ?
        assert(is.peek() == '\n');
      }

      std::cerr << "done.\n" << std::endl;
    }

    // Read mesh elements
    else if (line.substr(0, 9) == std::string("$Elements")) {
      std::cerr << "Parsing elements..." << std::endl;
      // Get number of elements
      int nb_elements;
      getline( is, line );
      if( parse<int>(nb_elements, line) )
        {
          std::cerr << "  number of elements (including boundaries): ";
          std::cerr << nb_elements << std::endl;
        }
      else
        {
          std::cerr << "  bad input (number of elements): ";
          std::cerr << "  expected an integer, got \'" << line <<
            "\'" << std::endl;
          exit(1);
        }

      // To avoid too many dynamic reallocations of the cells vector
      Cells.reserve( nb_elements );

       int element_counter = 0;
      if ( filetype == 0 ) {   // ascii
        //gmsh_parse_elements_ascii( is, nb_elements, elements );
      }

      else {                    // binary
        gmsh_parse_elements_binary<3,2>( is, nb_elements, BoundaryFaces, Cells );
        assert(is.peek() == '\n');
      }

      std::cerr << "done.\n" << std::endl;
    }

    else                        // ignore other lines
      {}

  }
  NPoints = Nodes.size();
  NVolumes = Cells.size();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int Mesh::readOceanMesh( std::istream &is ) {

  int filetype = 0;             // 0 --> ascii, 1 --> binary
  int datasize = sizeof( double ); // the gmsh standard for now

  if (!is) {
    std::cerr << "Error while reading mesh from stream. Aborting.\n";
    std::abort();
  }

  std::string line = "";                // string holding the current line

  while ( getline(is ,line) ) { // Iterate through stream, line by line
    // Read file header
    if (line.substr(0, 11) == std::string("OceanMesh2D")){

      std::cerr << "Parsing header..." << std::endl;
      // getline(is, line);
      getline( is, line );
			// Get number of nodes and elements
			int nb_nodes, nb_elements;

			std::istringstream isstream(line);
	    std::string token;
	    getline(isstream, token, ' ');
      if( parse<int>(nb_elements, token) )
				{
					std::cerr << "  number of elements: " << nb_elements << std::endl;
				}
			else
				{
					std::cerr << "  bad input (number of elements): ";
					std::cerr << "  expected an integer, got \'" << line <<
						"\'" << std::endl;
					exit(1);
				}
				getline(isstream, token, ' ');
				if( parse<int>(nb_nodes, token) )
					{
						std::cerr << "  number of nodes: " << nb_nodes << std::endl;
					}
				else
					{
						std::cerr << "  bad input (number of nodes): ";
						std::cerr << "  expected an integer, got \'" << line <<
							"\'" << std::endl;
						exit(1);
					}
			std::cerr << "Parsing nodes..." << std::endl;
			// Get node list
			getline(isstream, token, ' ');
			std::cerr << " \'" << token <<
				"\'" << std::endl;

			if (filetype == 0)   {// ascii
				oceanmesh_parse_nodes_ascii( is, nb_nodes, Nodes );
			} else {                    // binary
				// gmsh_parse_nodes_binary( is, nb_nodes, Nodes );
				// // Sanity check: is the next character a newline ?
				// assert(is.peek() == '\n');
				std::cerr << "Bad input for mesh type (expecting an ascii): ";
				exit(1);
			}

			std::cerr << "done parsing nodes.\n" << std::endl;

			std::cerr << "Parsing elements..." << std::endl;
			// To avoid too many dynamic reallocations of the cells vector
			Cells.reserve( nb_elements );
			// Get element list
			if (filetype == 0)  { // ascii
				oceanmesh_parse_elements_ascii( is, nb_elements, Cells );
			} else {                    // binary
				// gmsh_parse_nodes_binary( is, nb_nodes, Nodes );
				// // Sanity check: is the next character a newline ?
				// assert(is.peek() == '\n');
				std::cerr << "Bad input for mesh type (expecting an ascii): ";
				exit(1);
			}

			std::cerr << "done parsing elements.\n" << std::endl;
			//
      // // To avoid too many dynamic reallocations of the cells vector
      // Cells.reserve( nb_elements );
			//
      //  int element_counter = 0;
      // if ( filetype == 0 ) {   // ascii
      //   //gmsh_parse_elements_ascii( is, nb_elements, elements );
      // } else {                    // binary
      //   gmsh_parse_elements_binary<3,2>( is, nb_elements, BoundaryFaces, Cells );
      //   assert(is.peek() == '\n');
      // }
			//
      // std::cerr << "done.\n" << std::endl;
    } else   {                     // ignore other lines
    }

  NPoints = Nodes.size();
  NVolumes = Cells.size();

  return 0;
}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Reverse Cuthill McKee renumbering for bandwith reduction
void Mesh::RCMRenumbering() {

  std::vector<int> parent( Cells.size(), -1);
  unsigned int counter = 0;

  while ( counter < Cells.size() ) {

    // Choice of the initial cell : order 1 vertex, or boundary vertex
    int startCell = -1;

    for ( unsigned int i = 0; i < Cells.size(); ++i ) {

      if ( Cells.at(i).degenerate && parent[i] == -1 ) {
        startCell = i;
        break;
      }

      else if ( Cells.at(i).boundary && parent[i] == -1 )
        startCell = i;

      else if (parent[i] == -1 && startCell == -1 ) {
        startCell = i;
      }

    }

    std::queue<int> next;
    next.push( startCell );
    parent[ startCell ] = startCell;
    counter += 1;

    while ( !next.empty() ) {

      // first element of the queue
      int u = next.front();
      next.pop();

      // loop through neighbours
      boost::array<int,3> neighbors = Cells.at(u).neighbors();

      std::vector< Cell<3,2> > trueNeighbors;

      for (  int j = 0; j < 3; ++j ) {

	if ( neighbors.at(j) != -1 )
	  trueNeighbors.push_back( Cells.at(j) );
	
      }
   

      for ( int j = 0; j < 3; ++j ) {

        int neighbor = neighbors.at( 2 - j );
	
        if ( neighbor != -1 && parent[neighbor] == -1 ) {

          // neighbor is unvisited
          parent[neighbor] = counter;
          next.push( neighbor );
          counter +=1;
        }
      }
    }
  }

  // Reverse the renumbering vector
  //std::reverse( parent.begin(), parent.end() );

  // for (  int i = 0; i < Cells.size(); ++i ) {
  //   std::cerr << parent[i] << " ";
  // }

  std::cerr << std::endl;

  Cells_t NewCells( Cells.size() );

  for ( unsigned int i = 0; i < Cells.size(); ++i ) {

    Cell<3,2> oldCell = Cells.at( i );
    Cell<3,2> cell = Cell<3,2>( oldCell.vertices() );
    cell.boundary = oldCell.boundary;
    cell.degenerate = oldCell.degenerate;
    NewCells.at( parent[i] ) = cell;
  }

  // for (  int i = 0; i < Cells.size(); ++i )
  //   Cells.at(i) = NewCells.at(i);
  Cells.clear();
  Cells = NewCells;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void Mesh::WriteMeshBandwith( std::ofstream &os ) {

  // write mesh connectivity in PPM format.
  if ( Cells.size() > 5000 ) {
    std::cerr << "Mesh size too big, mesh bandwith not written.\n";
  }

  else {
    os << "P2\n";
    os << "# File generated by Volna: mesh bandwith\n";
    os << Cells.size() << "\n";
    os << Cells.size() << "\n";
    os << 255 << "\n";

    for ( unsigned int j = 0; j < Cells.size(); ++j ) {
      for ( unsigned int i = 0; i < Cells.size(); ++i ) {
        boost::array<int,3> neighbors = Cells.at(i).neighbors();
        if ( (unsigned int)neighbors[0] == j || (unsigned int)neighbors[1] == j ||
             (unsigned int)neighbors[2] == j )
          os << 255 << " ";
        else
          os << 0 << " ";
      }
      os << "\n";
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void Mesh::ComputeConnectivity() {
  std::cerr << "Computing mesh connectivity..." << std::endl;

  Facets.clear();
  BoundaryFaces.clear();
  tracked_faces_.clear();

  tracked_faces_.reserve( 3 * Cells.size() );

  // loop through elements and get all the faces (most them are thus duplicated)
  for ( unsigned int i = 0; i < Cells.size(); ++i ) {

    Triangle cell = Cells.at(i);
    boost::array< Face, 3 > cell_facets = cell.compute_facets();

    boost::array< std::pair<int, Face >, 3 > local_tracked_facets;

    for ( unsigned  int j = 0; j < 3; ++j) {
      cell_facets[j].add_cell(i);
      tracked_faces_.push_back (std::pair< int, Face >( i, cell_facets[j] ));
    }
  }

  std::sort( BoundaryFaces.begin(), BoundaryFaces.end() );


  std::sort( tracked_faces_.begin(), tracked_faces_.end(), PairCompare );

  int boundary_facets_counter = 0;
   int face_counter = 0;
  unsigned  int i = 0;

  while (1) {

    if ( i < tracked_faces_.size() - 1 ) {
      int triangle = tracked_faces_.at(i).first;
      Face face = tracked_faces_.at(i).second;
      int partition = face.partition();
      int next_triangle = tracked_faces_.at(i+1).first;
      Face next_face = tracked_faces_.at(i+1).second;
      int next_partition = next_face.partition();

      // This connectivity is to be added for both interior and boundary
      // faces
      Cells.at(triangle).add_facet(face_counter);

      // Interior face
      //
      // In this case, we update the connectivity for
      // both facets <--> triangles and triangles <--> triangles

      if ( face.has_the_same_vertices( next_face ) ){
        Cells.at(next_triangle).add_facet(face_counter);
        // face --> triangle connectivity
        face.add_cell(next_triangle);

        // triangle <--> triangle connectivity
        Cells.at(triangle).add_neighbor(next_triangle);
        Cells.at(next_triangle).add_neighbor(triangle);

        face.add_partition(next_partition);

        i += 2;
      }


      // Boundary face
      //
      // No more connectivity to add. However, we need to get the
      // physical and elementary region tags (for tracking different
      // boundary conditions)
      else {
        // is the current face in the list of boundary facets (read from
        // mesh) ?
        unsigned int j = boundary_facets_counter;

        while ( (j < BoundaryFaces.size()) && (BoundaryFaces.at(j) < face) )
          ++j;

        // if ( (j< boundary_faces_.size()) && (boundary_faces_.at(j).has_the_same_vertices(face))) {
        //   Facet<N-1,d-1> *pboundary_facet = &boundary_faces_[i];
        //   face.set_physical( pboundary_facet->physical() );
        // }
        boundary_facets_counter = j;
        ++i;
      }

      // Finally add the facet to the list of facets
      Facets.push_back( face );
      ++face_counter;

    }


    else if ( i == tracked_faces_.size() - 1) {
      // in this case, the only remaining facet in the list is
      // necessarily a boundary face
      int last_triangle = tracked_faces_.at(i).first;
      // triangle --> facet connectivity
      Cells.at(last_triangle).add_facet( face_counter );
      Face face = tracked_faces_.at(i).second;

      Facets.push_back( face );

      break;
    }

    else if ( i >= tracked_faces_.size() )
      break;
  }

  NFaces = Facets.size();
  std::cerr << "  number of facets: " << Facets.size() << std::endl;
  std::cerr << "done.\n" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::LegacyInterface() {

  for ( unsigned int i = 0; i < Cells.size(); ++i ) {

    // detect boundary Triangles
    if ( Cells.at(i).neighbors()[2] == -1 ) {
      
      Cells.at(i).boundary = true;
      Cells.at(i).boundary_type = -1;

      if ( Cells.at(i).neighbors()[1] == -1 ) {
        Cells.at(i).degenerate = true;
	Cells.at(i).boundary_type = -2;
      }
    }

  }


}


///////////////////////////////////////////////////////////////////////////////

// void Mesh::ComputeGradientInterpolator() {

  // std::cerr << "Computing least square matrices for gradient ";
  // std::cerr << "reconstruction... ";
	//
  // for ( int i = 0; i < NVolumes; ++i ) {
	//
  //   // We construct the least square interpolation matrix as
  //   // M = sum( Di * (Di)^t ), where Di = Xi - X0
	//
  //   LSQMatrix lsq = LSQMatrix::Zero();
	//
  //   for ( unsigned int j = 0; j < ncell; ++j )
	//
  //     if ( Cells.at(i).neighbors()[j] != -1 ) {
  //       Vector DeltaPoints;
	// DeltaPoints.x() =
	//   CellCenters.x( Cells.at(i).neighbors()[j] ) -
  //         CellCenters.x( i );
	// DeltaPoints.y() =
	//   CellCenters.y( Cells.at(i).neighbors()[j] ) -
  //         CellCenters.y( i );
	//
  //       DeltaPoints /= DeltaPoints.norm();
  //       lsq += (DeltaPoints * DeltaPoints.transpose()).block<2,2>(0,0);
  //     }
	//
  //   if ( !Cells.at(i).degenerate ) {
  //     LSQMatrix inverse = LSQMatrix::Zero();
  //     lsq.marked<Eigen::SelfAdjoint>().computeInverse( &inverse );
  //     GradientInterpolator.push_back( inverse );
  //   }
	//
  //   // In this degenerate case, the (boundary) cell has less neighbors than the
  //   // ambient space dimension, and we don't reconstruct the gradient.
  //   // TODO -- reconstruct the gradient using the triangles in the 2-ring
  //   //else
  //     GradientInterpolator.push_back( LSQMatrix::Zero() );
	//
  // }

//   std::cerr << "done.\n";
//
// }

void Mesh::ComputeGeometricQuantities() {

  std::cerr << "Computing mesh geometric quantities... ";

  CellCenters = GeomValues( NVolumes );
  CellVolumes = ScalarValue::Zero( NVolumes );
#pragma omp parallel for
  for (  int i = 0; i < NVolumes; ++i ) {

    // Cell barycenters
    Point center = Point::Zero();

    for (  int j = 0; j < 3; ++j )
      center += Nodes[ Cells.at(i).vertices()[j] ];

    center /= 3.;

    CellCenters.x( i ) = center( 0 );
    CellCenters.y( i ) = center( 1 );
    CellCenters.z( i ) = center( 2 );

    // Cell volumes
    Point p1 = Nodes[ Cells.at(i).vertices()[0] ];
    Point p2 = Nodes[ Cells.at(i).vertices()[1] ];
    Point p3 = Nodes[ Cells.at(i).vertices()[2] ];

    RealType volume = 0.5 * std::abs( orient2d( p1, p2, p3 ) );

    CellVolumes( i ) = volume;

  }

  FacetCenters = GeomValues::GeomValues( NFaces );
  FacetNormals = GeomValues::GeomValues( NFaces );
  FacetVolumes = ScalarValue::Zero( NFaces );
#pragma omp parallel for
  for (  int i = 0; i < NFaces; ++i ) {

    // Facet barycenter
    Point center = Point::Zero();

    for (  int j = 0; j < 2; ++j )
      center += Nodes[ Facets.at(i).vertices()[j] ];

    center /= 2.;

    FacetCenters.x(i) = center( 0 );
    FacetCenters.y(i) = center( 1 );
    FacetCenters.z(i) = center( 2 );

    // Facet volume (e.g. length)
    Point p1 = Nodes[ Facets.at(i).vertices()[0] ];
    Point p2 = Nodes[ Facets.at(i).vertices()[1] ];

    RealType volume = (p2 - p1).norm();

    FacetVolumes( i ) = volume;

    // Facet normal
    Vector normal = Vector::UnitX();
    normal.x() =  - ( p2 - p1 ).y();
    normal.y() = ( p2 - p1).x();
    normal.normalize();

    // Flip the normal if needed.
    // RealType
    //   ScalarProduct = ( FacetCenters.col(i) -
    //                     CellCenters.col( Facets.at(i).LeftCell() ) )
    //   .dot ( normal );

    RealType Orienter =
      ( FacetCenters.x(i) - CellCenters.x( Facets.at(i).LeftCell() ))
      * normal( 0 ) +
      ( FacetCenters.y(i) - CellCenters.y( Facets.at(i).LeftCell() ))
      * normal( 1 );

    if ( Orienter <= 0. )
      normal *= -1;

    FacetNormals.x( i ) = normal( 0 );
    FacetNormals.y( i ) = normal( 1 );
    FacetNormals.z( i ) = 0.;

  }

  std::cerr << "done.\n\n";

}

#endif // MESH_HPP
