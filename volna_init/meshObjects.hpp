#ifndef MESHOBJECTS_HPP
#define MESHOBJECTS_HPP

#include <iostream>
#include <algorithm>
#include <utility>
#include <boost/array.hpp>

#include "precomputed.hpp"

///////////////////////////////////////////////////////////////////////////////

// Generic mesh elements -- parametrized by N: number of element
// vertices, and d: mesh dimension.
//
// N.B. A (minor) problem with this parametrization: some standards (e.g.
// VTK) make a difference between quadrangle and pixels (also hex and
// voxels)
template <size_t N, int d>
class MeshObject {
private:
  static const int gmsh_type_; // class constant
  static const int vtk_type_ ; // class constant
  static const int oceanmesh_type_ ;
  int elementary_, physical_, partition_; // gmsh tags
public:
  int boundary_type;
  boost::array<int,N> vertices_;
  // Constructors
  MeshObject();
  MeshObject( boost::array<int,N> );
  MeshObject( boost::array<int,N>, int );
  // Getters
  boost::array<int,N> vertices();
  int partition();
  int physical() { return physical_; };
  int gmsh_type() { return gmsh_type_; };
  int vtk_type() { return vtk_type_; };
	int oceanmesh_type() { return oceanmesh_type_; };
  // Setters
  void set_partition ( unsigned int partition ) {
    partition_ = partition; };
  void set_physical ( unsigned int physical ) {
    physical_ = physical; };
  // I/O
  template <size_t N2, int d2>
  friend std::ostream& operator << ( std::ostream&,
                                     MeshObject<N,d> );
  // Other
  bool operator < ( const MeshObject<N,d> &) const;
  bool has_the_same_vertices( const MeshObject<N,d>& );
};

///////////////////////////////////////////////////////////////////////////////
// Class variables
template<size_t N, int d>
const int MeshObject<N,d>::gmsh_type_ = GMSH_TYPE<N,d>::value();

template<size_t N, int d>
const int MeshObject<N,d>::vtk_type_ = VTK_TYPE<N,d>::value();

template<size_t N, int d>
const int MeshObject<N,d>::oceanmesh_type_ = OceanMesh_TYPE<N,d>::value();

///////////////////////////////////////////////////////////////////////////////
// Constructors

template <size_t N, int d>
MeshObject<N, d>::MeshObject():
  elementary_(0), physical_(0), partition_(0), boundary_type(0) {}

template <size_t N, int d>
MeshObject<N,d>::MeshObject( boost::array<int,N> vertices ):
  elementary_(0), physical_(0), partition_(0), boundary_type(0),
  vertices_(vertices) {}

template <size_t N, int d>
MeshObject<N,d>::MeshObject( boost::array<int,N> vertices, int partition ):
  elementary_(0), physical_(0), partition_(partition), boundary_type(0),
  vertices_(vertices) {}

///////////////////////////////////////////////////////////////////////////////
// Getters

template <size_t N, int d>
boost::array<int,N> MeshObject<N,d>::vertices() {
  return vertices_;
}

template <size_t N, int d>
int MeshObject<N,d>::partition() {
  return partition_;
}

///////////////////////////////////////////////////////////////////////////////
// I/O

template<size_t N, int d>
std::ostream& operator << ( std::ostream& os,
			    MeshObject<N,d> object ) {

  os << "Vertices:";
  for ( unsigned int i = 0; i < N; ++i ) {
    os << " " << object.vertices()[i]-1;
  }
  os << "\n";
  return os;
}

///////////////////////////////////////////////////////////////////////////////
// Other

template<size_t N, int d>
bool MeshObject<N,d>::operator < ( const MeshObject<N,d> &other_object )
  const {
  boost::array<int,N> vertices1 = vertices_;
  std::sort( vertices1.begin(), vertices1.end() );
  boost::array<int,N> vertices2 = other_object.vertices_;
  std::sort(vertices2.begin(), vertices2.end());
  return ( vertices1 < vertices2 );
}

template<size_t N, int d>
bool MeshObject<N,d>::has_the_same_vertices(const MeshObject<N,d>& other_object) {
  boost::array<int, N> vertices1 = vertices_;
  std::sort(vertices1.begin(), vertices1.end());
  boost::array<int, N> vertices2 = other_object.vertices_;
  std::sort(vertices2.begin(), vertices2.end());
  return (vertices1 == vertices2);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<size_t N, int d>
class Facet : public MeshObject<N,d> {
private:
  boost::array<int,2> partitions_;
  int LeftCell_, RightCell_;
public:
  Facet();
  Facet ( boost::array<int,N> );
  Facet ( boost::array<int,N>, int );
  int LeftCell();
  int const LeftCell() const;
  int RightCell();
  int const RightCell() const;
  boost::array<int,2> partitions();
  // Setters
  void add_cell( int );
  void add_partition( int );
  // I/O
  template <size_t N2, int d2>  // This is not N --- not a bug, but I don't
  // really understand...
  friend std::ostream& operator << ( std::ostream&, const Facet<N,d>& );
  // Other
  bool on_boundary() { return ( RightCell_ < 0 ); };
  const bool on_boundary() const { return ( RightCell_ < 0 ); };
};

///////////////////////////////////////////////////////////////////////////////
// Constructors

template <size_t N, int d>
Facet<N, d>::Facet():
  MeshObject<N,d>(), LeftCell_(-1), RightCell_(-1) {}

template <size_t N, int d>
Facet<N,d>::Facet( boost::array<int,N> vertices ):
  MeshObject<N,d>( vertices ), LeftCell_(-1), RightCell_(-1) {}

template <size_t N, int d>
Facet<N,d>::Facet( boost::array<int,N> vertices, int partition ):
  MeshObject<N,d>( vertices, partition ),
  LeftCell_(-1), RightCell_(-1) { partitions_[0] = partition;}

///////////////////////////////////////////////////////////////////////////////
//
// Getters

template <size_t N, int d>
boost::array<int,2> Facet<N,d>::partitions() {
  return partitions_;
}

template <size_t N, int d>
int Facet<N,d>::LeftCell() {
  return LeftCell_;
}

template <size_t N, int d>
int const Facet<N,d>::LeftCell() const {
  return LeftCell_;
}

template <size_t N, int d>
int Facet<N,d>::RightCell() {
  return RightCell_;
}

template <size_t N, int d>
int const Facet<N,d>::RightCell() const {
  return RightCell_;
}

///////////////////////////////////////////////////////////////////////////////
// Setters

template <size_t N, int d>
void Facet<N, d>::add_cell( int cell_id ) {
  // by convention, LeftCell > RightCell
  if ( LeftCell_ < 0 )
    LeftCell_ = cell_id;
  else {
    RightCell_ = cell_id;
    if ( LeftCell_ < RightCell_ )
      std::swap( LeftCell_, RightCell_ );
  }
}

template <size_t N, int d>
void Facet<N,d>::add_partition( int partition_id ) {
  partitions_[1]=partition_id;
}

///////////////////////////////////////////////////////////////////////////////
// I/O

template <size_t N, int d>
std::ostream& operator << ( std::ostream& os, Facet<N,d> facet ) {

  os << (MeshObject<N,d>) facet;

  os << "Cells: " << facet.LeftCell() << " " << facet.RightCell() <<"\n";

  return os;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<size_t N, int d>
class Cell : public MeshObject<N,d> {
private:
  boost::array<int,N> neighbors_;
  boost::array<int,N> facets_;
public:
  // LEGACY INTERFACE
  bool boundary, degenerate;
  // END LEGACY INTERFACE

  // Constructor
  Cell();
  Cell ( boost::array<int,N> );
  Cell ( boost::array<int,N>, int );
  // Getters
  boost::array< int,N > neighbors() {
    return neighbors_;
  }
  boost::array< int,N > facets() {
    return facets_;
  }
  const boost::array< int,N > facets() const {
    return facets_;
  }
  // Setters
  void add_facet( int );
  void add_neighbor( int );
  void reorder_neighbors( );
  // I/O
  template <size_t N2, int d2>  // This is not N --- not a bug, but I don't
  // really understand...
  friend std::ostream& operator << ( std::ostream&, const Cell<N,d>& );
  // Other
  boost::array< Facet<N-1,d-1>, N >  compute_facets();
};

///////////////////////////////////////////////////////////////////////////////
// Constructors

template <size_t N, int d>
Cell<N,d>::Cell():
  MeshObject<N,d>(), boundary(false), degenerate(false) {
  facets_.assign(-1);
  neighbors_.assign(-1);
};

template <size_t N, int d>
Cell<N, d>::Cell( boost::array<int,N> vertices ):
  MeshObject<N,d>( vertices ), boundary(false), degenerate(false) {
  facets_.assign(-1);
  neighbors_.assign(-1);
};

template <size_t N, int d>
Cell<N, d>::Cell( boost::array<int,N> vertices, int partition ):
  MeshObject<N,d>( vertices, partition ), boundary(false), degenerate(false) {
  facets_.assign(-1);
  neighbors_.assign(-1);
};

///////////////////////////////////////////////////////////////////////////////
// Setters

template <size_t N, int d>
void Cell<N, d>::add_facet( int facet_id ) {
  unsigned int counter = 0;
  while ( facets_[counter] != -1 )
    ++counter;
  facets_[counter] = facet_id;
}

template <size_t N, int d>
void Cell<N, d>::add_neighbor( int neighbor_id ) {
  unsigned int counter = 0;
  while ( neighbors_[counter] != -1 )
    ++counter;
  neighbors_[counter] = neighbor_id;
}

///////////////////////////////////////////////////////////////////////////////
// I/O

template <size_t N, int d>
std::ostream& operator << ( std::ostream& os, Cell<N,d> cell ) {

  os << (MeshObject<N,d>) cell;

  os << "Neighbors:";
  for ( unsigned int i = 0; i < N; ++i ) {
    os << " " << cell.neighbors()[i];
  }

  os << "\nPartition: " << cell.partition() << "\n";

  return os;
}

///////////////////////////////////////////////////////////////////////////////
// Other

template <size_t N, int d>
boost::array< Facet<N-1,d-1>, N > Cell<N,d>::compute_facets() {

  boost::array< Facet<N-1,d-1>, N > faces;

  boost::array< int, N*(N-1) >
    faces_indices = ELEM_FACETS_INDICES<N,d>::value();

  for ( unsigned int i = 0; i < N; ++i ) {

    boost::array<int, N-1> face_indices;

    for ( unsigned int j = (N-1)*i; j < (N-1)*i + N-1; ++j ) {
      face_indices[j- ((N-1)*i)] =  (*this).vertices_[ (faces_indices)[j] ];
    }

    Facet<N-1,d-1> face = Facet<N-1,d-1> ( face_indices, (*this).partition() );
    faces[i] = face;
  }
  return faces;
}



typedef Facet<2,1> Face;
typedef Cell<3,2> Triangle;


#endif // MESHOBJECTS_HPP
