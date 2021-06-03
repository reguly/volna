#ifndef PRECOMPUTED_HPP
#define PRECOMPUTED_HPP

#include <boost/array.hpp>
// Precompute OceanMesh element type (as function of element dimension, and
// number of vertices) at compile time, using template metaprogramming
// to emulate a switch statement. See :
// http://ubiety.uwaterloo.ca/~tveldhui/papers/Template-Metaprograms/meta-art.html

// Class declarations -- n: number of vertices, d: dimension
template<int n, unsigned int d>
class OceanMesh_TYPE {
public:
  static inline int value() { return 0; }
};

template<>			// point
class OceanMesh_TYPE<1, 0> {
public:
  static inline int value() { return 15; }
};

template<>			// line
class OceanMesh_TYPE<2, 1> {
public:
  static inline int value() { return 1; }
};

template<>			// triangle
class OceanMesh_TYPE<3, 2> {
public:
  static inline int value() { return 2; }
};

template<>			// quadrangle
class OceanMesh_TYPE<4, 2> {
public:
  static inline int value() { return 3; }
};

template<>			// tetrahedron
class OceanMesh_TYPE<4, 3> {
public:
  static inline int value() { return 4; }
};

template<>			// hexahedron
class OceanMesh_TYPE<8, 3> {
public:
  static inline int value() { return 5; }
};

// To get the OceanMesh type for an element: OceanMesh_TYPE<N, d>::value().
// Examples : OceanMesh_TYPE<3, 2>::value() for a triangle, OceanMesh_TYPE<8, 3>
// for an hexahedron.


template<int n, unsigned int d>
class OceanMesh_FACE_TYPE {
public:
  static inline int value() { return n-1; }
};

template<>
class OceanMesh_FACE_TYPE<2, 1> {
public:				// line --> point
  static inline int value() { return 15; }
};

template<>
class OceanMesh_FACE_TYPE<3, 2> {
public:				// triangle --> line
  static inline int value() { return 1; }
};

template<>			// quadrangle --> line
class OceanMesh_FACE_TYPE<4, 2> {
public:
  static inline int value() { return 1; }
};

template<>			// tetrahedron --> triangle
class OceanMesh_FACE_TYPE<4, 3> {
public:
  static inline int value() { return 2; }
};

template<>			// hexahedron --> quadrangle
class OceanMesh_FACE_TYPE<8, 3> {
public:
  static inline int value() { return 4; }
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<int n, unsigned int d>
class OceanMesh_FACE_LENGTH {
public:
  static inline int value() { return n-1; }
};

template<>
class OceanMesh_FACE_LENGTH<2, 1> {
public:				// line --> point
  static inline int value() { return 0; }
};

template<>
class OceanMesh_FACE_LENGTH<3, 2> {
public:				// triangle --> line
  static inline int value() { return 2; }
};

template<>			// quadrangle --> line
class OceanMesh_FACE_LENGTH<4, 2> {
public:
  static inline int value() { return 2; }
};

template<>			// tetrahedron --> triangle
class OceanMesh_FACE_LENGTH<4, 3> {
public:
  static inline int value() { return 3; }
};

template<>			// hexahedron --> quadrangle
class OceanMesh_FACE_LENGTH<8, 3> {
public:
  static inline int value() { return 4; }
};


///////////////////////////////////////////////////////////////////////////////

// Precompute gmsh element type (as function of element dimension, and
// number of vertices) at compile time, using template metaprogramming
// to emulate a switch statement. See :
// http://ubiety.uwaterloo.ca/~tveldhui/papers/Template-Metaprograms/meta-art.html

// Class declarations -- n: number of vertices, d: dimension
template<int n, unsigned int d>
class GMSH_TYPE {
public:
  static inline int value() { return 0; }
};

template<>			// point
class GMSH_TYPE<1, 0> {
public:
  static inline int value() { return 15; }
};

template<>			// line
class GMSH_TYPE<2, 1> {
public:
  static inline int value() { return 1; }
};

template<>			// triangle
class GMSH_TYPE<3, 2> {
public:
  static inline int value() { return 2; }
};

template<>			// quadrangle
class GMSH_TYPE<4, 2> {
public:
  static inline int value() { return 3; }
};

template<>			// tetrahedron
class GMSH_TYPE<4, 3> {
public:
  static inline int value() { return 4; }
};

template<>			// hexahedron
class GMSH_TYPE<8, 3> {
public:
  static inline int value() { return 5; }
};

// To get the gmsh type for an element: GMSH_TYPE<N, d>::value().
// Examples : GMSH_TYPE<3, 2>::value() for a triangle, GMSH_TYPE<8, 3>
// for an hexahedron.


template<int n, unsigned int d>
class GMSH_FACE_TYPE {
public:
  static inline int value() { return n-1; }
};

template<>
class GMSH_FACE_TYPE<2, 1> {
public:				// line --> point
  static inline int value() { return 15; }
};

template<>
class GMSH_FACE_TYPE<3, 2> {
public:				// triangle --> line
  static inline int value() { return 1; }
};

template<>			// quadrangle --> line
class GMSH_FACE_TYPE<4, 2> {
public:
  static inline int value() { return 1; }
};

template<>			// tetrahedron --> triangle
class GMSH_FACE_TYPE<4, 3> {
public:
  static inline int value() { return 2; }
};

template<>			// hexahedron --> quadrangle
class GMSH_FACE_TYPE<8, 3> {
public:
  static inline int value() { return 4; }
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<int n, unsigned int d>
class GMSH_FACE_LENGTH {
public:
  static inline int value() { return n-1; }
};

template<>
class GMSH_FACE_LENGTH<2, 1> {
public:				// line --> point
  static inline int value() { return 0; }
};

template<>
class GMSH_FACE_LENGTH<3, 2> {
public:				// triangle --> line
  static inline int value() { return 2; }
};

template<>			// quadrangle --> line
class GMSH_FACE_LENGTH<4, 2> {
public:
  static inline int value() { return 2; }
};

template<>			// tetrahedron --> triangle
class GMSH_FACE_LENGTH<4, 3> {
public:
  static inline int value() { return 3; }
};

template<>			// hexahedron --> quadrangle
class GMSH_FACE_LENGTH<8, 3> {
public:
  static inline int value() { return 4; }
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// Precompute VTK element type (as function of element dimension, and
// number of vertices) at compile time, using template metaprogramming
// dark magic (emulating a switch statement)
// http://ubiety.uwaterloo.ca/~tveldhui/papers/Template-Metaprograms/meta-art.html

// Class declarations -- n: number of vertices, d: dimension
template<int n, unsigned int d>
class VTK_TYPE {
public:
  static inline int value() { return 0; }
};

template<>			// point
class VTK_TYPE<1, 2> {
public:
  static inline int value() { return 1; }
};

template<>			// line
class VTK_TYPE<2, 2> {
public:
  static inline int value() { return 3; }
};

template<>			// triangle
class VTK_TYPE<3, 2> {
public:
  static inline int value() { return 5; }
};

template<>			// quadrangle
class VTK_TYPE<4, 2> {
public:
  static inline int value() { return 9; }
};

template<>			// tetrahedron
class VTK_TYPE<4, 3> {
public:
  static inline int value() { return 10; }
};

template<>			// hexahedron
class VTK_TYPE<8, 3> {
public:
  static inline int value() { return 12; }
};

template<>			// prism (aka wedge in vtk lingo)
class VTK_TYPE<6, 3> {
public:
  static inline int value() { return 13; }
};

template<>			// pyramid
class VTK_TYPE<5, 3> {
public:
  static inline int value() { return 14; }
};

// To get the vtk type for an element: VTK_TYPE<N, d>::value(),
// where N is the number of vertices, and d the dimension (for
// instance, VTK_TYPE<3, 2>::value() for a triangle, or
// VTK_TYPE<8, 3> for an hexahedron).

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// Precompute METIS element type (as function of element dimension, and
// number of vertices) at compile time, using template metaprogramming
// dark magic (emulating a switch statement)
// http://ubiety.uwaterloo.ca/~tveldhui/papers/Template-Metaprograms/meta-art.html

// Class declarations -- n: number of vertices, d: dimension
template<int n, unsigned int d>
class METIS_TYPE {
public:
  static inline int value() { return 0; }
};

template<>			// triangle
class METIS_TYPE<3, 2> {
public:
  static inline int value() { return 1; }
};

template<>			// quadrangle
class METIS_TYPE<4, 2> {
public:
  static inline int value() { return 4; }
};

template<>			// tetrahedron
class METIS_TYPE<4, 3> {
public:
  static inline int value() { return 2; }
};

template<>			// hexahedron
class METIS_TYPE<8, 3> {
public:
  static inline int value() { return 3; }
};

// To get the vtk type for an element: METIS_TYPE<N, d>::value(),
// where N is the number of vertices, and d the dimension (for
// instance, METIS_TYPE<3, 2>::value() for a triangle, or
// METIS_TYPE<8, 3> for an hexahedron).

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// Precompute elements facets -- return a serialized list of facet indices

// Class declarations -- n: number of vertices, d: dimension
template<size_t N, unsigned int d>
class ELEM_FACETS_INDICES {
public:
  static inline boost::array<int, N> value() { 
    boost::array<int, N> foo;
    return foo; }
};

template<>			// triangle
class ELEM_FACETS_INDICES<3, 2> {
public:
  static inline boost::array<int, 6> value() {
    boost::array<int, 6> result = 
      {{ 0, 1, 1, 2, 2, 0 }};
    return result;
  }
};

template<>			// quadrangle
class ELEM_FACETS_INDICES<4, 2> {
public:
  static inline boost::array<int, 8> value() {
    boost::array<int, 8> result = 
      {{ 0, 1, 1, 2, 2, 3, 3, 0 }};
    return result;
  }
};

template<>			// tetrahedron
class ELEM_FACETS_INDICES<4, 3> {
public:
  static inline boost::array<int, 12> value() {
    boost::array<int, 12> result = 
      {{ 0, 1, 2, 0, 3, 2, 0, 3, 1, 1, 2, 3}};
    return result;
  }
};

template<>			// hexahedron
class ELEM_FACETS_INDICES<8, 3> {
public:
  static inline boost::array<int, 24> value() {
    boost::array<int, 24> result  = 
      { { 0, 3, 2, 1, 0, 3, 5, 4, 0, 4, 7, 3, 
	  1, 2, 6, 5, 2, 3, 7, 6, 4, 5, 6, 7 }};
    return result;
  }
};

#endif // PRECOMPUTED_HPP
