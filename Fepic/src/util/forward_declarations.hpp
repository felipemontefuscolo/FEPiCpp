// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.

#ifndef FEPIC_FORWARD_DECLARATIONS_HPP
#define FEPIC_FORWARD_DECLARATIONS_HPP



class _Labelable;

// volume
template<class _Traits> class _Poly3d;
template<class _Traits> class Tetrahedron;
template<class _Traits> class Hexahedron;

// surface
template<class _Traits> class _Poly2d;
template<class _Traits> class Triangle;
template<class _Traits> class Quadrangle;

// line
template<class _Traits> class Edge;

//Point
template<class _Traits> class Point;

// Halfs
template<class _Traits> class HalfEdge;
template<class _Traits> class HalfFace;

// Halfs
template<class _Traits> class HalfEdgeLab;
template<class _Traits> class HalfFaceLab;


template<class ElmType, class _Traits> class _MetaCellOf;

template<typename Dim, class _Traits> class _MetaHalfOf;

template<typename Dim, class _Traits> class _MetaHalfLabOf;


template<int dim> class Polytope;

template<int dim> class Simplex;

template<int dim> class Hypercube;



/* mesh */
class _iMeshNameHandler;
template<class _Traits> class _MeshIoMsh;
template<class _Traits> class _MeshIoVtk;
template<class _Traits> class iMesh;

/* Cell Functions */

template<class CellType> int numNodes(int order);
template<class CellType> int numBubbles(int order);


#endif

