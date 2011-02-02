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



template<class Traits> class iMesh;

class _Labelable;

// volume
template<class Traits> class _Poly3d;
template<class Traits> class Tetrahedron;
template<class Traits> class iHexahedron;

// surface
template<class Traits> class _Poly2d;
template<class Traits> class Triangle;
template<class Traits> class iQuadrangle;

// line
template<class Traits> class Edge;

//Point
template<class Traits> class Point;

// Halfs
template<class Traits> class HalfEdge;
template<class Traits> class HalfFace;

// Halfls
template<class Traits> class HalfEdgeLab;
template<class Traits> class HalfFaceLab;


template<class ElmType, class Traits> class ElementProperties;

template<int CellDim, class Cell> class VolumeDef; 

template<int Dim, class Cell> class FaceDef;

template<int Dim, class Traits> class HalfDef;

template<int Dim, class Traits> class HalflDef;





/* Empty classes */

template<int dim> class Simplex {};

template<int dim> class Hypercube {};







#endif

