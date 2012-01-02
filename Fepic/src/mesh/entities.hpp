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

#ifndef FEPIC_ENTITIES_HPP
#define FEPIC_ENTITIES_HPP


template<int dim>
class Polytope { public:
  typedef Polytope<dim-1> Derived;

};

template<int dim>
class Simplex { public:
  typedef Simplex<dim-1> Derived;

};

template<int dim>
class Hypercube { public:
  typedef Hypercube<dim-1> Derived;

};

class UndefinedCell
{ public:
  enum { dim=-1,
         n_vertices=0,
         n_nodes=0,
         n_facets=0,
         n_corners=0,
         n_vertices_per_facet=0};

  int _icells_pos[0];
  int _nodes[0];     
  int _facets[0];    
  int _icells[0];    
  int _icells_anchors[0];
  int _corners[0];         
};

#endif





