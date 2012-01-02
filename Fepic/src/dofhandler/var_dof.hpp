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

#ifndef FEPIC_VAR_DOF_HPP
#define FEPIC_VAR_DOF_HPP

#include "dof_traits.hpp"
#include "../shapefunctions/shape_functions.hpp"
#include "../custom_eigen/custom_eigen.hpp"
#include <vector>
#include <string>

class Mesh;
class ShapeFunction;
class Cell;
class Facet;
class Corner;
class Point;

class VarDofs
{
  friend class DofHandler;

  //typedef std::vector<int> MiniContainer;
  //typedef std::vector<MiniContainer> Container;
  typedef Eigen::Map<Eigen::ArrayXXi> Container;

  void setMesh(Mesh *m) {_mesh_ptr = m;}
  void setInitialDofId(int fid) {_initial_dof_id = fid;}
  void setInitialDofAddress(int* a) {_initial_dof_address = a;}
  void setType(ShapeFunction * sf, int dim=1);
  void setType(int ndpv, int ndpr, int ndpf, int ndpc);
  void getDivisions(int*& vertices_beg, int*& corners_beg, int*& facets_beg, int*& cells_beg) const;

  void setUp();


public:

  VarDofs(const char* name, Mesh * m=NULL, int ndpv=0, int ndpr=0, int ndpf=0, int ndpc=0, int fdi=0, int* a=NULL)
    : _name(name), _mesh_ptr(m), _vertices_dofs(NULL,0,0), _corners_dofs(NULL,0,0), _facets_dofs(NULL,0,0), _cells_dofs(NULL,0,0)
  {
    _n_dof_within_vertice = ndpv; // interior
    _n_dof_within_corner = ndpr;  // interior
    _n_dof_within_facet = ndpf;   // interior
    _n_dof_within_cell = ndpc;    // interior
    _initial_dof_id = fdi;
    _initial_dof_address = a;

  }

  // users
  int numDofs() const;
  int numDofsPerCell() const;
  int numDofsPerFacet() const;
  int numDofsPerCorner() const;

  void getCellDofs(int *dofs, Cell const*) const;
  void getFacetDofs(int *dofs, Facet const*) const;
  void getCornerDofs(int *dofs, Corner const*) const;
  void getVertexDofs(int *dofs, Point const*) const;
  
  void getCellAssociatedDofs(int* dofs, Cell const*) const;
  void getFacetAssociatedDofs(int* dofs, Facet const*) const;
  void getCornerAssociatedDofs(int* dofs, Corner const*) const;
  void getVertexAssociatedDofs(int* dofs, Point const*) const;
  
  char const* getName() const {return _name.c_str();};

  int const* data() const
  {
    if (_n_dof_within_vertice > 0)
      {return _vertices_dofs.data(); }
    else if (_n_dof_within_corner > 0)
      {return _corners_dofs.data();  }
    else if (_n_dof_within_facet > 0)
      {return _facets_dofs.data();  }
    else
      {return _cells_dofs.data();  }
  }
  int* data()
  {
    if (_n_dof_within_vertice > 0)
      {return _vertices_dofs.data(); }
    else if (_n_dof_within_corner > 0)
      {return _corners_dofs.data();  }
    else if (_n_dof_within_facet > 0)
      {return _facets_dofs.data();  }
    else
      {return _cells_dofs.data();  }
  }
  int  totalSize() const;

  float getGrowFactor() const {return _grow_factor;}

protected:
  std::string _name;
  Mesh*       _mesh_ptr;
  int         _n_dof_within_vertice;
  int         _n_dof_within_corner;
  int         _n_dof_within_facet;
  int         _n_dof_within_cell;
  int         _initial_dof_id;
  int*        _initial_dof_address;
  float       _grow_factor;

  Container _vertices_dofs;  // 0
  Container _corners_dofs;   // 1
  Container _facets_dofs;    // 2
  Container _cells_dofs;     // 3


};

#endif
