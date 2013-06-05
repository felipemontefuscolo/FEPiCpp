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
#include "Fepic/src/mesh/elements/cell_element.hpp"
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

  void setMesh(Mesh *m) {m_mesh_ptr = m;}
  void setInitialDofId(int fid) {m_initial_dof_id = fid;}
  void setInitialDofAddress(int* a) {m_initial_dof_address = a;}
  void setType(ShapeFunction * sf, int dim=1, int ntags=0, int const*tags=NULL);
  void setType(int ndpv, int ndpr, int ndpf, int ndpc, int ntags=0, int const*tags=NULL);
  void getDivisions(int*& vertices_beg, int*& corners_beg, int*& facets_beg, int*& cells_beg) const;

  void setUp();
  // call this routine when only m_initial_dof_address was changed
  void updateFromInitialDofAddres();

public:

  VarDofs(const char* name, Mesh * m=NULL, int ndpv=0, int ndpr=0, int ndpf=0, int ndpc=0, int fdi=0, int* a=NULL, int ntags=0, int const*tags=NULL)
    : m_name(name), m_mesh_ptr(m), m_vertices_dofs(NULL,0,0), m_corners_dofs(NULL,0,0), m_facets_dofs(NULL,0,0), m_cells_dofs(NULL,0,0)
  {
    m_n_dof_within_vertice = ndpv; // interior
    m_n_dof_within_corner = ndpr;  // interior
    m_n_dof_within_facet = ndpf;   // interior
    m_n_dof_within_cell = ndpc;    // interior
    m_n_dofs = 0;
    m_n_links = 0;
    m_initial_dof_id = fdi;
    m_initial_dof_address = a;

    if (ntags>0)
    {
      FEPIC_CHECK(tags!=NULL, "tags NULL pointer", std::runtime_error);
      m_considered_tags.resize(ntags);
    }
    for (int i = 0; i < ntags; ++i)
      m_considered_tags[i] = tags[i];
    
    // nao precisa disso aqui, pode ser feito no setUp()
    //if (m!=NULL && m->cellDim() < 3)
    //  m_n_dof_within_corner=0;
    
  }

  // users
  int numDofs() const;
  int numDofsPerVertex() const;
  int numDofsPerCell() const;
  int numDofsPerFacet() const;
  int numDofsPerCorner() const;

  void getCellDofs(int *dofs, Cell const*) const;
  void getFacetDofs(int *dofs, CellElement const*) const;
  void getCornerDofs(int *dofs, CellElement const*) const;
  void getVertexDofs(int *dofs, CellElement const*) const;
  
  void getCellAssociatedDofs(int* dofs, Cell const*) const;
  void getFacetAssociatedDofs(int* dofs, CellElement const*) const;
  void getCornerAssociatedDofs(int* dofs, CellElement const*) const;
  void getVertexAssociatedDofs(int* dofs, CellElement const*) const;
  
  char const* getName() const {return m_name.c_str();};

  void linkVertexDofs(Point const* point1, Point const* point2);

  int const* data() const
  {
    if (m_n_dof_within_vertice > 0)
      {return m_vertices_dofs.data(); }
    else if (m_n_dof_within_corner > 0)
      {return m_corners_dofs.data();  }
    else if (m_n_dof_within_facet > 0)
      {return m_facets_dofs.data();  }
    else
      {return m_cells_dofs.data();  }
  }
  int* data()
  {
    if (m_n_dof_within_vertice > 0)
      {return m_vertices_dofs.data(); }
    else if (m_n_dof_within_corner > 0)
      {return m_corners_dofs.data();  }
    else if (m_n_dof_within_facet > 0)
      {return m_facets_dofs.data();  }
    else
      {return m_cells_dofs.data();  }
  }
  int  totalSize() const;

  float getGrowFactor() const {return m_grow_factor;}

protected:
  std::string m_name;
  Mesh*       m_mesh_ptr;
  int         m_n_dof_within_vertice;
  int         m_n_dof_within_corner;
  int         m_n_dof_within_facet;
  int         m_n_dof_within_cell;
  int         m_n_dofs;
  int         m_n_links;
  int         m_initial_dof_id;
  int*        m_initial_dof_address;
  float       m_grow_factor;
  
  std::vector<int> m_considered_tags; /* if m_considered_tags.size()==0, then all tags are considered. */

  Container m_vertices_dofs;  // 0
  Container m_corners_dofs;   // 1
  Container m_facets_dofs;    // 2
  Container m_cells_dofs;     // 3
  

};

#endif
