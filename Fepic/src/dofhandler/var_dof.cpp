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

#include "../mesh/mesh.hpp"
#include "var_dof.hpp"
#include "../util/assert.hpp"

void VarDofs::setType(ShapeFunction * sf, int dim)
{
  setType(sf->numDofsAssociatedToVertice()*dim, sf->numDofsAssociatedToCorner()*dim, sf->numDofsAssociatedToFacet()*dim, sf->numDofsAssociatedToCell()*dim);
}

void VarDofs::setType(int ndpv, int ndpr, int ndpf, int ndpc)
{
  _n_dof_within_vertice = ndpv;
  _n_dof_within_corner = ndpr;
  _n_dof_within_facet = ndpf;
  _n_dof_within_cell = ndpc;
}

int VarDofs::totalSize() const
{
  unsigned const n_nodes_total = _mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = _mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = _mesh_ptr->numFacetsTotal();
  unsigned const n_cells_total = _mesh_ptr->numCellsTotal();
  
  int 
  total  = n_nodes_total*_n_dof_within_vertice;
  total += n_corners_total*_n_dof_within_corner;
  total += n_facets_total*_n_dof_within_facet;
  total += n_cells_total*_n_dof_within_cell;
  
  return total;
}

void VarDofs::getDivisions(int*& vertices_beg, int*& corners_beg, int*& facets_beg, int*& cells_beg) const
{
  unsigned const n_nodes_total = _mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = _mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = _mesh_ptr->numFacetsTotal();
  //unsigned const n_cells_total = _mesh_ptr->numCellsTotal();
  
  vertices_beg = _initial_dof_address;
  corners_beg  = vertices_beg  +  n_nodes_total  *_n_dof_within_vertice;
  facets_beg   = corners_beg   +  n_corners_total*_n_dof_within_corner;
  cells_beg    = facets_beg    +  n_facets_total *_n_dof_within_facet;
  
}

void VarDofs::setUp()
{
  unsigned const n_nodes_total = _mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = _mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = _mesh_ptr->numFacetsTotal();
  unsigned const n_cells_total = _mesh_ptr->numCellsTotal();

  unsigned const n_vertices = _mesh_ptr->numVertices();
  unsigned const n_corners = _mesh_ptr->numCorners();
  unsigned const n_facets  = _mesh_ptr->numFacets();
  //unsigned const n_cells  = _mesh_ptr->numCells();

  int* vertices_beg;
  int* corners_beg;
  int* facets_beg;
  int* cells_beg;

  getDivisions(vertices_beg, corners_beg, facets_beg, cells_beg);

  Mesh * mesh = _mesh_ptr;

  // assigns the degrees of freedom to each entity.

  // vertices dof
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      unsigned dof_counter = _initial_dof_id;
      if (_n_dof_within_vertice > 0)
      {
        new (&_vertices_dofs) Container(vertices_beg, n_nodes_total, _n_dof_within_vertice);

        Point const*p;
        for (unsigned i = 0; i < n_nodes_total; ++i)
        {
          p = mesh->getNode(i);
          if (p->disabled() || !(mesh->isVertex(p)))
            continue;
          for (unsigned j = 0; j < _vertices_dofs.cols(); ++j)
            _vertices_dofs(i,j) = dof_counter++;
        }
      }
    }

    // corners dof
    #pragma omp section
    {
      unsigned dof_counter = n_vertices*_n_dof_within_vertice + _initial_dof_id;
      if (_n_dof_within_corner > 0)
      {
        new (&_corners_dofs) Container(corners_beg, n_corners_total, _n_dof_within_corner);

        Corner const*p;
        for (unsigned i = 0; i < n_corners_total; ++i)
        {
          p = mesh->getCorner(i);
          if (p->disabled())
            continue;
          for (unsigned j = 0; j < _corners_dofs.cols(); ++j)
            _corners_dofs(i,j) = dof_counter++;
        }
      }
    }

    // facets dof
    #pragma omp section
    {
      unsigned dof_counter = n_vertices*_n_dof_within_vertice + n_corners*_n_dof_within_corner + _initial_dof_id;
      if (_n_dof_within_facet > 0)
      {
        new (&_facets_dofs) Container(facets_beg, n_facets_total, _n_dof_within_facet);

        Facet const*p;
        for (unsigned i = 0; i < n_facets_total; ++i)
        {
          p = mesh->getFacet(i);
          if (p->disabled())
            continue;
          for (unsigned j = 0; j < _facets_dofs.cols(); ++j)
            _facets_dofs(i,j) = dof_counter++;
        }
      }

    }

    // cells dof
    #pragma omp section
    {
      unsigned dof_counter = n_vertices*_n_dof_within_vertice +
                             n_corners*_n_dof_within_corner +
                             n_facets*_n_dof_within_facet +
                             _initial_dof_id;
      if (_n_dof_within_cell > 0)
      {
        new (&_cells_dofs) Container(cells_beg, n_cells_total, _n_dof_within_cell);

        Cell const*p;
        for (unsigned i = 0; i < n_cells_total; ++i)
        {
          p = mesh->getCell(i);
          if (p->disabled())
            continue;
          for (unsigned j = 0; j < _cells_dofs.cols(); ++j)
            _cells_dofs(i,j) = dof_counter++;
        }
      }

    }

  }

}


int VarDofs::numDofs() const
{
  unsigned const n_vertices = _mesh_ptr->numVertices();
  unsigned const n_corners = _mesh_ptr->numCorners();
  unsigned const n_facets  = _mesh_ptr->numFacets();
  unsigned const n_cells  = _mesh_ptr->numCells();

  return n_vertices*_n_dof_within_vertice +
         n_corners*_n_dof_within_corner +
         n_facets*_n_dof_within_facet +
         n_cells*_n_dof_within_cell;
}



int VarDofs::numDofsPerCell() const
{
  return _n_dof_within_vertice*_mesh_ptr->numVerticesPerCell() +
         _n_dof_within_corner *_mesh_ptr->numCornersPerCell() +
         _n_dof_within_facet  *_mesh_ptr->numFacetsPerCell() +
         _n_dof_within_cell;

}


int VarDofs::numDofsPerFacet() const
{
  // 3D
  int const num_corners_per_facet = _mesh_ptr->numVerticesPerFacet();

  return _n_dof_within_vertice*_mesh_ptr->numVerticesPerFacet() +
         _n_dof_within_corner *         num_corners_per_facet +
         _n_dof_within_facet;

}

int VarDofs::numDofsPerCorner() const
{
  return _n_dof_within_vertice*_mesh_ptr->numVerticesPerCorner() +
         _n_dof_within_corner;

}

void VarDofs::getCellDofs(int *dofs, Cell const* cell) const
{
  int const n_vtcs_p_cell = _mesh_ptr->numVerticesPerCell();
  int const n_crns_p_cell = _mesh_ptr->numCornersPerCell();
  int const n_fcts_p_cell = _mesh_ptr->numFacetsPerCell();

  // vertices
  for (int i = 0; i < n_vtcs_p_cell; ++i)
  {
    const int vtx_id = cell->getNodeId(i);

    for (int j = 0; j < _n_dof_within_vertice; ++j)
    {
      *dofs++ = _vertices_dofs(vtx_id,j);
    }
  }

  // corners
  for (int i = 0; i < n_crns_p_cell; ++i)
  {
    const int crn_id = cell->getCornerId(i);

    for (int j = 0; j < _n_dof_within_corner; ++j)
    {
      *dofs++ = _corners_dofs(crn_id,j);
    }
  }

  // facets
  for (int i = 0; i < n_fcts_p_cell; ++i)
  {
    const int fct_id = cell->getFacetId(i);

    for (int j = 0; j < _n_dof_within_facet; ++j)
    {
      *dofs++ = _facets_dofs(fct_id,j);
    }
  }

  // cells
  {
    const int cell_id = _mesh_ptr->getCellId(cell);

    for (int j = 0; j < _n_dof_within_cell; ++j)
    {
      *dofs++ = _cells_dofs(cell_id,j);
    }
  }

}

void VarDofs::getFacetDofs(int *dofs, Facet const* facet) const
{
  
  int const n_vtcs_p_facet = _mesh_ptr->numVerticesPerFacet();
  int const n_crns_p_facet = n_vtcs_p_facet; // belive

  Cell const* icell = _mesh_ptr->getCell(facet->getIncidCell());
  const int pos = facet->getPosition();

  int vtcs_ids[n_vtcs_p_facet];
  icell->getFacetVerticesId(pos, vtcs_ids);

  // vertices
  for (int i = 0; i < n_vtcs_p_facet; ++i)
  {
    for (int j = 0; j < _n_dof_within_vertice; ++j)
    {
      *dofs++ = _vertices_dofs(vtcs_ids[i],j);
    }
  }

  int crns_ids[n_crns_p_facet];
  icell->getFacetCornersId(pos, crns_ids);
  
  // corners
  for (int i = 0; i < n_crns_p_facet; ++i)
  {
    for (int j = 0; j < _n_dof_within_corner; ++j)
    {
      *dofs++ = _corners_dofs(crns_ids[i],j);
    }
  }

  // facets
  {
    const int fct_id = _mesh_ptr->getFacetId(facet);

    for (int j = 0; j < _n_dof_within_facet; ++j)
    {
      *dofs++ = _facets_dofs(fct_id,j);
    }
  }

}

void VarDofs::getCornerDofs(int *dofs, Corner const* corner) const
{
  int const n_vtcs_p_corner = _mesh_ptr->numVerticesPerCorner();

  Cell const* icell = _mesh_ptr->getCell(corner->getIncidCell());
  const int pos = corner->getPosition();

  int vtcs_ids[n_vtcs_p_corner];
  icell->getCornerVerticesId(pos, vtcs_ids);


  // vertices
  for (int i = 0; i < n_vtcs_p_corner; ++i)
  {
    for (int j = 0; j < _n_dof_within_vertice; ++j)
    {
      *dofs++ = _vertices_dofs(vtcs_ids[i],j);
    }
  }

  // corners
  {
    const int crn_id = _mesh_ptr->getCornerId(corner);

    for (int j = 0; j < _n_dof_within_corner; ++j)
    {
      *dofs++ = _corners_dofs(crn_id,j);
    }
  }

}

void VarDofs::getVertexDofs(int *dofs, Point const* point) const
{
  // vertices
  {
    const int pt_id = _mesh_ptr->getPointId(point);
    
    for (int j = 0; j < _n_dof_within_vertice; ++j)
    {
      *dofs++ = _vertices_dofs(pt_id,j);
    }
  }

}



void VarDofs::getCellAssociatedDofs(int* dofs, Cell const* cell) const
{
  const int cell_id = _mesh_ptr->getCellId(cell);
  
  for (int j = 0; j < _n_dof_within_cell; ++j)
    *dofs++ = _cells_dofs(cell_id,j);
}
void VarDofs::getFacetAssociatedDofs(int* dofs, Facet const* facet) const
{
  const int fct_id = _mesh_ptr->getFacetId(facet);
  
  for (int j = 0; j < _n_dof_within_facet; ++j)
    *dofs++ = _facets_dofs(fct_id,j);
}
void VarDofs::getCornerAssociatedDofs(int* dofs, Corner const* corner) const
{
  const int crn_id = _mesh_ptr->getCornerId(corner);
  
  for (int j = 0; j < _n_dof_within_corner; ++j)
    *dofs++ = _corners_dofs(crn_id,j);
}
void VarDofs::getVertexAssociatedDofs(int* dofs, Point const* point) const
{
  getVertexDofs(dofs, point);
}



























