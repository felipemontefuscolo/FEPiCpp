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

void VarDofs::setType(ShapeFunction * sf, int dim, int ntags, int const*tags)
{
  setType(sf->numDofsAssociatedToVertice()*dim, sf->numDofsAssociatedToCorner()*dim, sf->numDofsAssociatedToFacet()*dim, sf->numDofsAssociatedToCell()*dim,ntags,tags);
}

void VarDofs::setType(int ndpv, int ndpr, int ndpf, int ndpc, int ntags, int const*tags)
{
  m_n_dof_within_vertice = ndpv;
  m_n_dof_within_corner = ndpr;
  m_n_dof_within_facet = ndpf;
  m_n_dof_within_cell = ndpc;

  if (ntags>0)
  {
    FEPIC_CHECK(tags!=NULL, "tags NULL pointer", std::runtime_error);
    m_considered_tags.resize(ntags);
  }
  for (int i = 0; i < ntags; ++i)
    m_considered_tags[i] = tags[i];  
}

int VarDofs::totalSize() const
{
  unsigned const n_nodes_total = m_mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = m_mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = m_mesh_ptr->numFacetsTotal();
  unsigned const n_cells_total = m_mesh_ptr->numCellsTotal();
  
  int 
  total  = n_nodes_total*m_n_dof_within_vertice;
  total += n_corners_total*m_n_dof_within_corner;
  total += n_facets_total*m_n_dof_within_facet;
  total += n_cells_total*m_n_dof_within_cell;
  
  return total;
}

void VarDofs::getDivisions(int*& vertices_beg, int*& corners_beg, int*& facets_beg, int*& cells_beg) const
{
  unsigned const n_nodes_total = m_mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = m_mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = m_mesh_ptr->numFacetsTotal();
  //unsigned const n_cells_total = m_mesh_ptr->numCellsTotal();
  
  vertices_beg = m_initial_dof_address;
  corners_beg  = vertices_beg  +  n_nodes_total  *m_n_dof_within_vertice;
  facets_beg   = corners_beg   +  n_corners_total*m_n_dof_within_corner;
  cells_beg    = facets_beg    +  n_facets_total *m_n_dof_within_facet;
  
}

void VarDofs::setUp()
{
  unsigned const n_nodes_total = m_mesh_ptr->numNodesTotal();
  unsigned const n_corners_total = m_mesh_ptr->numCornersTotal();
  unsigned const n_facets_total = m_mesh_ptr->numFacetsTotal();
  unsigned const n_cells_total = m_mesh_ptr->numCellsTotal();

  //unsigned const n_vertices = m_mesh_ptr->numVertices();
  //unsigned const n_corners = m_mesh_ptr->numCorners();
  //unsigned const n_facets  = m_mesh_ptr->numFacets();
  //unsigned const n_cells  = m_mesh_ptr->numCells();

  int* vertices_beg;
  int* corners_beg;
  int* facets_beg;
  int* cells_beg;
  
  int tag;
  bool is_considered;

  getDivisions(vertices_beg, corners_beg, facets_beg, cells_beg);

  Mesh * mesh = m_mesh_ptr;

  if (m_mesh_ptr->cellDim() < 3)
    m_n_dof_within_corner = 0;

  unsigned dof_counter = m_initial_dof_id;

  // vertices dof
  if (m_n_dof_within_vertice > 0)
  {
    new (&m_vertices_dofs) Container(vertices_beg, n_nodes_total, m_n_dof_within_vertice);

    Point const*p;
    for (unsigned i = 0; i < n_nodes_total; ++i)
    {
      p = mesh->getNodePtr(i);
      tag = p->getTag();
      
      // check for tag
      is_considered = m_considered_tags.empty() ? true : checkValue(m_considered_tags.begin(), m_considered_tags.end(), tag);
      
      if (!(mesh->isVertex(p)) || !is_considered || p->isDisabled() )
        continue;
      for (unsigned j = 0; j < m_vertices_dofs.cols(); ++j)
        m_vertices_dofs(i,j) = dof_counter++;
    }
  }

  // corners dof
  if (m_n_dof_within_corner > 0)
  {
    new (&m_corners_dofs) Container(corners_beg, n_corners_total, m_n_dof_within_corner);

    Corner const*p;
    for (unsigned i = 0; i < n_corners_total; ++i)
    {
      p = mesh->getCornerPtr(i);
      tag = p->getTag();
      
      // check for tag
      is_considered = m_considered_tags.empty() ? true : checkValue(m_considered_tags.begin(), m_considered_tags.end(), tag);
      
      if (!is_considered || p->isDisabled())
        continue;
      for (unsigned j = 0; j < m_corners_dofs.cols(); ++j)
        m_corners_dofs(i,j) = dof_counter++;
    }
  }

  // facets dof
  if (m_n_dof_within_facet > 0)
  {
    new (&m_facets_dofs) Container(facets_beg, n_facets_total, m_n_dof_within_facet);

    Facet const*p;
    for (unsigned i = 0; i < n_facets_total; ++i)
    {
      p = mesh->getFacetPtr(i);
      tag = p->getTag();
      
      // check for tag
      is_considered = m_considered_tags.empty() ? true : checkValue(m_considered_tags.begin(), m_considered_tags.end(), tag);
      
      if (!is_considered || p->isDisabled())
        continue;
      for (unsigned j = 0; j < m_facets_dofs.cols(); ++j)
        m_facets_dofs(i,j) = dof_counter++;
    }
  }


  // cells dof
  if (m_n_dof_within_cell > 0)
  {
    new (&m_cells_dofs) Container(cells_beg, n_cells_total, m_n_dof_within_cell);

    Cell const*p;
    for (unsigned i = 0; i < n_cells_total; ++i)
    {
      p = mesh->getCellPtr(i);
      tag = p->getTag();
      
      // check for tag
      is_considered = m_considered_tags.empty() ? true : checkValue(m_considered_tags.begin(), m_considered_tags.end(), tag);
      
      if (!is_considered || p->isDisabled())
        continue;
      for (unsigned j = 0; j < m_cells_dofs.cols(); ++j)
        m_cells_dofs(i,j) = dof_counter++;
    }
  }

  m_n_dofs = dof_counter - m_initial_dof_id;

}


int VarDofs::numDofs() const
{
  return m_n_dofs;
}



int VarDofs::numDofsPerCell() const
{
  return m_n_dof_within_vertice*m_mesh_ptr->numVerticesPerCell() +
         m_n_dof_within_corner *m_mesh_ptr->numCornersPerCell() +
         m_n_dof_within_facet  *m_mesh_ptr->numFacetsPerCell() +
         m_n_dof_within_cell;

}


int VarDofs::numDofsPerFacet() const
{
  // 3D
  int const numm_corners_per_facet = m_mesh_ptr->numVerticesPerFacet();

  return m_n_dof_within_vertice*m_mesh_ptr->numVerticesPerFacet() +
         m_n_dof_within_corner *         numm_corners_per_facet +
         m_n_dof_within_facet;

}

int VarDofs::numDofsPerCorner() const
{
  return m_n_dof_within_vertice*m_mesh_ptr->numVerticesPerCorner() +
         m_n_dof_within_corner;

}

void VarDofs::getCellDofs(int *dofs, Cell const* cell) const
{
  int const n_vtcs_p_cell = m_mesh_ptr->numVerticesPerCell();
  int const n_crns_p_cell = m_mesh_ptr->numCornersPerCell();
  int const n_fcts_p_cell = m_mesh_ptr->numFacetsPerCell();

  // vertices
  for (int i = 0; i < n_vtcs_p_cell; ++i)
  {
    const int vtx_id = cell->getNodeId(i);

    for (int j = 0; j < m_n_dof_within_vertice; ++j)
    {
      *dofs++ = m_vertices_dofs(vtx_id,j);
    }
  }

  // corners
  for (int i = 0; i < n_crns_p_cell; ++i)
  {
    const int crn_id = cell->getCornerId(i);

    for (int j = 0; j < m_n_dof_within_corner; ++j)
    {
      *dofs++ = m_corners_dofs(crn_id,j);
    }
  }

  // facets
  for (int i = 0; i < n_fcts_p_cell; ++i)
  {
    const int fct_id = cell->getFacetId(i);

    for (int j = 0; j < m_n_dof_within_facet; ++j)
    {
      *dofs++ = m_facets_dofs(fct_id,j);
    }
  }

  // cells
  {
    const int cell_id = m_mesh_ptr->getCellId(cell);

    for (int j = 0; j < m_n_dof_within_cell; ++j)
    {
      *dofs++ = m_cells_dofs(cell_id,j);
    }
  }

}

void VarDofs::getFacetDofs(int *dofs, CellElement const* facet) const
{
  
  int const n_vtcs_p_facet = m_mesh_ptr->numVerticesPerFacet();
  int const n_crns_p_facet = n_vtcs_p_facet; // belive

  Cell const* icell = m_mesh_ptr->getCellPtr(facet->getIncidCell());
  const int pos = facet->getPosition();

  //int vtcs_ids[n_vtcs_p_facet];
  int *vtcs_ids = new int [n_vtcs_p_facet];
  icell->getFacetVerticesId(pos, vtcs_ids);

  // vertices
  for (int i = 0; i < n_vtcs_p_facet; ++i)
  {
    for (int j = 0; j < m_n_dof_within_vertice; ++j)
    {
      *dofs++ = m_vertices_dofs(vtcs_ids[i],j);
    }
  }

  //int crns_ids[n_crns_p_facet];
  int *crns_ids = new int [n_crns_p_facet];
  icell->getFacetCornersId(pos, crns_ids);
  
  // corners
  for (int i = 0; i < n_crns_p_facet; ++i)
  {
    for (int j = 0; j < m_n_dof_within_corner; ++j)
    {
      *dofs++ = m_corners_dofs(crns_ids[i],j);
    }
  }

  // facets
  {
    //const int fct_id = m_mesh_ptr->getFacetId(facet);
    const int fct_id = m_mesh_ptr->getCellPtr(facet->getIncidCell())->getFacetId(facet->getPosition());

    for (int j = 0; j < m_n_dof_within_facet; ++j)
    {
      *dofs++ = m_facets_dofs(fct_id,j);
    }
  }

  delete [] vtcs_ids;
  vtcs_ids = NULL;
  delete [] crns_ids;
  crns_ids = NULL;
}

void VarDofs::getCornerDofs(int *dofs, CellElement const* corner) const
{
  int const n_vtcs_p_corner = m_mesh_ptr->numVerticesPerCorner();

  Cell const* icell = m_mesh_ptr->getCellPtr(corner->getIncidCell());
  const int pos = corner->getPosition();

  //int vtcs_ids[n_vtcs_p_corner];
  int *vtcs_ids = new int [n_vtcs_p_corner];
  icell->getCornerVerticesId(pos, vtcs_ids);


  // vertices
  for (int i = 0; i < n_vtcs_p_corner; ++i)
  {
    for (int j = 0; j < m_n_dof_within_vertice; ++j)
    {
      *dofs++ = m_vertices_dofs(vtcs_ids[i],j);
    }
  }

  // corners
  {
    //const int crn_id = m_mesh_ptr->getCornerId(corner);
    const int crn_id = m_mesh_ptr->getCellPtr(corner->getIncidCell())->getCornerId(corner->getPosition());

    for (int j = 0; j < m_n_dof_within_corner; ++j)
    {
      *dofs++ = m_corners_dofs(crn_id,j);
    }
  }
  
  delete [] vtcs_ids;
  vtcs_ids = NULL;
}

void VarDofs::getVertexDofs(int *dofs, CellElement const* point) const
{
  // vertices
  {
    //const int pt_id = m_mesh_ptr->getPointId(point);
    const int pt_id = m_mesh_ptr->getCellPtr(point->getIncidCell())->getNodeId(point->getPosition());
    
    for (int j = 0; j < m_n_dof_within_vertice; ++j)
    {
      *dofs++ = m_vertices_dofs(pt_id,j);
    }
  }

}



void VarDofs::getCellAssociatedDofs(int* dofs, Cell const* cell) const
{
  const int cell_id = m_mesh_ptr->getCellId(cell);
  
  for (int j = 0; j < m_n_dof_within_cell; ++j)
    *dofs++ = m_cells_dofs(cell_id,j);
}
void VarDofs::getFacetAssociatedDofs(int* dofs, CellElement const* facet) const
{
  //const int fct_id = m_mesh_ptr->getFacetId(facet);
  const int fct_id = m_mesh_ptr->getCellPtr(facet->getIncidCell())->getFacetId(facet->getPosition());
  
  for (int j = 0; j < m_n_dof_within_facet; ++j)
    *dofs++ = m_facets_dofs(fct_id,j);
}
void VarDofs::getCornerAssociatedDofs(int* dofs, CellElement const* corner) const
{
  //const int crn_id = m_mesh_ptr->getCornerId(corner);
  const int crn_id = m_mesh_ptr->getCellPtr(corner->getIncidCell())->getCornerId(corner->getPosition());
  
  for (int j = 0; j < m_n_dof_within_corner; ++j)
    *dofs++ = m_corners_dofs(crn_id,j);
}
void VarDofs::getVertexAssociatedDofs(int* dofs, CellElement const* point) const
{
  getVertexDofs(dofs, point);
}




void VarDofs::linkVertexDofs(Point const* point1, Point const* point2)
{
  const int pt_id1 = m_mesh_ptr->getCellPtr(point1->getIncidCell())->getNodeId(point1->getPosition());
  const int pt_id2 = m_mesh_ptr->getCellPtr(point2->getIncidCell())->getNodeId(point2->getPosition());
  
  for (int j = 0; j < m_n_dof_within_vertice; ++j)
  {
    FEPIC_CHECK(m_vertices_dofs(pt_id2,j) >= 0, "Point has no associated dofs", std::invalid_argument);
    m_vertices_dofs(pt_id1,j) = m_vertices_dofs(pt_id2,j);
  }
  
}






















