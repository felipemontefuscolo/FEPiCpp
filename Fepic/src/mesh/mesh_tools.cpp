#include "mesh_tools.hpp"
#include "mesh.hpp"
#include <vector>

// da pra melhorar
void MeshTools::removeCell(Cell * cell, Mesh *mesh)
{
  const int cell_dim = mesh->cellDim();
  const int nds_per_cell = mesh->numNodesPerCell();

  const int n_facets_per_cell = mesh->numFacetsPerCell();
  const int n_corners_per_cell = mesh->numCornersPerCell();

  typedef unsigned char MyBool;



  int id;

  // killing nodes
  for (int i = 0; i < nds_per_cell; ++i)
  {
    id = cell->getNodeId(i);
    if (mesh->inSingleCell(mesh->getNode(id)))
      mesh->disablePoint(id);
  }

  // killing facets
  if (cell_dim > 1)
    for (int i = 0; i < n_facets_per_cell; ++i)
    {
      id = cell->getFacetId(i);
      if (mesh->inSingleCell(mesh->getFacet(id)))
        mesh->disableFacet(id);
    }

  // killing corners
  if (cell_dim > 1)
    for (int i = 0; i < n_corners_per_cell; ++i)
    {
      id = cell->getCornerId(i);
      if (mesh->inSingleCell(mesh->getCorner(id)))
        mesh->disableCorner(id);
    }



  Point  * pt;
  Facet  * ft;
  Corner * cr;

  Cell * neib;
  for (int i = 0; i < n_facets_per_cell; ++i)
  {
    const int neib_id = cell->getIncidCell(i);
    if (neib_id < 0)
      continue;
    neib = mesh->getCell( neib_id );

    // nodes
    for (int j = 0; j < nds_per_cell; ++j)
    {
      pt = mesh->getNode(neib->getNodeId(j));
      pt->setIncidCell(neib_id);
      pt->setPosition(j);
    }

    //facets
    if (cell_dim > 1)
      for (int j = 0; j < n_facets_per_cell; ++j)
      {
        ft = mesh->getFacet(neib->getFacetId(j));
        ft->setIncidCell(neib_id);
        ft->setPosition(j);
      }

    //corners
    if (cell_dim > 1)
      for (int j = 0; j < n_corners_per_cell; ++j)
      {
        cr = mesh->getCorner(neib->getCornerId(j));
        cr->setIncidCell(neib_id);
        cr->setPosition(j);
      }
      
    // neib
    int const ith = cell->getIncidCellPos(i);
    neib->setIncidCell(ith, -1);
    
  }


  // killing cells
  id = mesh->getCellId(cell);
  mesh->disableCell(id);
}


