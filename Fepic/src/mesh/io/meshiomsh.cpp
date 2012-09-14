#include "meshiomsh.hpp"
#include <cstdio>
#include "../../util/assert.hpp"



ECellType MeshIoMsh::identifiesMshMeshType(const char* filename, int &space_dim) const
{
  /*
   * O que é feito:
   *
   * 1) primeiramente lê-se apenas as células (conectividade clássica);
   * 2) depois são construídos as facet-edges e(ou) facet-faces, lendo-se
   *    os elementos de dimensões menores.
   *
   * */

  FILE * file_ptr = fopen(filename, "r");

  FEPIC_ASSERT(file_ptr, "can not find mesh file", std::invalid_argument);

  double  coord[3];
  int     type_tag;
  char    buffer[256];

  long    nodes_file_pos;  // $Nodes
  long    elems_file_pos;  // $Elements

  // nós
  nodes_file_pos = find_keyword("$Nodes", 6, file_ptr);
  //nodes_file_pos++;  // only to avoid gcc warning: variable ‘nodes_file_pos’ set but not used [-Wunused-but-set-variable]

  FEPIC_ASSERT(nodes_file_pos>0, "invalid file format", std::invalid_argument);

  space_dim = 1;

  int num_pts(0);
  int node_number(0);
  fscanf(file_ptr, "%d", &num_pts);

  fgets(buffer, sizeof(buffer), file_ptr); // escapa do \n
  for (int i=0; i< num_pts; ++i)
  {
    fgets(buffer, sizeof(buffer), file_ptr);
    sscanf(buffer, "%d %lf %lf %lf", &node_number, &coord[0], &coord[1], &coord[2]);
    FEPIC_ASSERT(node_number==i+1, "wrong file format", std::invalid_argument);
    if (coord[1] != 0)
      space_dim = 2;
    if (coord[2] != 0)
    {
      space_dim = 3;
      break;
    }
  }

  // contagem de elementos e alocação
  elems_file_pos = find_keyword("$Elements", 9, file_ptr);

  FEPIC_ASSERT(elems_file_pos>0, "invalid file format", std::invalid_argument);

//  int num_cells=0;
  int num_elms;

  fscanf(file_ptr, "%d", &num_elms);

  /* ---------------------------------------
   * Detectando a ordem da malha, verificando sequencia dos elementos,
   * e contando o número de células.
   * --------------------------------------- */
//  bool wrong_file_err=true;
  int  elem_number, elm_dim, num_tags, physical;
  int current_elm_dim = 0;
  ECellType current_elm_type = UNDEFINED_CELLT;
  fgets(buffer, sizeof(buffer), file_ptr); // escapa do \n
  for (int k=0; k < num_elms; ++k)
  {

    fscanf(file_ptr, "%d %d %d %d", &elem_number, &type_tag, &num_tags, &physical);

    //// sincronização
    FEPIC_ASSERT(elem_number==k+1, "invalid file format", std::invalid_argument);

    for (int j=1; j<num_tags; ++j)
    {
      fscanf(file_ptr, "%s", buffer);
    }

    elm_dim = dimForMshTag(EMshTag(type_tag));

    if(elm_dim > current_elm_dim)
    {
      current_elm_type = mshTag2ctype(EMshTag(type_tag));
      current_elm_dim = elm_dim;
      if (elm_dim==3)
        break;
    }
  
    fgets(buffer, sizeof(buffer), file_ptr);

  }

  fclose(file_ptr);
  
  return current_elm_type;

}




/* A malha deve estar alocada.
 * */
void MeshIoMsh::readFileMsh(const char* filename, Mesh * mesh)
{
  /*
   * O que é feito:
   *
   * 1) primeiramente lê-se apenas as células (conectividade clássica);
   * 2) depois são construídos as facet-edges e(ou) facet-faces, lendo-se
   *    os elementos de dimensões menores.
   *
   * */

  FEPIC_ASSERT(mesh, "invalid mesh pointer", std::invalid_argument);


  this->_registerFile(filename, ".msh");

  FILE * file_ptr = fopen(filename, "r");

  double  coord[3];
  int     type_tag;
  char    buffer[256];

  int const spacedim = mesh->spaceDim();
  FEPIC_ASSERT(spacedim>0 && spacedim < 4, "mesh has invalid spacedim", std::invalid_argument);

  long    nodes_file_pos;  // $Nodes
  long    elems_file_pos;  // $Elements

  // nós
  nodes_file_pos = find_keyword("$Nodes", 6, file_ptr);
  //nodes_file_pos++;  // only to avoid gcc warning: variable ‘nodes_file_pos’ set but not used [-Wunused-but-set-variable]

  FEPIC_ASSERT(nodes_file_pos>0, "invalid file format", std::invalid_argument);

  int num_pts(0);
  int node_number(0);
  fscanf(file_ptr, "%d", &num_pts);

  //mesh->_pointL_original_size = num_pts;
  mesh->resizePointL(num_pts);
  
  //mesh->printInfo();
  //std::cout << "DEBUGGGGGGGGGGGGGGGGGGG:  "<<mesh << std::endl;
  //printf("DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGG num_pts=%d; numNodesTotal()=%d; numNodes()=%d\n",num_pts,mesh->numNodesTotal(), mesh->numNodes());
  
  fgets(buffer, sizeof(buffer), file_ptr); // escapa do \n
  for (int i=0; i< num_pts; ++i)
  {
    fgets(buffer, sizeof(buffer), file_ptr);
    sscanf(buffer, "%d %lf %lf %lf", &node_number, &coord[0], &coord[1], &coord[2]);
    FEPIC_ASSERT(node_number==i+1, "wrong file format", std::invalid_argument);
    mesh->getNode(i)->setCoord(coord);
  }
  // os pontos não estão completas: falta atribuir os labels

  // contagem de elementos e alocação
  
  elems_file_pos = find_keyword("$Elements", 9, file_ptr);

  FEPIC_ASSERT(elems_file_pos>0, "invalid file format", std::invalid_argument);

  int num_cells=0;
  int num_elms;

  fscanf(file_ptr, "%d", &num_elms);

  /* ---------------------------------------
   * Detectando a ordem da malha, verificando sequencia dos elementos,
   * e contando o número de células.
   * --------------------------------------- */
  const int mesh_cell_msh_tag = mesh->cellMshTag();
  bool wrong_file_err=true;
  int  elem_number;               // contagem para verificação de erros.
  for (int k = 0; k < num_elms; ++k)
  {
    fscanf(file_ptr, "%d %d", &elem_number, &type_tag);
  
    // check sequence
    if (elem_number != k+1)
    {
      wrong_file_err=true;
      break;
    }

    // check type
    if (type_tag == mesh_cell_msh_tag)
    {
      wrong_file_err=false;
      ++num_cells;
    }

    if (!fgets(buffer, sizeof(buffer), file_ptr) || feof(file_ptr))
    {
      wrong_file_err=true;
      break;
    }
    
  }
  FEPIC_ASSERT(!wrong_file_err, "Wrong file format. Make sure you created the mesh with correct file format. ", std::invalid_argument);

  Cell*  cell;
  //cell = Cell::create(mesh->cellType());

  {
    unsigned aux;
    //mesh->_cellL_original_size = num_cells;
    mesh->resizeCellL(num_cells);
    aux = mesh->estimateNumFacets(num_cells, mesh->cellType());
    //mesh->_facetL_original_size = aux;
    mesh->reserveFacetL(aux);
    if (mesh->cellDim() > 2)
    {
      aux = mesh->estimateNumCorners(num_cells, mesh->cellType());
      //mesh->_cornerL_original_size = aux;
      mesh->reserveCornerL(aux);
    }
  }

  /* --------------------------------------
   * Lendo as células
   * -------------------------------------- */
  fseek (file_ptr , elems_file_pos , SEEK_SET );
  fscanf(file_ptr, "%d", &num_elms);

  this->timer.restart();

  int  inc(0);
  int  nodes_per_cell = mesh->nodesPerCell();
  int  id_aux;
  int  num_tags;
  int  elm_dim;
  int  physical;
  int const cell_dim = mesh->cellDim();
  int const cell_msh_tag = mesh->cellMshTag();
  for (int k=0; k < num_elms; ++k)
  {
    fscanf(file_ptr, "%d %d %d %d", &elem_number, &type_tag, &num_tags, &physical);

    // sincronização
    FEPIC_ASSERT(elem_number==k+1, "invalid file format", std::invalid_argument);

    for (int j=1; j<num_tags; ++j)
    {
      fscanf(file_ptr, "%s", buffer);
    }

    elm_dim = dimForMshTag(EMshTag(type_tag));

    if (elm_dim==0)
    {
      fscanf(file_ptr, "%d", &id_aux);
      --id_aux;
      mesh->getNode(id_aux)->setTag(physical);
    }
    else if (elm_dim == cell_dim)
    {
      cell = mesh->getCell(inc++);
      FEPIC_ASSERT(cell_msh_tag == type_tag, "Invalid cell or invalid mesh", std::runtime_error);
      for (int i=0; i< nodes_per_cell; ++i)
      {
        fscanf(file_ptr, "%d", &id_aux);
        cell->setNodeId(i, id_aux-1);
      }
      cell->setTag(physical);
      //mesh->pushCell(cell);
    }
    else
    {
      fgets(buffer, sizeof(buffer), file_ptr);
    }
  }// end for k


  this->timer.elapsed("readFileMsh(): read connectivity");
  // até aqui, apenas foi lido a conectividade
  //

  /* constroi as facets e Corners */
  if (mesh->qBuildAdjacency())
    mesh->buildAdjacency();
  else
  {
    fclose(file_ptr);
    return;
  }

  /*
  ___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__
  _|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___| */

  this->timer.restart();
  /* Procura por elementos no contorno que coincidam com as facets. Quando encontrados,
   * os labels são propagados para as respectivas facet que por sua vez propagam os
   * labels para seus nós. Os labels só são propagados se eles não foram definidos
   * nos nós anteriormente.
   */
  fseek (file_ptr , elems_file_pos , SEEK_SET );
  fscanf(file_ptr, "%d", &num_elms);

  int n_nodes               = mesh->nodesPerCell();
  //int n_vertices_per_facet  = mesh->verticesPerFacet();
  int n_nodes_per_facet     = mesh->nodesPerFacet();
  //int n_vertices_per_corner = mesh->verticesPerCorner();
  int n_nodes_per_corner    = mesh->nodesPerCorner();


  int* nodes  = new int[n_nodes_per_facet];      // facet nodes
  //int* vtcs   = new int[n_vertices_per_facet];
  int* bnodes = new int[n_nodes_per_corner];     // corner nodes
  //int* bvtcs  = new int[n_vertices_per_corner];
  int facet_id;
  int corner_id;

  for (int k=0; k < num_elms; ++k)
  {

    fscanf(file_ptr, "%d %d %d %d", &elem_number, &type_tag, &num_tags, &physical);

    //// sincronização
    //FEPIC_ASSERT(elem_number==k+1, "invalid file format", std::invalid_argument);

    for (int j=1; j<num_tags; ++j)
    {
      fscanf(file_ptr, "%s", buffer);
    }

    elm_dim = dimForMshTag(EMshTag(type_tag));

    if ((elm_dim == 0) && (cell_dim!=2))
    {
      fscanf(file_ptr, "%d", &id_aux);
    }
    else if (elm_dim == cell_dim-1) // facets
    {
      for (int i=0; i<n_nodes_per_facet; ++i)
      {
        fscanf(file_ptr, "%d", &nodes[i]);
        --nodes[i];
        if (mesh->getNode(nodes[i])->getTag() == 0)
          mesh->getNode(nodes[i])->setTag(physical);
      }
      //std::copy( nodes, nodes + n_vertices_per_facet, vtcs );
      if (cell_dim > 1)
      {
        if (mesh->getFacetIdFromVertices(nodes, facet_id))
          mesh->getFacet(abs(facet_id))->setTag(physical); //std::cout << (++TESTE) << std::endl;
        else
        {
          printf("WARNING: INVALID FACET IN INPUT MESH! vtcs: ");
          for (int zz = 0; zz < n_nodes_per_facet; ++zz)
            printf("%d ", nodes[zz]);
          printf("\n");
        }
      }
    }
    else if (elm_dim == cell_dim-2) // corners
    {
      for (int i=0; i<n_nodes_per_corner; ++i)
      {
        fscanf(file_ptr, "%d", &bnodes[i]);
        --bnodes[i];
        if (mesh->getNode(bnodes[i])->getTag() == 0)
          mesh->getNode(bnodes[i])->setTag(physical);
      }
      //std::copy( bnodes, bnodes + n_vertices_per_corner, bvtcs );
      if (cell_dim>2)
      {
        if (mesh->getCornerIdFromVertices(bnodes, corner_id))
        {
          mesh->getCorner(abs(corner_id))->setTag(physical); //std::cout << (++TESTE) << std::endl;
        }
        else if (mesh->isVertex(mesh->getNode(bnodes[0])) ) // if is vertex
            printf("WARNING: INVALID CORNER IN INPUT MESH!\n");
      }
    }
    else
    {
      for (int i=0; i<n_nodes; ++i)
      {
        fscanf(file_ptr, "%d", &id_aux);
        --id_aux;

        if ((mesh->getNode(id_aux)->getTag()) == 0)
          mesh->getNode(id_aux)->setTag(physical);
      }
    }


  }

  this->timer.elapsed("readFileMsh(): search for boundary elements");

  if (mesh->cellDim()>2)
  {
    const int n_corners_per_facet = mesh->numCornersPerFacet();
    // int facet_facets[n_corners_per_facet];
    int *facet_facets = new int [n_corners_per_facet];
    Facet const* facet;
    Corner* corner;
    for (int i = 0; i < mesh->numFacetsTotal(); ++i)
    {
      facet = mesh->getFacet(i);
      if (facet->disabled())
        continue;
      mesh->getCell(facet->getIncidCell())->getFacetCornersId(facet->getPosition(), facet_facets);

      for (int j = 0; j < n_corners_per_facet; ++j)
      {
        corner = mesh->getCorner(facet_facets[j]);
        if (corner->getTag() == 0)
          corner->setTag(facet->getTag());
      }
    }
    
    delete [] facet_facets;
    facet_facets = NULL;
  }

  fclose(file_ptr);
  //File.close();

  delete [] nodes;
  //delete [] vtcs;
  delete [] bnodes;
  //delete [] bvtcs;

  mesh->timer.addItems(this->timer);
}




