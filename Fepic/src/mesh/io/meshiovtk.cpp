#include "meshiovtk.hpp"
#include <tr1/array>
#include "../../util/misc2.hpp"


void MeshIoVtk::attachMesh(Mesh const* mesh)
{
  FEPIC_ASSERT(mesh != NULL, "invalid mesh", std::invalid_argument);
  _mesh = mesh;
  _spacedim = _mesh->spaceDim();

  ECellType s = _mesh->cellType();
  
  switch (s)
  {
    case EDGE2:        _c_printer = &MeshIoVtk::_printCellVtk_Edge2;         break;
    case EDGE3:        _c_printer = &MeshIoVtk::_printCellVtk_Edge3;         break;
    case TRIANGLE3:    _c_printer = &MeshIoVtk::_printCellVtk_Triangle3;     break;
    case TRIANGLE6:    _c_printer = &MeshIoVtk::_printCellVtk_Triangle6;     break;
    case QUADRANGLE4:  _c_printer = &MeshIoVtk::_printCellVtk_Quadrangle4;   break;
    case QUADRANGLE8:  _c_printer = &MeshIoVtk::_printCellVtk_Quadrangle8;   break;
    case QUADRANGLE9:  _c_printer = &MeshIoVtk::_printCellVtk_Quadrangle9;   break;
    case TETRAHEDRON4: _c_printer = &MeshIoVtk::_printCellVtk_Tetrahedron4;  break;
    case TETRAHEDRON10:_c_printer = &MeshIoVtk::_printCellVtk_Tetrahedron10; break;
    case HEXAHEDRON8:  _c_printer = &MeshIoVtk::_printCellVtk_Hexahedron8;   break;
    case HEXAHEDRON20: _c_printer = &MeshIoVtk::_printCellVtk_Hexahedron20;  break;
    case HEXAHEDRON27: _c_printer = &MeshIoVtk::_printCellVtk_Hexahedron27;  break;
    default:
      FEPIC_CHECK(false, "invalid mesh cell type", std::invalid_argument);
  }
  

  switch (_spacedim)
  {
    case 1:  _p_printer = &MeshIoVtk::_printPointVtk_1d; break;
    case 2:  _p_printer = &MeshIoVtk::_printPointVtk_2d; break;
    case 3:  _p_printer = &MeshIoVtk::_printPointVtk_3d; break;
    default:
      FEPIC_CHECK(false, "invalid mesh dimension", std::invalid_argument);
  }
  
}


void MeshIoVtk::_printPointVtk_1d(Point const* p, FILE *fp) const
{
  fprintf(fp, "%f %d %d\n", static_cast<float>(p->getCoord(0)), 0, 0);
}
void MeshIoVtk::_printPointVtk_2d(Point const* p, FILE *fp) const
{
  fprintf(fp, "%f %f %d\n", static_cast<float>(p->getCoord()[0]), static_cast<float>(p->getCoord()[1]), 0);
}
void MeshIoVtk::_printPointVtk_3d(Point const* p, FILE *fp) const
{
  fprintf(fp, "%f %f %f\n", static_cast<float>(p->getCoord()[0]), static_cast<float>(p->getCoord()[1]), static_cast<float>(p->getCoord()[2]));
}


void MeshIoVtk::_printCellVtk_Edge2(int const* ids, FILE *fp) const
{
    fprintf(fp,"2 %d %d\n", ids[0], ids[1]);
}

void MeshIoVtk::_printCellVtk_Edge3(int const* ids, FILE *fp) const
{
    fprintf(fp,"2 %d %d\n", ids[0], ids[2]);
    fprintf(fp,"2 %d %d\n", ids[2], ids[1]);
}

void MeshIoVtk::_printCellVtk_Triangle3(int const* ids, FILE *fp) const
{
  fprintf(fp,"3 %d %d %d\n", *ids, ids[1], ids[2]);
}

void MeshIoVtk::_printCellVtk_Triangle6(int const* ids, FILE *fp) const
{
  fprintf(fp,"6 %d %d %d %d %d %d\n", *ids, ids[1], ids[2], ids[3], ids[4], ids[5]);
}

void MeshIoVtk::_printCellVtk_Quadrangle4(int const* ids, FILE *fp) const
{
  fprintf(fp,"4 %d %d %d %d\n", *ids, ids[1], ids[2], ids[3]);
}

void MeshIoVtk::_printCellVtk_Quadrangle8(int const* ids, FILE *fp) const
{
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", *ids, ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7]);
}

void MeshIoVtk::_printCellVtk_Quadrangle9(int const* ids, FILE *fp) const
{
  fprintf(fp,"4 %d %d %d %d\n", ids[8], ids[7], ids[0], ids[4]);
  fprintf(fp,"4 %d %d %d %d\n", ids[8], ids[4], ids[1], ids[5]);
  fprintf(fp,"4 %d %d %d %d\n", ids[8], ids[5], ids[2], ids[6]);
  fprintf(fp,"4 %d %d %d %d\n", ids[8], ids[6], ids[3], ids[7]);
}

void MeshIoVtk::_printCellVtk_Tetrahedron4(int const* ids, FILE *fp) const
{
  fprintf(fp,"4 %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
}

void MeshIoVtk::_printCellVtk_Tetrahedron10(int const* ids, FILE *fp) const
{
  fprintf(fp,"10 %d %d %d %d %d %d %d %d %d %d\n",
          ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7], ids[9], ids[8]);
}

void MeshIoVtk::_printCellVtk_Hexahedron8(int const* ids, FILE *fp) const
{
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7]);
}

void MeshIoVtk::_printCellVtk_Hexahedron20(int const* ids, FILE *fp) const
{
  fprintf(fp,"20 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
         ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7], ids[8], ids[11],
         ids[13], ids[9], ids[16], ids[18], ids[19], ids[17], ids[10], ids[12], ids[14], ids[15]);
}

void MeshIoVtk::_printCellVtk_Hexahedron27(int const* ids, FILE *fp) const
{
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[ 0],ids[ 8],ids[20],ids[ 9],ids[10],ids[21],ids[26],ids[22]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[10],ids[21],ids[26],ids[22],ids[ 4],ids[16],ids[25],ids[17]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[ 8],ids[ 1],ids[11],ids[20],ids[21],ids[12],ids[23],ids[26]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[21],ids[12],ids[23],ids[26],ids[16],ids[ 5],ids[18],ids[25]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[ 9],ids[20],ids[13],ids[ 3],ids[22],ids[26],ids[24],ids[15]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[22],ids[26],ids[24],ids[15],ids[17],ids[25],ids[19],ids[ 7]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[20],ids[11],ids[ 2],ids[13],ids[26],ids[23],ids[14],ids[24]);
  fprintf(fp,"8 %d %d %d %d %d %d %d %d\n", ids[26],ids[23],ids[14],ids[24],ids[25],ids[18],ids[ 6],ids[19]);
}


namespace _VtkTagsInitializers
{
  std::tr1::array<int , N_CELL_TYPES> _tags()
  {
    std::tr1::array<int , N_CELL_TYPES> tab;
    
    tab[log2_i32(EDGE2        ) ] = 3 ;
    tab[log2_i32(EDGE3        ) ] = 3 ;
    tab[log2_i32(TRIANGLE3    ) ] = 5 ;
    tab[log2_i32(TRIANGLE6    ) ] = 22;
    tab[log2_i32(QUADRANGLE4  ) ] = 9 ;
    tab[log2_i32(QUADRANGLE8  ) ] = 23;
    tab[log2_i32(QUADRANGLE9  ) ] = 9 ;
    tab[log2_i32(TETRAHEDRON4 ) ] = 10;
    tab[log2_i32(TETRAHEDRON10) ] = 24;
    tab[log2_i32(HEXAHEDRON8  ) ] = 12;
    tab[log2_i32(HEXAHEDRON20 ) ] = 25;
    tab[log2_i32(HEXAHEDRON27 ) ] = 12;
    
    return tab;
  }
  
  std::tr1::array<int , N_CELL_TYPES> _n_divs()
  {
    std::tr1::array<int , N_CELL_TYPES> tab;
    
    tab[log2_i32(EDGE2        ) ] = 1;
    tab[log2_i32(EDGE3        ) ] = 2;
    tab[log2_i32(TRIANGLE3    ) ] = 1;
    tab[log2_i32(TRIANGLE6    ) ] = 1;
    tab[log2_i32(QUADRANGLE4  ) ] = 1;
    tab[log2_i32(QUADRANGLE8  ) ] = 1;
    tab[log2_i32(QUADRANGLE9  ) ] = 4;
    tab[log2_i32(TETRAHEDRON4 ) ] = 1;
    tab[log2_i32(TETRAHEDRON10) ] = 1;
    tab[log2_i32(HEXAHEDRON8  ) ] = 1;
    tab[log2_i32(HEXAHEDRON20 ) ] = 1;
    tab[log2_i32(HEXAHEDRON27 ) ] = 8;
    
    return tab;
  }
  
}


int MeshIoVtk::getVtkTag(ECellType type)
{
  std::tr1::array<int , N_CELL_TYPES> vtk_tag = _VtkTagsInitializers::_tags();
      
  unsigned idx = log2_i32(type);
    
  FEPIC_CHECK(idx < N_CELL_TYPES, "invalid or not supported cell type", std::invalid_argument);
  
  return vtk_tag[idx];
}

int MeshIoVtk::getNumDivisions(ECellType type)
{
  std::tr1::array<int , N_CELL_TYPES> divs = _VtkTagsInitializers::_n_divs();
      
  unsigned idx = log2_i32(type);
    
  FEPIC_CHECK(idx < N_CELL_TYPES, "invalid or not supported cell type", std::invalid_argument);
  
  return divs[idx];
}


/** Print mesh in vtk file format.
 * @param flinear if the mesh has higher order elements and this flag is false, than
 * the mesh print each higher order cell as a compound of linear cells, otherwise the
 * cells are printed as if lienar cells, disregarding high order nodes..
 */
void MeshIoVtk::writeVtk(std::string outname)
{
  /*
  *  DISABLE entities are not printed.
  */

  this->_add_node_scalar_vtk_n_calls=0;
  this->_add_cell_scalar_vtk_n_calls=0;

  if (!(this->_sofn_already_called))
    this->setOutputFileName(outname.c_str());

  if (!strcmp(outname.c_str(), "")) // if equal
    outname = this->_popNextName(this->_filenumVtk, ".vtk");
  
  ++_filenumVtk;
  
  FILE *file_ptr = fopen(outname.c_str(), "w");

  fprintf(file_ptr, "# vtk DataFile Version 2.0\n"
                    "unstructured grid\n"
                    "ASCII\n"
                    "DATASET UNSTRUCTURED_GRID\n"
                    "POINTS %d float\n", _mesh->numNodes());

  Point const* point;
  int n_points_total = _mesh->numNodesTotal();
  
  /* printing the points */
  for (int i = 0; i < n_points_total; ++i)
  {
    point = _mesh->getNode(i);
    if (point->disabled())
      continue;
    CALL_MEMBER_FN(*this, _p_printer)(point, file_ptr);
  }
  
  
  int  n_cd   = getNumDivisions(_mesh->cellType());    // number os cell's division
  int  n_pseudo_cells = _mesh->numCells() * n_cd;
  int  n_cells_total  = _mesh->numCellsTotal();
  int  nodes[FEPIC_MAX_N_NODES_P_CELL];
  Cell const *cell;
  
  /* imprimindo células */
  fprintf(file_ptr,"\nCELLS %d %d\n", n_pseudo_cells, n_cd>1 ? n_pseudo_cells*(1+_mesh->verticesPerCell()) :
                                                               n_pseudo_cells*(1+_mesh->nodesPerCell()));
                                                               
  //printf("DEBUG n_cells_total = %d\n", n_cells_total);
  for (int k = 0; k < n_cells_total; ++k)
  {
    cell = _mesh->getCell(k);
    if (cell->disabled())
      continue;
    _mesh->getCellNodesContigId(cell,nodes);
    CALL_MEMBER_FN(*this, _c_printer)(nodes, file_ptr);
  }
  
  // printing types
  int type = MeshIoVtk::getVtkTag(_mesh->cellType());
  fprintf(file_ptr,"\nCELL_TYPES %d\n", n_pseudo_cells);
  for (int k = 0; k < n_cells_total; ++k)
  {
    cell = _mesh->getCell(k);
    if (cell->disabled())
      continue;
    for (int i = 0; i < n_cd; ++i)
    {
      fprintf(file_ptr, "%d\n", type);
    }
  }
  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}

void MeshIoVtk::addNodeScalarVtk(const char* nome_var, DefaultGetDataVtk const& data)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  const int num_pts_total = _mesh->numNodesTotal();
  const int num_pts       = _mesh->numNodes();
  if (_add_node_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"POINT_DATA %d\n\n", num_pts);
  }
  _add_node_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s float\nLOOKUP_TABLE default\n", nome_var);

  
  for (int i=0; i<num_pts_total; ++i)
    if (!(_mesh->getNode(i)->disabled()))
      fprintf(file_ptr,"%f\n", static_cast<float>(data.get_data_r(i)));

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}


void MeshIoVtk::addCellScalarVtk(const char* nome_var, DefaultGetDataVtk const& data)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  int  n_cd   = getNumDivisions(_mesh->cellType()); // num divisions
  const int num_cells_total = _mesh->numCellsTotal();
  const int num_cells = _mesh->numCells();
  if (_add_cell_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"CELL_DATA %d\n\n", num_cells*n_cd);
  }
  _add_cell_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s float\nLOOKUP_TABLE default\n", nome_var);

  for (int i=0; i<num_cells_total; ++i)
  {
    if (_mesh->getCell(i)->disabled())
      continue;
    for (int j = 0; j < n_cd; ++j)
      fprintf(file_ptr,"%f\n", static_cast<float>(data.get_data_r(i)));
  }

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}




// ======================= Interger version =====================================

void MeshIoVtk::addNodeIntVtk(const char* nome_var, DefaultGetDataVtk const& data)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  const int num_pts_total = _mesh->numNodesTotal();
  const int num_pts = _mesh->numNodes();
  
  if (_add_node_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"POINT_DATA %d\n\n", num_pts);
  }
  _add_node_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s int\nLOOKUP_TABLE default\n", nome_var);

  
  for (int i=0; i<num_pts_total; ++i)
    if (!(_mesh->getNode(i)->disabled()))
      fprintf(file_ptr,"%d\n", data.get_data_i(i));

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}


void MeshIoVtk::addCellIntVtk(const char* nome_var, DefaultGetDataVtk const& data)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  int  n_cd   = getNumDivisions(_mesh->cellType()); // num divisions
  const int num_cells_total = _mesh->numCellsTotal();
  const int num_cells = _mesh->numCells();
  
  if (_add_cell_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"CELL_DATA %d\n\n", num_cells*n_cd);
  }
  _add_cell_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s int\nLOOKUP_TABLE default\n", nome_var);

  for (int i=0; i<num_cells_total; ++i)
  {
    if (_mesh->getCell(i)->disabled())
      continue;
    for (int j = 0; j < n_cd; ++j)
      fprintf(file_ptr,"%d\n", data.get_data_i(i));
  }

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}



void MeshIoVtk::printPointTagVtk(const char* nome_var)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  const int num_pts = _mesh->numNodes();
  if (_add_node_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"POINT_DATA %d\n\n", num_pts);
  }
  _add_node_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s int\nLOOKUP_TABLE default\n", nome_var);

  int const num_pts_total = _mesh->numNodesTotal();
  Point const* point;
  for (int i=0; i<num_pts_total; ++i)
  {
    point = _mesh->getNode(i);
    if (point->disabled())
      continue;
    fprintf(file_ptr,"%d\n", point->getTag());
  }

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}


void MeshIoVtk::printPointIcellVtk(const char* nome_var)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  const int num_pts = _mesh->numNodes();
  if (_add_node_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"POINT_DATA %d\n\n", num_pts);
  }
  _add_node_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s int\nLOOKUP_TABLE default\n", nome_var);

  int const num_pts_total = _mesh->numNodesTotal();
  Point const* point;
  for (int i=0; i<num_pts_total; ++i)
  {
    point = _mesh->getNode(i);
    if (point->disabled())
      continue;
    fprintf(file_ptr,"%d\n", _mesh->getCellContigId(point->getIncidCell()));
  }

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}


void MeshIoVtk::printPointPositionVtk(const char* nome_var)
{

  std::string ss = this->_popNextName(this->_filenumVtk-1, ".vtk");

  FILE * file_ptr = fopen(ss.c_str(), "a");
  
  FEPIC_ASSERT(file_ptr != NULL, "could not open the file", std::runtime_error);
  
  const int num_pts = _mesh->numNodes();
  if (_add_node_scalar_vtk_n_calls==0)
  {
    fprintf(file_ptr,"POINT_DATA %d\n\n", num_pts);
  }
  _add_node_scalar_vtk_n_calls++;

  fprintf(file_ptr,"SCALARS %s int\nLOOKUP_TABLE default\n", nome_var);

  int const num_pts_total = _mesh->numNodesTotal();
  Point const* point;
  for (int i=0; i<num_pts_total; ++i)
  {
    point = _mesh->getNode(i);
    if (point->disabled())
      continue;
    fprintf(file_ptr,"%d\n", point->getPosition());
  }

  fprintf(file_ptr,"\n");

  fclose(file_ptr);
}



