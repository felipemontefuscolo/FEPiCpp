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

#ifndef FEPIC_MESHIOVTK_HPP
#define FEPIC_MESHIOVTK_HPP

#include <cstdio>
#include "../../util/assert.hpp"
#include "../mesh.hpp"
#include "meshnamehandler.hpp"
#include "../enums.hpp"

class DefaultGetDataVtk
{
public:

  /** @brief Returns data for contiguous id.
   *  @warning I said for contiguous id.
   *  @note warning: for contiguous id.
   */  

  DefaultGetDataVtk(Real * r = NULL, int * i = NULL) : data_r(r), data_i(i) {}
  //DefaultGetDataVtk(int * i = NULL, Real * r = NULL) : data_r(r), data_i(i) {}
  
  virtual Real get_data_r(int id) const
  {
    return data_r[id];
  }
  virtual int get_data_i(int id) const
  {
    return data_i[id];
  }
  
  virtual ~DefaultGetDataVtk() {}
  
protected:
  Real * data_r;
  int  * data_i;
};


class MeshIoVtk : public iMeshNameHandler
{

public:

  MeshIoVtk() :m_filenumVtk(0), m_add_node_scalar_vtk_n_calls(0), m_mesh(NULL) {};

  explicit MeshIoVtk(Mesh const* mesh, int filenum=0) : m_filenumVtk(filenum), m_add_node_scalar_vtk_n_calls(0)
  {
    attachMesh(mesh);
  }

  MeshIoVtk(MeshIoVtk const&) {};

  typedef void (MeshIoVtk::*CprinterMemFunPtr)(int const*, FILE*) const;
  typedef void (MeshIoVtk::*PprinterMemFunPtr)(Point const*, FILE*) const;

  void attachMesh(Mesh const* mesh);

  std::string outputFileName()
  {
    return this->fi_popNextName(this->m_filenumVtk, ".vtk");
  }
  
  void writeVtk(std::string outname = "");
  void addNodeScalarVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addCellScalarVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addNodeIntVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addCellIntVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void printPointTagVtk(const char* nome_var="node_labels"); //  debug
  void printPointIcellVtk(const char* nome_var="pt_icell"); //  debug
  void printPointPositionVtk(const char* nome_var="pt_ic_pos"); //  debug


  void fi_printPointVtk_1d(Point const* p, FILE *fp) const;
  void fi_printPointVtk_2d(Point const* p, FILE *fp) const;
  void fi_printPointVtk_3d(Point const* p, FILE *fp) const;


  void fi_printCellVtk_Edge2(int const* ids, FILE *fp) const;
  
  void fi_printCellVtk_Edge3(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Triangle3(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Triangle6(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Quadrangle4(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Quadrangle8(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Quadrangle9(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Tetrahedron4(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Tetrahedron10(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Hexahedron8(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Hexahedron20(int const* ids, FILE *fp) const;

  void fi_printCellVtk_Hexahedron27(int const* ids, FILE *fp) const;


  static int getVtkTag(ECellType type);

  static int getNumDivisions(ECellType type);

protected:
  int  m_filenumVtk;
  int  m_add_node_scalar_vtk_n_calls;
  int  m_add_cell_scalar_vtk_n_calls;
  int  m_spacedim;
  ECellType m_fep_tag;
  Mesh const* m_mesh;
  CprinterMemFunPtr m_c_printer;
  PprinterMemFunPtr m_p_printer;

};




#endif
