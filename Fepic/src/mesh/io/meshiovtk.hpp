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


class MeshIoVtk : public _MeshNameHandler
{

public:

  MeshIoVtk() :_filenumVtk(0), _add_node_scalar_vtk_n_calls(0), _mesh(NULL) {};

  explicit MeshIoVtk(Mesh const* mesh, int filenum=0) : _filenumVtk(filenum), _add_node_scalar_vtk_n_calls(0)
  {
    attachMesh(mesh);
  }

  MeshIoVtk(MeshIoVtk const&) {};

  typedef void (MeshIoVtk::*CprinterMemFunPtr)(int const*, FILE*) const;
  typedef void (MeshIoVtk::*PprinterMemFunPtr)(Point const*, FILE*) const;

  void attachMesh(Mesh const* mesh);

  std::string outputFileName()
  {
    return this->_popNextName(this->_filenumVtk, ".vtk");
  }
  
  void writeVtk(std::string outname = "");
  void addNodeScalarVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addCellScalarVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addNodeIntVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void addCellIntVtk(const char* nome_var, DefaultGetDataVtk const& data);
  void printPointTagVtk(const char* nome_var="node_labels"); //  debug
  void printPointIcellVtk(const char* nome_var="pt_icell"); //  debug
  void printPointPositionVtk(const char* nome_var="pt_ic_pos"); //  debug


  void _printPointVtk_1d(Point const* p, FILE *fp) const;
  void _printPointVtk_2d(Point const* p, FILE *fp) const;
  void _printPointVtk_3d(Point const* p, FILE *fp) const;


  void _printCellVtk_Edge2(int const* ids, FILE *fp) const;
  
  void _printCellVtk_Edge3(int const* ids, FILE *fp) const;

  void _printCellVtk_Triangle3(int const* ids, FILE *fp) const;

  void _printCellVtk_Triangle6(int const* ids, FILE *fp) const;

  void _printCellVtk_Quadrangle4(int const* ids, FILE *fp) const;

  void _printCellVtk_Quadrangle8(int const* ids, FILE *fp) const;

  void _printCellVtk_Quadrangle9(int const* ids, FILE *fp) const;

  void _printCellVtk_Tetrahedron4(int const* ids, FILE *fp) const;

  void _printCellVtk_Tetrahedron10(int const* ids, FILE *fp) const;

  void _printCellVtk_Hexahedron8(int const* ids, FILE *fp) const;

  void _printCellVtk_Hexahedron20(int const* ids, FILE *fp) const;

  void _printCellVtk_Hexahedron27(int const* ids, FILE *fp) const;


  static int getVtkTag(ECellType type);

  static int getNumDivisions(ECellType type);

protected:
  int  _filenumVtk;
  int  _add_node_scalar_vtk_n_calls;
  int  _add_cell_scalar_vtk_n_calls;
  int  _spacedim;
  ECellType _fep_tag;
  Mesh const* _mesh;
  CprinterMemFunPtr _c_printer;
  PprinterMemFunPtr _p_printer;

};




#endif
