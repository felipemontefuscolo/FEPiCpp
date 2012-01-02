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

#ifndef FEPIC_DOF_HANDLER_HPP
#define FEPIC_DOF_HANDLER_HPP

#include "../shapefunctions/shape_declarations.hpp"
#include "var_dof.hpp"
#include <vector>
#include <set>
#include <string>


class DofHandler
{
  typedef std::vector<int>                   DataContainer;
  typedef std::vector<VarDofs>               DofsContainer;
  typedef Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixBool;
public:  
  
  DofHandler(Mesh* mesh=NULL, float gf=0.07) : _mesh_ptr(mesh), _grow_factor(gf) {}
  
  void setMesh(Mesh * mesh) {_mesh_ptr=mesh;};
  void addVariable(const char* var_name, ShapeFunction *sf, int dim=1);
  void addVariable(const char* var_name, int ndpv, int ndpr, int ndpf, int ndpc);
  VarDofs const& getVariable(int i) const {return _vars[i];}
  VarDofs & getVariable(int i) {return _vars[i];}
  
  void setVariablesRelationship(bool const* v);
  
  void printSparsityMatlab(std::string matname, std::string filename = "");
  void getSparsityTable(std::vector<std::set<int> > & table);
  void getMetisCSRfromTable(int *, int *, std::vector<std::set<int> > const& table);
  void getCSR_adjacency(std::vector<int> &adjncy, std::vector<int> &xadj);

  void metisRenumber();
  void boostMinimumDegreeRenumber();
  void boostCuthillMcKeeRenumber();
  void CuthillMcKeeRenumber();
  
  int numVars() const {return _vars.size();};
  int numDofs() const;
  void SetUp();
  int const* data() const {return _data.data();};
  int* data() {return _data.data();};
  int totalSize() const {return _data.size();};
  
private:
  Mesh*         _mesh_ptr;
  float         _grow_factor;
  MatrixBool    _relations;
  DofsContainer _vars;
  DataContainer _data; // single block with dofs
};



#endif
