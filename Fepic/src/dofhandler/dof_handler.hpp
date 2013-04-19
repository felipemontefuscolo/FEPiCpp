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
  
  DofHandler(Mesh* mesh=NULL, float gf=0.07) : m_mesh_ptr(mesh), m_grow_factor(gf) {}
  
  void setMesh(Mesh * mesh) {m_mesh_ptr=mesh;};
  
  /* Add a variable
   * @param var_name name of the variable
   * @param sf pointer to a ShapeFunction object. This object contains info about dofs.
   * @param ncomps number of components of the variable.
   * @param ntags number of different mesh elements tags which will be considered. 
   * @param tags only mesh elements with these tags will be considered.
   */ 
  void addVariable(const char* var_name, ShapeFunction *sf, int ncomps=1, int ntags=0, int const* tags=NULL);
  
  /*  Add a variable.
   * @param var_name name of the variable
   * @param ndpv number of dofs per vertice
   * @param ndpr number of dofs per corner (ignored for 1D and 2D cells)
   * @param ndpf number of dofs per facet (ignored for 1D cells)
   * @param ndpc number of dofs per cell
   * @param ntags number of different mesh elements tags which will be considered. 
   * @param tags only mesh elements with these tags will be considered.
   */ 
  void addVariable(const char* var_name, int ndpv, int ndpr, int ndpf, int ndpc, int ntags=0, int const* tags=NULL);
  VarDofs const& getVariable(int i) const {return m_vars[i];}
  VarDofs & getVariable(int i) {return m_vars[i];}
  
  void setVariablesRelationship(bool const* v);
  
  void printSparsityMatlab(std::string matname, std::string filename = "");
  void getSparsityTable(std::vector<std::set<int> > & table);
  void getMetisCSRfromTable(int *, int *, std::vector<std::set<int> > const& table);
  void getCSR_adjacency(std::vector<int> &adjncy, std::vector<int> &xadj);

  void metisRenumber();
  void boostMinimumDegreeRenumber();
  void boostCuthillMcKeeRenumber();
  void CuthillMcKeeRenumber();
  
  void removeDofsGaps();
  
  
  int numVars() const {return m_vars.size();};
  int numDofs() const;
  void SetUp();
  int const* data() const {return m_data.data();};
  int* data() {return m_data.data();};
  int totalSize() const {return m_data.size();};
  
private:
  Mesh*         m_mesh_ptr;
  float         m_grow_factor;
  MatrixBool    m_relations;
  DofsContainer m_vars;
  DataContainer m_data; // single block with dofs
};



#endif
