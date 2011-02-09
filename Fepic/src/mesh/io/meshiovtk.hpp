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

#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<typename _Traits::MeshT*>(this)
  #define CONST_THIS static_cast<const typename _Traits::MeshT*>(this)
#endif 





template<class _Traits>
class _MeshIoVtk
{
  
public:  
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::MeshT  MeshT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoVtk() : _filenumVtk(0), _add_scalar_vtk_n_calls(0) {};
  _MeshIoVtk(_MeshIoVtk const&) {};
  
public:

  void writeVtk(bool flinear=false);
  template<class T>
  void addScalarVtk(const char* nome_var, T&& scalar, uint num_pts);
  template<class T>
  void addVectorVtk(const char* nome_var, T&& arrayos, int dim, uint num_pts);
  void addPointTagVtk(const char* nome_var); // para debug
  void addPointHalfVtk(const char* nome_var);  // para debug  
  
protected:  
  uint _filenumVtk;
  uint _add_scalar_vtk_n_calls;
};


/** Imprime os pontos e células no arquivo vtk.
 * @param flinear se true, força a malha ser impressa como linear, se false (default), a malha é impressa da ordem que ela é.
 *
 * ESTA DESCRIÇÃO PRECISA SER MELHORADA
 */
template<class _Traits>
void _MeshIoVtk<_Traits>::writeVtk(bool flinear)
{
  /*
  *
  * NOTA: TODOS os pontos são impressos. APENAS AS CÉLULAS VIVAS
  * são impressas.
  */

  int  order  = flinear ? 1 : THIS->_order;
  int  nc_mm  = CellT::getNumSubdivisions(order);    // num de mini-células por mini-malha
  uint ncells = THIS->getNumCells() * nc_mm;

  THIS->_add_scalar_vtk_n_calls=0;

	std::string ss = THIS->_popNextName(this->_filenumVtk, ".vtk");
  ++_filenumVtk;

  std::ofstream Fout( ss.data() );

  Fout << "# vtk DataFile Version 2.0" << std::endl
       << "unstructured grid" << std::endl
       << "ASCII" << std::endl
       << "DATASET UNSTRUCTURED_GRID" << std::endl
       << "POINTS " << THIS->_pointL.size() << " float" << std::endl;

  /* imprimindo pontos */
  auto pit = THIS->_pointL.begin();
  for (auto pend=THIS->_pointL.end(); pit != pend ; ++pit)
  {
    pit->printSelfVtk(Fout);
    Fout << std::endl;
  }

  /* imprimindo células */
  Fout << std::endl;
  Fout << "CELLS " << ncells << " " << ncells*(1+CellT::n_vertices) << std::endl;
  auto cit = THIS->_cellL.begin();
  for (auto cellend=THIS->_cellL.end(); cit != cellend ; ++cit)
  {
    if(!cit->disabled())
    {
      cit->printSelfVtk(Fout, order);
      Fout << std::endl;
    }
  }
  int type = CellT::getCellTypeVtk();
  Fout << std::endl; cit = THIS->_cellL.begin();
  Fout << "CELL_TYPES " << ncells << std::endl;
  for (auto cellend=THIS->_cellL.end(); cit != cellend ; ++cit)
  {
    if(!cit->disabled())
    {
			for (int i = 0; i < nc_mm; i++)
				Fout << type << std::endl;
    }
  }

  Fout << std::endl;


  Fout.close();
}

/** T: qualquer objeto chamável
 */
template<class _Traits>
template<class T>
void _MeshIoVtk<_Traits>::addScalarVtk(const char* nome_var, T&& scalar, uint num_pts)
{

	std::string ss = THIS->_popNextName(this->_filenumVtk, ".vtk");

  std::ofstream Fout;
  Fout.open( ss.data(), std::ios::app);

  Fout << std::setprecision(8);
  Fout << std::setiosflags(std::ios::showpoint);

  if (_add_scalar_vtk_n_calls==0)
  {
    Fout << "POINT_DATA " << num_pts << std::endl;
    Fout << std::endl;
  }
  _add_scalar_vtk_n_calls++;

	Fout << "SCALARS " << nome_var << " float"   << std::endl;
  Fout << "LOOKUP_TABLE default"            << std::endl;

  for (uint i=0; i<num_pts; ++i)
    Fout << scalar[i] << std::endl;

  Fout << std::endl << std::endl;

	Fout.close();
};


template<class _Traits>
void _MeshIoVtk<_Traits>::addPointTagVtk(const char* nome_var="node labels")
{

  uint num_pts = THIS->getNumNodes();

	std::string ss = THIS->_popNextName(_filenumVtk, ".vtk");

  std::ofstream Fout;
  Fout.open( ss.data(), std::ios::app);

  Fout << std::setprecision(8);
  Fout << std::setiosflags(std::ios::showpoint);

  if (_add_scalar_vtk_n_calls==0)
  {
    Fout << "POINT_DATA " << num_pts << std::endl;
    Fout << std::endl;
  }
  _add_scalar_vtk_n_calls++;

	Fout << "SCALARS " << nome_var << " int"   << std::endl;
  Fout << "LOOKUP_TABLE default"            << std::endl;

  for (uint i=0; i<num_pts; ++i) {
    Fout << (THIS->getNode(i)->getTag()) << std::endl;
  };

  Fout << std::endl << std::endl;

	Fout.close();
};


template<class _Traits>
void _MeshIoVtk<_Traits>::addPointHalfVtk(const char* nome_var="node labels")
{

  uint num_pts = THIS->getNumNodes();

	std::string ss = THIS->_popNextName(_filenumVtk, ".vtk");

  std::ofstream Fout;
  Fout.open( ss.data(), std::ios::app);

  Fout << std::setprecision(8);
  Fout << std::setiosflags(std::ios::showpoint);

  if (_add_scalar_vtk_n_calls==0)
  {
    Fout << "POINT_DATA " << num_pts << std::endl;
    Fout << std::endl;
  }
  _add_scalar_vtk_n_calls++;

	Fout << "SCALARS " << nome_var << " int"   << std::endl;
  Fout << "LOOKUP_TABLE default"            << std::endl;

  for (uint i=0; i<num_pts; ++i) {
    Fout << (THIS->getNode(i)->getHalf()->getIDCell()) << std::endl;
  };

  Fout << std::endl << std::endl;

	Fout.close();
};





#undef THIS
#undef CONST_THIS

#endif
