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

/* please add posfix */



template<class _Traits>
class _MeshIoVtk
{

public:
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::MeshT  MeshT;
  typedef typename _Traits::PointT  PointT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoVtk() : _filenumVtk(0), _add_scalar_vtk_n_calls(0) {};
  _MeshIoVtk(_MeshIoVtk const&) {};

public:

  void writeVtk(std::string outname, bool flinear);
  template<class T>
  void addScalarVtk(const char* nome_var, T& scalar, int num_pts);
  template<class T>
  void addVectorVtk(const char* nome_var, T& arrayos, int dim, int num_pts);
  void addPointTagVtk(const char* nome_var); // para debug
  void addPointHalfVtk(const char* nome_var);  // para debug

  /** Imprime as coordenadas de um ponto em um stream dado.
   * @param o ponto
  *  @param o stream onde se vai imprimir.
  *  @param space Espaço entre a impressão de cada dimensão da coordenada.
  */
  static void printPointVtk(PointT const& p, std::ostream &o, int space = 22)
  {
    switch (_Traits::spacedim) {
      case 1:
      {
        o << std::left << std::setw(space) << p.getCoord(0) << " 0.0 0.0";
        break;
      }
      case 2:
      {
        o << std::left << std::setw(space) << p.getCoord(0)
                       << std::setw(space) << p.getCoord(1) << " 0.0";
        break;
      }
      case 3:
      {
        o << std::left << std::setw(space) << p.getCoord(0)
                       << std::setw(space) << p.getCoord(1)
                       << std::setw(space) << p.getCoord(2);
        break;
      }
      default:
      {
        std::cout << "Error: invalid dimension\n";
        break;
      }
    }
  }


  /// Tetrahedron Cell printer
  template<class Cell_T>
  static void printCellVtk(Cell_T const& cell, std::ostream &o, int order=1,
  typename std::enable_if<std::is_same<Simplex<1>,  typename Cell_T::PolytopeT>::value ||
                          std::is_same<Hypercube<1>,typename Cell_T::PolytopeT>::value ||
                          std::is_same<Polytope<1>, typename Cell_T::PolytopeT>::value >::type* = NULL)
  {
    FEPIC_CHECK( (order==1) || (order==cell.getOrder()), "invalid order", std::invalid_argument);

    if (order<=1)
      o << "2 " << cell.getNodeIdx(0) << " " << cell.getNodeIdx(1);
    else
    {
      o << "2 " << cell.getNodeIdx(0) << " " << cell.getNodeIdx(2);
      for (int i = 0; i < order-2; i++)
      {
        o << std::endl;
        o << "2 " << cell.getNodeIdx(2+i) << " " << cell.getNodeIdx(3+i);
      }
      o << std::endl;
      o << "2 " << cell.getNodeIdx(order) << " " << cell.getNodeIdx(1);
    }
  }

  /// Tetrahedron Cell printer
  template<class Cell_T>
  static void printCellVtk(Cell_T const& cell, std::ostream &o, int order=1,
  typename std::enable_if<std::is_same<Simplex<2>, typename Cell_T::PolytopeT>::value>::type* = NULL)
  {
    FEPIC_CHECK( (order==1) || (order==cell.getOrder()), "invalid order", std::invalid_argument);
    matrixi const& minimesh = Cell_T::getMinimesh(order);

    o <<"3 "<<  cell.getNodeIdx(minimesh[0][0]) << " " <<
                cell.getNodeIdx(minimesh[0][1]) << " " <<
                cell.getNodeIdx(minimesh[0][2]);
    for (int i = 1; i < minimesh.size(); ++i)
    {
      o << std::endl;
      o <<"3 "<<  cell.getNodeIdx(minimesh[i][0]) << " " <<
                  cell.getNodeIdx(minimesh[i][1]) << " " <<
                  cell.getNodeIdx(minimesh[i][2]);
    }
  }


  /// Tetrahedron Cell printer
  template<class Cell_T>
  static void printCellVtk(Cell_T const& cell, std::ostream &o, int order=1,
  typename std::enable_if<std::is_same<Simplex<3>, typename Cell_T::PolytopeT>::value>::type* = NULL)
  {
    FEPIC_CHECK( (order==1) || (order==cell.getOrder()), "invalid order", std::invalid_argument);
    matrixi const& minimesh = Cell_T::getMinimesh(order);

    o <<"4 "<< cell.getNodeIdx(minimesh[0][0]) << " "
            << cell.getNodeIdx(minimesh[0][1]) << " "
            << cell.getNodeIdx(minimesh[0][2]) << " "
            << cell.getNodeIdx(minimesh[0][3]);
    for (int i = 1; i < minimesh.size(); ++i)
    {
      o << std::endl;
      o <<"4 "<< cell.getNodeIdx(minimesh[i][0]) << " "
              << cell.getNodeIdx(minimesh[i][1]) << " "
              << cell.getNodeIdx(minimesh[i][2]) << " "
              << cell.getNodeIdx(minimesh[i][3]);
    }

  }


  /* ---- Edge ----*/
  template<class Cell_T>
  static int getTagVtk(typename std::enable_if<std::is_same<Simplex<1>,  typename Cell_T::PolytopeT>::value ||
                                               std::is_same<Hypercube<1>,typename Cell_T::PolytopeT>::value ||
                                               std::is_same<Polytope<1>, typename Cell_T::PolytopeT>::value >::type* = NULL)
  {return 3; }


  // triangle
  template<class Cell_T>
  static int getTagVtk(typename std::enable_if<std::is_same<Simplex<2>, typename Cell_T::PolytopeT>::value>::type* = NULL)
  {return 5; }

  // tetrahedron
  template<class Cell_T>
  static int getTagVtk(typename std::enable_if<std::is_same<Simplex<3>, typename Cell_T::PolytopeT>::value>::type* = NULL)
  {return 10; }




protected:
  int _filenumVtk;
  int _add_scalar_vtk_n_calls;
};


/** Print mesh in vtk file format.
 * @param flinear if the mesh has higher order elements and this flag is false, than
 * the mesh print each higher order cell as a compound of linear cells, otherwise the
 * cells are printed as if lienar cells, disregarding high order nodes..
 */
template<class _Traits>
void _MeshIoVtk<_Traits>::writeVtk(std::string outname = "", bool flinear = false)
{
  /*
  * NOTA: TODOS os pontos são impressos. APENAS AS CÉLULAS VIVAS
  * são impressas.
  */

  int  order  = flinear ? 1 : THIS->_order;
  int  nc_mm  = CellT::getNumSubdivisions(order);    // num de mini-células por mini-malha
  int ncells = THIS->getNumCells() * nc_mm;

  THIS->_add_scalar_vtk_n_calls=0;

  std::string ss = outname=="" ? THIS->_popNextName(this->_filenumVtk, ".vtk") : outname;
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
    _MeshIoVtk::printPointVtk(*pit, Fout);
    Fout << std::endl;
  }

  /* imprimindo células */
  Fout << std::endl;
  Fout << "CELLS " << ncells << " " << ncells*(1+CellT::n_vertices) << std::endl;
  auto cit = THIS->_cellL.begin();
  for (auto cellend=THIS->_cellL.end(); cit != cellend ; ++cit)
  {
    if(1)//if(!cit->disabled())
    {
      _MeshIoVtk::printCellVtk(*cit, Fout, order);
      Fout << std::endl;
    }
  }
  int type = _MeshIoVtk::getTagVtk<CellT>();
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
void _MeshIoVtk<_Traits>::addScalarVtk(const char* nome_var, T& scalar, int num_pts)
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

  for (int i=0; i<num_pts; ++i)
    Fout << scalar[i] << std::endl;

  Fout << std::endl << std::endl;

  Fout.close();
};


template<class _Traits>
void _MeshIoVtk<_Traits>::addPointTagVtk(const char* nome_var="node labels")
{

  int num_pts = THIS->getNumNodes();

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

  for (int i=0; i<num_pts; ++i) {
    Fout << (THIS->getNode(i)->getTag()) << std::endl;
  };

  Fout << std::endl << std::endl;

  Fout.close();
};


template<class _Traits>
void _MeshIoVtk<_Traits>::addPointHalfVtk(const char* nome_var="node labels")
{

  int num_pts = THIS->getNumNodes();

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

  for (int i=0; i<num_pts; ++i) {
    Fout << (THIS->getNode(i)->getHalf()->getIDCell()) << std::endl;
  };

  Fout << std::endl << std::endl;

  Fout.close();
};





#undef THIS
#undef CONST_THIS

#endif
