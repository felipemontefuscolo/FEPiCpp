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

#ifndef FEPIC_MESHIOFSF_HPP
#define FEPIC_MESHIOFSF_HPP


#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<typename _Traits::MeshT*>(this)
  #define CONST_THIS static_cast<const typename _Traits::MeshT*>(this)
#endif 


template<class _Traits>
class _MeshIoFsf
{
  
public:  
  typedef typename _Traits::CellT   CellT;
  typedef typename _Traits::MeshT   MeshT;
  typedef typename _Traits::HalfLT  HalfLT;
  typedef typename _Traits::PointT  PointT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoFsf() : _filenumFsf(0), _add_scalar_fsf_n_calls(0) {};
  _MeshIoFsf(_MeshIoFsf const&) {};

public:

  void readFileFsf(const char* filename);
  void writeFsf();
  template<class T>
  void addScalarFsf(const char* nome_var, T& scalar, uint num_pts);
  template<class T>
  void addVectorFsf(const char* nome_var, T& arrayos, int dim, uint num_pts);
  void addPointTagFsf(const char* nome_var); // para debug
  void addPointHalfFsf(const char* nome_var);  // para debug     




  /** Imprime as coordenadas de um ponto em um stream dado.
   * @param o ponto
  *  @param o stream onde se vai imprimir.
  *  @param space Espaço entre a impressão de cada dimensão da coordenada.
  */ 
  static void printPointFsf(PointT const& p, std::ostream &o, int space = 22)
  {
    switch (_Traits::spacedim) {
      case 1:
      {
        o << std::left << std::setw(space) << p.getCoord(0) << " 0.0 0.0 "
          << p.getTag() << " " << p.getFlags();
        break;
      }
      case 2:
      {
        o << std::left << std::setw(space) << p.getCoord(0)
                       << std::setw(space) << p.getCoord(1) << " 0.0 "
                       << p.getTag() << " " << p.getFlags();
        break;
      }
      case 3:
      {
        o << std::left << std::setw(space) << p.getCoord(0)
                       << std::setw(space) << p.getCoord(1)
                       << std::setw(space) << p.getCoord(2) << " "
                       << p.getTag() << " " << p.getFlags();
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
  static void printCellFsf(Cell_T const& cell, std::ostream &o,
  typename std::enable_if<std::is_same<Simplex<1>,  typename Cell_T::PolytopeT>::value ||
                          std::is_same<Hypercube<1>,typename Cell_T::PolytopeT>::value ||
                          std::is_same<Polytope<1>, typename Cell_T::PolytopeT>::value >::type* = NULL)
  {
    /* print:
     * 
     * _nodes | _order | _tag | _flags
     *  */
    
    for (int i = 0; i < cell.getNumNodes(); i++)
    {   
      o << cell.getNodeIdx(i) << " ";
    }
    
    o << cell.getOrder() << " " << cell.getTag() << " " << cell.getFlags();
  }

  /// Tetrahedron Cell printer
  template<class Cell_T>
  static void printCellFsf(Cell_T const& cell, std::ostream &o,
  typename std::enable_if<std::is_same<Simplex<2>, typename Cell_T::PolytopeT>::value ||
                          std::is_same<Simplex<3>, typename Cell_T::PolytopeT>::value >::type* = NULL)
  {
    /* print:
     * 
     * _nodes | _order | _tag | _flags | half1<c, ith, anc> | h2<> | ...
     *  */
    
    for (int i = 0; i < cell.getNumNodes(); i++)
    {   
      o << cell.getNodeIdx(i) << " ";
    }
    
    o << cell.getOrder() << " " << cell.getTag() << " " << cell.getFlags();
    
    for (int i = 0; i < cell.n_borders; i++)
    {   
      auto h = cell.getHalf(i);
      o<<" "<< h->getIncidCell() <<" "<< h->getPosition() <<" "<< h->getAnchor();
    }
  }
 
  static void printHalflsFsf(HalfLT const& hl, std::ostream &o)
  {
    o << hl.getIncidCell() << " " << hl.getPosition() << " " << hl.getAnchor() << " " << hl.getTag() << " " << hl.getFlags();
  }
  
protected:  
  uint _filenumFsf;
  uint _add_scalar_fsf_n_calls;
};

template<class _Traits>
void _MeshIoFsf<_Traits>::readFileFsf(const char* filename)
{
  
  
}



template<class _Traits>
void _MeshIoFsf<_Traits>::writeFsf()
{
  /*
  * NOTA: TODOS os pontos são impressos. APENAS AS CÉLULAS VIVAS
  * são impressas.
  */
  uint order  = THIS->getOrder();
  uint ncells = THIS->getNumCellsTotal();
  uint nhalfls= THIS->getNumHalfLTotal();

  THIS->_add_scalar_fsf_n_calls=0;

	std::string ss = THIS->_popNextName(this->_filenumFsf, ".fsf");
  ++_filenumFsf;

  std::ofstream Fout( ss.data() );

  Fout << "# FEPiC++ file format : alpha version\n\n"
       << "SPACEDIM " << _Traits::spacedim << "\n"
       << "CELLTYPE " << CellT::PolytopeT::name() << "\n"
       << "MESHORDER " << THIS->getOrder() << "\n\n"
       
       << "# x | y | z | tag | flags \n"
       << "POINTS " << THIS->_pointL.size() << "\n";

  /* imprimindo pontos */
  auto pit = THIS->_pointL.begin();
  for (auto pend=THIS->_pointL.end(); pit != pend ; ++pit)
  {
    _MeshIoFsf::printPointFsf(*pit, Fout);
    Fout << "\n";
  }

  /* imprimindo células */
  Fout << "\n #_nodes | _order | _tag | _flags | [half1(c, ith, anc) | half2(...) | ...]\n";
  Fout << "CELLS " << ncells << "\n";
  auto cit = THIS->_cellL.begin();
  for (auto cellend=THIS->_cellL.end(); cit != cellend ; ++cit)
  {
    _MeshIoFsf::printCellFsf(*cit, Fout);
    Fout << std::endl;
  }
  Fout << std::endl;

  Fout << "#incident cell | position | anchor | tag | flags\n";
  Fout << "HALFLS " << nhalfls << "\n";
  for (auto hit = THIS->_mhalfL.begin(); hit != THIS->_mhalfL.end(); ++hit)
  {
    _MeshIoFsf::printHalflsFsf(*hit, Fout);
    Fout << std::endl;
  }
  Fout << std::endl;

  Fout.close();  
  
}

#undef THIS
#undef CONST_THIS

#endif


