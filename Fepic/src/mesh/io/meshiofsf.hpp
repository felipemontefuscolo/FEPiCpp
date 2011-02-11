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
  typedef typename _Traits::HalfT  HalfT;
  typedef typename _Traits::PointT  PointT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoFsf() : _filenumFsf(0), _add_scalar_fsf_n_calls(0) {};
  _MeshIoFsf(_MeshIoFsf const&) {};

public:

  void readFileFsf(const char* filename);
  void writeFsf(std::string outname);
  template<class T>
  void addScalarFsf(const char* nome_var, T& scalar, int num_pts);
  template<class T>
  void addVectorFsf(const char* nome_var, T& arrayos, int dim, int num_pts);
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
    
    o << cell.getTag() << " " << cell.getFlags();
    
    for (int i = 0; i < cell.n_borders; i++)
    {   
      auto h = cell.getHalf(i);
      o<<" "<< h->getIncidCell() <<" "<< h->getPosition() <<" "<< h->getAnchor();
    }
  }
 
  static void printHalfsFsf(HalfT const& hl, std::ostream &o)
  {
    o << hl.getIncidCell() << " " << hl.getPosition() << " " << hl.getAnchor() << " " << hl.getTag() << " " << hl.getFlags();
  }
  
protected:  
  int _filenumFsf;
  int _add_scalar_fsf_n_calls;
};

template<class _Traits>
void _MeshIoFsf<_Traits>::readFileFsf(const char* filename)
{
  
  THIS->_registerFile(filename, ".fsf");

  std::ifstream File(filename);

  int     type_aux;
  char    lixo[256];
  int     order;
  size_t  pos;
  std::string s;
  
  const char* keys[] = {"POINTS", "CELLS", "HALFLS"};
  int keys_size = std::extent<decltype(keys)>::value;

  pos = search_word(File, "CELLTYPE");
  FEPIC_ASSERT(pos != static_cast<size_t>(-1), "invalid file format", std::invalid_argument);
  File >> s;
  FEPIC_ASSERT(CellPolytopeT::name() == s, "invalid mesh cell type", std::invalid_argument);
  
  pos = search_word(File, "MESHORDER");
  FEPIC_ASSERT(pos != static_cast<size_t>(-1), "invalid file format", std::invalid_argument);
  File >> order;
  FEPIC_ASSERT((order > 0) && (order < MAX_CELLS_ORDER), "invalid or not supported order", std::invalid_argument);
  

  /*
   *      READING POINTS
   * 
   * */
  pos = search_word(File, "POINTS");
  FEPIC_ASSERT(pos != static_cast<size_t>(-1), "invalid file format: can not find POINTS", std::invalid_argument);
  int n_pts;
  int tag, flags;
  File >> n_pts;
  
  THIS->_pointL.resize(n_pts);
  
  int num;
  double coord[3];
  for (int i = 0; i < n_pts; ++i)
  {
    File >> coord[0] >> coord[1] >> coord[2] >> tag >> flags;

    THIS->_pointL[i].setCoord(coord);
    THIS->_pointL[i].setTag(tag);
    THIS->_pointL[i].setFlags(flags);
  }
  File >> s;
  FEPIC_ASSERT(s == std::string("END_POINTS"), "can not find the key END_CELLS", std::invalid_argument);


  /*
   *      READING CELLS
   * 
   * */
  pos = search_word(File, "CELLS");
  FEPIC_ASSERT(pos != static_cast<size_t>(-1), "invalid file format: can not find CELLS", std::invalid_argument);
  int n_cells, nodeid;
  File >> n_cells;
  
  THIS->_cellL.resize(n_cells);
  
  int n_nodes_per_cell = ::numNodes<CellPolytopeT>(order);
  int n_borders = CellT::n_borders;
  int incid_cell;
  int position;
  int anchor;
  num=0;
  for (int i = 0; i < n_cells; i++)
  {
    THIS->_cellL[i].setOrder(order);
    for (int n = 0; n < n_nodes_per_cell; ++n)
    {
      File >> nodeid; // node
      THIS->_cellL[i].setNode(n, nodeid);
    }
    File >> tag >> flags;
    THIS->_cellL[i].setTag(tag);
    THIS->_cellL[i].setFlags(flags);
    
    for (int b = 0; b < n_borders; ++b)
    {
      File >> incid_cell >> position >> anchor;
      THIS->_cellL[i].getHalf(b)->setIncidCell(incid_cell);
      THIS->_cellL[i].getHalf(b)->setPosition(position);
      THIS->_cellL[i].getHalf(b)->setAnchor(anchor);
    }
    
  }
  File >> s;
  FEPIC_ASSERT(s == std::string("END_CELLS"), "can not find the key END_CELLS", std::invalid_argument);
    


  /*
   *      READING HALFLS
   * 
   * */
  pos = search_word(File, "HALFLS");
  FEPIC_ASSERT(pos != static_cast<size_t>(-1), "invalid file format: can not find HALFLS", std::invalid_argument);
  int n_halfs;
  File >> n_halfs;
  
  THIS->_halfL.resize(n_halfs);
  
  num=0;
  for (int i = 0; i < n_halfs; i++)
  {
    File >> incid_cell >> position >> anchor >> tag >> flags;
    THIS->_halfL[i].setIncidCell(incid_cell);
    THIS->_halfL[i].setPosition(position);
    THIS->_halfL[i].setAnchor(anchor);
    THIS->_halfL[i].setTag(tag);
    THIS->_halfL[i].setFlags(flags);
  }
  File >> s;
  FEPIC_ASSERT(s == std::string("END_HALFLS"), "can not find the key END_CELLS", std::invalid_argument);
  
  
  File.close();
  
}



template<class _Traits>
void _MeshIoFsf<_Traits>::writeFsf(std::string outname = "")
{
  /*
  * NOTA: TODOS os pontos são impressos. APENAS AS CÉLULAS VIVAS
  * são impressas.
  */
  int order  = THIS->getOrder();
  int ncells = THIS->getNumCellsTotal();
  int nhalfs= THIS->getNumHalfTotal();

  THIS->_add_scalar_fsf_n_calls=0;

  
	std::string ss = outname=="" ? THIS->_popNextName(this->_filenumFsf, ".fsf") : outname;
  ++_filenumFsf;

  std::ofstream Fout( ss.data() );

  Fout << "# FEPiC++ file format : alpha version\n\n"
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
  Fout << "END_POINTS\n\n";

  /* imprimindo células */
  Fout << "\n# nodes | tag | flags | [half1(c, position, anc) | half2(...) | ...]\n";
  Fout << "CELLS " << ncells << "\n";
  auto cit = THIS->_cellL.begin();
  for (auto cellend=THIS->_cellL.end(); cit != cellend ; ++cit)
  {
    _MeshIoFsf::printCellFsf(*cit, Fout);
    Fout << std::endl;
  }
  Fout << "END_CELLS\n\n";

  Fout << "#incident cell | position | anchor | tag | flags\n";
  Fout << "HALFLS " << nhalfs << "\n";
  for (auto hit = THIS->_halfL.begin(); hit != THIS->_halfL.end(); ++hit)
  {
    _MeshIoFsf::printHalfsFsf(*hit, Fout);
    Fout << std::endl;
  }
  Fout << "END_HALFLS\n\n";

  Fout << std::endl;

  Fout.close();  
  
}

#undef THIS
#undef CONST_THIS

#endif


