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

#ifndef FEPIC_IMESH_HPP
#define FEPIC_IMESH_HPP






template<class _Traits>
class iMesh : public _iMeshNameHandler,
              public _MeshIoMsh<_Traits>,
              public _MeshIoVtk<_Traits>,
              public _MeshIoFsf<_Traits>
{

public:

  friend class _MeshIoMsh<_Traits>;
  friend class _MeshIoVtk<_Traits>;
  friend class _MeshIoFsf<_Traits>;
  friend class _iMeshNameHandler;

  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::PointT PointT;

  typedef typename _Traits::HalfT  HalfT;
  typedef typename _Traits::HalfLT HalfLT;

  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

  typedef std::deque<CellT>     CellList;
  typedef std::deque<PointT>    PointList;
  typedef std::deque<HalfLT>    HalflList;

  typedef typename CellList::iterator  CellIterator;
  typedef typename PointList::iterator PointIterator;
  typedef typename HalflList::iterator HalflIterator;

  typedef typename CellList::const_iterator  CellConstIterator;
  typedef typename PointList::const_iterator PointConstIterator;
  typedef typename HalflList::const_iterator HalflConstIterator;

  typedef Eigen::Matrix<double, _Traits::spacedim, 1> VecT;


  iMesh()
  {
    _order = 1;
  };

  iMesh(iMesh const&) = delete;
  ~iMesh() = default;

  /* ------------------------------------------------------------
  * IO members
  * ------------------------------------------------------------*/

  void setFamilyFiles()
  {
    _is_family=true;
  }

  /* Em todas as funções de leitura de arquivo, DEVE ser chamado
     * a função setOrder() */

  //void readFileMsh(const char *filename);

  void buildAdjacency4face();
  void buildAdjacency4volume();

  void writeFileState();

        /*-------------------------------------------------------------*/
        /*-------------------------------------------------------------*/

  /** Imprime algumas informações da malha.
  */
  void printInfo(std::ostream &o = std::cout) const
  {
    //std::cout << "mesh file: " << _meshfile << std::endl;
    std::cout << "cell type: " << CellT::name() << std::endl;
    std::cout << "order:     " << this->_order << std::endl;
    std::cout << "# nodes:   " << getNumNodes() << std::endl;
    std::cout << "# cells:   " << getNumCells() << std::endl;
    std::cout << "# mhalfs:  " << getNumHalfls() << std::endl;
  }


  /** Altera a ordem da malha.
  *  @param order a nova ordem da malha
  *  @warning NÂO TERMINADO
  */
  void setOrder(int order)
  {
    if ((_order != order) && (order>0))
    {
      _order = order;
      for (auto it= _cellL.begin(), cend=_cellL.end(); it != cend; ++it)
        it->setOrder(_order);

      /* (re)cria os nós que são compartilhados entre as células. */
      _MeshMethods<iMesh<_Traits>, CellT::dim>::remodelCellsNodes(*this, order);
    }
  }

  /** NOT FOR USERS
  * Depois de mudada a ordem de cada célula da malha e atualizada a matriz que contém a numeração
  * local dos nós celulares, chama-se esta função. Ela (re)cria os nós que são compartilhados entre
  * as células.
  */
  void remodelCellsNodes4face(int n)
  {

  } // end remodelCellsNodes4face()

  /** NOT FOR USERS
  * Depois de mudada a ordem da malha, chama-se esta função. Ela recria os nós
  * e determina suas posições.
  */
  void remodelCellsNodes4volume(int n)
  {

  }

  /** Verifica se os nós passados formam alguma Mhalf já existente.
  * @param[in] nodes os nós.
  * @param[out] id o id da Mhalf.
  * @param[in] dead_mh bool indicando se as Mhalf (dead) incluem na pesquisa.
  * @note o critério de existência é se os nós são ciclicamente iguais aos da Mhalf.
  */
  bool theseVerticesFormAHalfl(vectorui const& vtx, uint &half_id, bool dead_mh = false)
  {
    // OTIMIZAR

    bool RET=false;
    uint nummhalf = this->getNumHalfLTotal();
    // IMPLEMENTAR ITERADOR
    for(uint i=0; i<nummhalf; ++i)
    {
      bool b = this->getHalfl(i)->hasTheseVertices(vtx, *this);

      if(b && (dead_mh || (!this->getHalfl(i)->disabled())) )
      {
        half_id = i;
        RET = true;
      }
    }
    return RET;
  }

  /** Retorna o n-ésimo nó da malha.
  */
  PointT* getNode(uint nth)
  {
    FEPIC_CHECK(nth<this->_pointL.size(), "invalid index", std::out_of_range);
    return &_pointL[nth];
  }

  const PointT* getNode(uint nth) const
  {
    FEPIC_CHECK(nth<this->_pointL.size(), "invalid index", std::out_of_range);
    return &_pointL[nth];
  }

  /** Retorna a n-ésima mhalf-edge/face (adivinha pelo tipo da malha)
  */
  HalfLT* getHalfl(uint nth)
  {
    FEPIC_CHECK(nth<this->_mhalfL.size(), "invalid index", std::out_of_range);
    return &_mhalfL[nth];
  }

  const HalfLT* getHalfl(uint nth) const
  {
    FEPIC_CHECK(nth<this->_mhalfL.size(), "invalid index", std::out_of_range);
    return &_mhalfL[nth];
  }


  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  CellT* getCell(uint nth)
  {
    FEPIC_CHECK(nth<this->_cellL.size(), "invalid index", std::out_of_range);
    return &_cellL[nth];
  }

  const CellT* getCell(uint nth) const
  {
    FEPIC_CHECK(nth<this->_cellL.size(), "invalid index", std::out_of_range);
    return &_cellL[nth];
  }

  /** Retorna na matriz X as coordenadas dos nós passados em map.
  *  As colunas de X correspondem a dimensão enquanto as linhas
  *  correspondem aos graus de liberdade.
  * */
  void getCellCoords(Eigen::MatrixXd &X, Eigen::VectorXi const& map)
  {
    X.resize(map.size(), _Traits::spacedim);

    for (int n = 0; n < map.size(); ++n)
      X.row(n) = this->getNode(map(n))->getCoord().transpose();
  }

  /** Adiciona uma célula e retorna seu id.
  */
  uint addCell(CellT const& C)
  {
    if (_dead_cells.empty())
    {
      _cellL.push_back(C);
      return _cellL.size()-1;
    }
    else
    {
      uint id = _dead_cells.back();
      _cellL.at(id) = C;
      _dead_cells.pop_back();
      return id;
    }
  }

  /** Adiciona um ponto e retorna seu id.
  */
  uint addPoint(PointT const& P)
  {
    if (_dead_points.empty())
    {
      _pointL.push_back(P);
      return _pointL.size()-1;
    }
    else
    {
      uint id = _dead_points.back();
      _pointL.at(id) = P;
      _dead_points.pop_back();
      return id;
    }
  }

  /** Adiciona uma half-xxxx
  *  @param h A half-xxxx a ser adicionada.
  *  @return A posição da half-xxxx na lista
  */
  uint addHalfl(HalfLT const& h)
  {
    if (_dead_mhalf.empty())
    {
      _mhalfL.push_back(h);
      return _mhalfL.size()-1;
    }
    else
    {
      uint id = _dead_mhalf.back();
      _mhalfL.at(id) = h;
      _dead_mhalf.pop_back();
      return id;
    }
  }

  /** Retorna o número células
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  uint getNumCells() const
  {
    return _cellL.size() - _dead_cells.size();
  }

  /** Retorna o número de células.
  * @note incluindo o/a(s) marcado/a(s) como killed.
  */
  uint getNumCellsTotal() const{
    return _cellL.size();
  }

  /** Retorna o número de nós.
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  uint getNumNodes() const
  {
    return _pointL.size() - _dead_points.size();
  }

  /** Retorna no número de nós.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  uint getNumNodesTotal() const{
    return _pointL.size();
  }

  /** Retorna número de mhalfs.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  uint getNumHalfls() const
  {
    return _mhalfL.size() - _dead_mhalf.size();
  }

  /** Retorna o número de mhalfs.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  uint getNumHalfLTotal() const
  {
    return _mhalfL.size();
  }

  /** apenas para 2D ainda.
  *  retorna o comprimento do contorno.
  */
  double getPerimeter()
  {
    auto it = _mhalfL.begin();
    vectorui vtx;
    double sum=0.;

    std::cout << _mhalfL.size() << std::endl;

    for (; it!= _mhalfL.end(); ++it)
      sum += it->HalfT::getLenght(*this);

    return sum;
  }

  /** Retorna a ordem da malha.
  */
  int getOrder() const
  {
    return _order;
  }

  /** Retorna um iterador apontando para o começo da lista de células.
  */
  CellIterator cellBegin()
  {
    return _cellL.begin();
  }

  /** Retorna um iterador apontando para o depois-do-final da lista de células.
  */
  CellIterator cellEnd()
  {
    return _cellL.end();
  }

  /** Retorna um iterador apontando para o começo da lista de pontos.
  */
  PointIterator pointBegin()
  {
    return _pointL.begin();
  }

  /** Retorna um iterador apontando para o depois-do-final da lista de pontos.
  */
  PointIterator pointEnd()
  {
    return _pointL.end();
  }

  /** Retorna um iterador apontando para o começo da lista de mhalfs.
  */
  HalflIterator mhalfBegin()
  {
    return _mhalfL.begin();
  }

  /** Retorna um iterador apontando para o depois-do-final da lista de mhalfs.
  */
  HalflIterator mhalfEnd()
  {
    return _mhalfL.end();
  }


  // propriedades da célula da malha e outros atributos auxiliares
  //matrixi edges_local_nodes; // face, volume
  //matrixi borders_local_nodes; // volume
  vectori opp_eln; // opposite edges_local_nodes
  matrixi opp_fln; // opposite borders_local_nodes
                          // ex: borders_local_nodes[f][opp_fln[i]] := é o que a célula oposta enxerga

  // entities
  CellList      _cellL;
  PointList     _pointL;
  HalflList     _mhalfL;

private:

  // deleted entities id's
  dequeui  _dead_cells;
  dequeui  _dead_points;
  dequeui  _dead_mhalf;

  int         _order;



};

#endif // IMESH
