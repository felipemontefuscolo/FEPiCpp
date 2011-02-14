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


#ifndef FEPIC_HALFCORE_HPP
#define FEPIC_HALFCORE_HPP

template<class _Traits>
class _HalfCore
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<HalfT*>(this)
  #define CONST_THIS static_cast<const HalfT*>(this)
#endif

protected:
  _HalfCore() {};
  _HalfCore(_HalfCore const&) {};

  /** HalfT class has:
    * int _incid_cell :   ??;
    * int _position   :   ??;
    * int _anchor     :   ??;
    * enum {cell_id_limit= ??};
    * enum {position_limit=??};
    * enum {anchor_limit=  ??};
    * */

public:
  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::HalfT HalfT;
  typedef typename _Traits::CellT CellT;

  int getIncidCell() const
  {
    return CONST_THIS->_incid_cell;
  }

  int getPosition() const
  {
    return CONST_THIS->_position;
  }

  int getAnchor() const
  {
    return CONST_THIS->_anchor;
  }

  void setIncidCell(int cellid)
  {
    THIS->_incid_cell = cellid;
  }

  void setPosition(int pos)
  {
    THIS->_position = pos;
  }

  void setAnchor(int anchor)
  {
    THIS->_anchor = anchor;
  }

  void setCompleteId(int cellid, int pos, int anchor=0)
  {
    THIS->_incid_cell = cellid;
    THIS->_position = pos;
    THIS->_anchor = anchor;
  }

  bool isBoundary() const
  {
    return CONST_THIS->_anchor==-1 ? true : false;
  }

  /** Retorna um vetor com os índices dos vértices que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  Eigen::VectorXi getVertices(MeshT const& mesh) const
  {
    return mesh.getCell(this->getIncidCell()) // CELL->

           ->getBorderVertices(this->getPosition(), mesh);
  }


  /** Retorna um vetor com os índices dos nós que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  Eigen::VectorXi getNodes(MeshT const& mesh) const
  {
      return  mesh.getCell(this->getIncidCell()) // CELL->

              ->cell->getBorderNodes(this->getPosition(), mesh);
  }


  /** Verifica se esta HalfFace contém exatamente os vértices (índices) passados
  *  em v.
  *  @param v vetor com os índices dos vértices.
  *  @param mesh a malha na qual esta HalfFace está contida.
  *  @return true se os nós formam a HalfFace e false caso contrário.
  *  @note Os nós podem estar em qualquer orientação cíclica da original.
  */
  bool hasTheseVertices(Eigen::VectorXi const& v, MeshT const& mesh) const
  {
    FEPIC_CHECK(v.size()==CellT::n_vertices_per_border, "", std::invalid_argument);
    const CellT * const cell = mesh.getCell(CONST_THIS->getIncidCell());

    const Eigen::VectorXi vtx (cell->getBorderVertices(CONST_THIS->getPosition()));

    return arrayIsCyclicallyEqual(v, vtx);
  }

  /** Faz com que cada nó desta HalfEdgeLab aponte para ela.
  *  @param mesh a malha na qual a HalfEdgeLab está contida.
  */
  void broadcastHalf2Nodes(MeshT & mesh) const
  {
    const Eigen::VectorXi v (mesh.getCell( CONST_THIS->getIncidCell() )->getBorderNodes( CONST_THIS->getPosition()));

    for (int i = 0; i < v.size(); ++i)
      mesh.getNode(v[i])->setHalf(*CONST_THIS);

  }

    /** Atribui o rótulo desta HalfEdgeLab a seus nós.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem tag=0;
  */
  void broadcastTag2Nodes(MeshT & mesh, bool force=false) const
  {
    const Eigen::VectorXi v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition(), mesh));

    if (force)
      for (int i = 0; i < v.size(); ++i)
        mesh.getNode(v[i])->setTag(CONST_THIS->getTag());
    else
      for (int i = 0; i < v.size(); i++)
        if (mesh.getNode(v[i])->getTag() == 0)
          mesh.getNode(v[i])->setTag(CONST_THIS->getTag());
  }

#undef THIS
#undef CONST_THIS

};



#endif

