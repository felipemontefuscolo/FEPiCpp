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

template<class _Der>
class _HalfCore
{
public:

#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<_Der*>(this)
  #define CONST_THIS static_cast<const _Der*>(this)
#endif
  

  uint getIncidCell() const
  {
    //return THIS->_incid_cell;
    CONST_THIS->_incid_cell;
  }

  int getPosition() const
  {
    return CONST_THIS->_position-1;
  }  

  uint getAnchor() const
  {
    return CONST_THIS->_anchor;
  }
  
  void setIncidCell(uint cellid)
  {
    FEPIC_ASSERT(cellid<=_Der::cell_id_limit, "cell id limit exceeded");
    THIS->_incid_cell = cellid;
  }

  void setPosition(int pos)
  {
    FEPIC_ASSERT((pos>=-1)||(pos<=_Der::position_limit), "position limit exceeded");
    THIS->_position = pos;
  }
  
  void setAnchor(uint anchor)
  {
    FEPIC_ASSERT(anchor<=_Der::anchor_limit, "anchor limit exceeded");
    THIS->_anchor = anchor;
  }

  void setCompleteId(uint cellid, int pos, uint anchor=0)
  {
    FEPIC_ASSERT(cellid<=_Der::cell_id_limit, "cell id limit exceeded");
    FEPIC_ASSERT((pos>=-1)||(pos<=_Der::position_limit), "position limit exceeded");
    FEPIC_ASSERT(anchor<=_Der::anchor_limit, "anchor limit exceeded");
    THIS->_incid_cell = cellid;
    THIS->_position = pos;
    THIS->_anchor = anchor;
  }

  /** Retorna um vetor com os índices dos vértices que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  template<class Mesh>
  Fepic::vectorui getVertices(Mesh const& mesh) const
  {
    return mesh.getCell(this->getIncidCell()) // CELL->
                
           ->getBorderVertices(this->getPosition(), mesh);
  }


  /** Retorna um vetor com os índices dos nós que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  template<class Mesh>
  Fepic::vectorui getNodes(Mesh const& mesh) const
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
  template<class Mesh>
  bool hasTheseVertices(Fepic::vectorui const& v, Mesh const& mesh) const
  {
    const auto cell = mesh.getCell(this->getIncidCell());
    
    Fepic::vectorui vtx (cell->getBorderVertices(this->getPosition()));
    
    return arrayIsCyclicallyEqual(v, vtx);
  }
  
#undef THIS
#undef CONST_THIS
protected:
  _HalfCore() {};
  _HalfCore(_HalfCore const&) {};
  
  /** _Der class has:
    * uint _incid_cell : ??;
    * uint _position   : ??;
    * uint _anchor     : ??;
    * enum {cell_id_limit=??};
    * enum {position_limit=??};
    * enum {anchor_limit=??};
    * */
    
    
    
};



#endif

