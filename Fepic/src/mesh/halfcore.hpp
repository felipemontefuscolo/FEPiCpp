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
    * uint _incid_cell :   ??;
    * uint _position   :   ??;
    * uint _anchor     :   ??;
    * enum {cell_id_limit= ??};
    * enum {position_limit=??};
    * enum {anchor_limit=  ??};
    * */
    
public:
  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::HalfT HalfT;
  typedef typename _Traits::CellT CellT;

  uint getIncidCell() const
  {
    return CONST_THIS->_incid_cell;
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
    FEPIC_CHECK(cellid<=HalfT::cell_id_limit, "cell id limit exceeded", std::out_of_range);
    THIS->_incid_cell = cellid;
  }

  void setPosition(int pos)
  {
    FEPIC_CHECK((pos>=-1)||(pos<=HalfT::position_limit), "position limit exceeded", std::out_of_range);
    THIS->_position = pos+1;
  }
  
  void setAnchor(uint anchor)
  {
    FEPIC_CHECK(anchor<=HalfT::anchor_limit, "anchor limit exceeded", std::out_of_range);
    THIS->_anchor = anchor;
  }

  void setCompleteId(uint cellid, int pos, uint anchor=0)
  {
    FEPIC_CHECK(cellid<=HalfT::cell_id_limit, "cell id limit exceeded", std::out_of_range);
    FEPIC_CHECK((pos>=-1)||(pos<=HalfT::position_limit), "position limit exceeded", std::out_of_range);
    FEPIC_CHECK(anchor<=HalfT::anchor_limit, "anchor limit exceeded", std::out_of_range);
    THIS->_incid_cell = cellid;
    THIS->_position = pos+1;
    THIS->_anchor = anchor;
  }

  /** Retorna um vetor com os índices dos vértices que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  vectorui getVertices(MeshT const& mesh) const
  {
    return mesh.getCell(this->getIncidCell()) // CELL->
                
           ->getBorderVertices(this->getPosition(), mesh);
  }


  /** Retorna um vetor com os índices dos nós que esta HalfFace contém.
  *  @param mesh a malha na qual a HalfFace está contida.
  */
  vectorui getNodes(MeshT const& mesh) const
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
  bool hasTheseVertices(vectorui const& v, MeshT const& mesh) const
  {
    FEPIC_CHECK(v.size()==CellT::n_vertices_per_border, "", std::invalid_argument);
    const auto cell = mesh.getCell(this->getIncidCell());
    
    vectorui vtx (cell->getBorderVertices(this->getPosition()));
    
    return arrayIsCyclicallyEqual(v, vtx);
  }
  
#undef THIS
#undef CONST_THIS

};



#endif

