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

#ifndef FEPIC_HALFEDGE_HPP
#define FEPIC_HALFEDGE_HPP


template<class _Traits>
class HalfEdge : public _HalfCore<_Traits>, public _Labelable
{
public:
  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::CellT CellT;

  enum {cell_id_limit=134217727};
  enum {position_limit=6};
  enum {anchor_limit=3};
  
  friend class _HalfCore<_Traits>;

  /** 0 <= incid_cell <= HalfEdge::cell_id_limit (134217727)
   * -1 <= position   <= HalfEdge::position_limit (6)
   *  0 <= anchor     <= HalfEdge::anchor_limit (3)
   */ 
  template<class... LabelArgs>
  HalfEdge(int incid_cell, int position, int anchor, LabelArgs... args) :
                                                         _Labelable(args...),
                                                         _incid_cell(incid_cell),
                                                         _position(position),
                                                         _anchor(anchor)
  {
    FEPIC_CHECK((incid_cell<=cell_id_limit)&&
                 (position<=position_limit && position>-2)&&
                 (anchor<=anchor_limit), "", std::out_of_range);
  }
  
  HalfEdge(HalfEdge const&) = default;
  HalfEdge() : _Labelable(), _incid_cell(0), _position(0), _anchor(0) {}
  ~HalfEdge() = default;

  /** Imprime a composição do iD desta HalfEdge.
  *  @param o o stream onde se vai escrever, e.g., std::cout.
  */ 
  void printSelf(std::ostream& o)
  {
    o << this->getIncidCell() << " " << this->getPosition();
  }
  
  /** Retorna o comprimento da HalfEdge, i.e., da aresta que
  *  ela representa.
  */ 
  double getLenght(MeshT& mesh) const
  {
    vectori nodes(this->getNodes(mesh));
    double sum=0;
    for (int i = 0; i < nodes.size()-1; i++)
      sum += mesh.getNode(nodes[i])->getDistance(*mesh.getNode(nodes[i+1]));
    
    return sum;
  }
  
  /** @return uma string contendo "Half-Edge"
  */ 
  static std::string getName()
  {
    return std::string("Half-Edge");
  }

  
protected:
  int8_t _position;
  int8_t _anchor;
  int    _incid_cell;
};      






#endif
