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

#ifndef FEPIC_HALFFACE_HPP
#define FEPIC_HALFFACE_HPP


template<class _Traits>
class HalfFace : public _HalfCore<_Traits>, public _Labelable
{
public:
  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::CellT CellT;

  enum {cell_id_limit=134217727};
  enum {position_limit=6};
  enum {anchor_limit=3};
  
  friend class _HalfCore<_Traits>;

  template<class... LabelArgs>
  HalfFace(int incid_cell, int position, int anchor, LabelArgs... args) :
                                                         _Labelable(args...),
                                                         _incid_cell(incid_cell),
                                                         _position(position),
                                                         _anchor(anchor) {}
  
  HalfFace(HalfFace const&) = default;
  HalfFace() : _Labelable(), _incid_cell(-1), _position(-1), _anchor(-1) {}
  ~HalfFace() = default;
  
  /** Imprime em um stream a composição do iD, i.e, imprime \n
  *  getIncidCell() << " " << getPosition() << " " << getAnchor()
  *  @param o o stream onde se vai escrever.
  */ 
  void printSelf(std::ostream& o) const
  {
    o << this->getIncidCell() << " " << this->getPosition() << " " << this->getAnchor();
  }
  
  /** @return uma string com o nome Half-Face.
  */ 
  static std::string getName()
  {
    return std::string("Half-Face");
  }
  
  /** Faz com que cada nó desta HalfFaceLab aponte para ela.
  *  @param mesh a malha na qual a HalfFaceLab está contida.
  */ 
  void broadcastHalf2Nodes(MeshT& mesh) const
  {
    vectori v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition()));
    
    for (int i = 0; i < v.size(); i++)
    {
      mesh.getNode(v[i])->setHalf(*this);
    }
  }
  
  /** Atribui o rótulo desta _MetaHalfLabOf a seus nós.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem tag=0;
  */ 
  void broadcastTag2Nodes(MeshT& mesh, bool force=false) const
  {
    vectori v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition(),mesh ));
    
    if (force)
      for (int i = 0; i < v.size(); ++i)
        mesh.getNode(v[i])->setTag(this->getTag());
    else
      for (int i = 0; i < v.size(); i++)
        if (mesh.getNode(v[i])->getTag() == 0)
          mesh.getNode(v[i])->setTag(this->getTag());
  } 

protected:
  int8_t  _position;
  int8_t  _anchor;
  int     _incid_cell;
};      







#endif
