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

#ifndef FEPIC_HALFFACELAB_HPP
#define FEPIC_HALFFACELAB_HPP




/** @class HalfFaceLab
 * 
 * Trata-se da classe HalfFace com herança de classe _Labelable.
 */ 
template<class _Traits>
class HalfFaceLab : public HalfFace<_Traits>, public _Labelable
{
public:

  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::CellT CellT;
  
  
  /** Construtor.
  *  @param incid_cell o iD da célula incidente.
  *  @param ith a posição da aresta nesta célula.
  *  @param anchor o índice âncora.
  *  @param tag o rótulo.
  */
  template<class... LabelArgs>
  HalfFaceLab(uint incid_cell, int ith, int anchor, LabelArgs... args) : HalfFace<_Traits>(incid_cell, ith, anchor),
                                                                         _Labelable(args...) {}
  HalfFaceLab(HalfFaceLab const&) = default;
  HalfFaceLab() = default;
  ~HalfFaceLab() = default;
  /** Faz com que cada nó desta HalfFaceLab aponte para ela.
  *  @param mesh a malha na qual a HalfFaceLab está contida.
  */ 
  void broadcastHalf2Nodes(MeshT& mesh) const
  {
    vectorui v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition()));
    
    for (uint i = 0; i < v.size(); i++)
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
  

};






#endif
