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

#ifndef FEPIC_HALFEDGELAB_HPP
#define FEPIC_HALFEDGELAB_HPP


/** @class HalfEdgeLab
 * 
 * Trata-se da classe HalfEdge com herança de classe _Labelable.
 */ 
template<class Traits>
class HalfEdgeLab : public HalfEdge<Traits>, public _Labelable
{
public:

  typedef typename Traits::MeshT MeshT;    
  
  /** Construtor.
  *  @param cellid O iD da célula incidente.
  *  @param ith A posição da aresta nesta célula.
  *  @param label o rótulo.
  */ 
  HalfEdgeLab(uint cellid, int ith, int label=0) : HalfEdge<Traits>(cellid, ith), _Labelable(label) {}
  HalfEdgeLab() : HalfEdge<Traits>(), _Labelable() {};
  
  /** Faz com que cada nó desta HalfEdgeLab aponte para ela.
  *  @param mesh a malha na qual a HalfEdgeLab está contida.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
    Fepic::vectorui v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition(), mesh ));
    
    for (uint i = 0; i < v.size(); i++)
      mesh.getNode(v[i])->setHalf(*this);
    
  }
  
  /** Atribui o rótulo desta HalfEdgeLab a seus nós.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem label=0;
  */ 
  void propagateTag(MeshT & mesh, bool force=false) const
  {
    Fepic::vectori v (mesh.getCell( this->getIncidCell() )->getBorderNodes( this->getPosition(), mesh));
    
    if (force)
      for (int i = 0; i < v.size(); ++i)
        mesh.getNode(v[i])->setTag(this->getTag());
    else
      for (int i = 0; i < v.size(); i++)
        if (mesh.getNode(v[i])->getTag() == 0)
          mesh.getNode(v[i])->setTag(this->getTag());
  }                                                                                                                        
  
  /** Destrutor.
  */ 
  ~HalfEdgeLab() {}
};







#endif
