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
class HalfEdge : public _HalfCore<_Traits>
{
public:
  typedef typename _Traits::MeshT MeshT;
  typedef typename _Traits::CellT CellT;

  enum {cell_id_limit=134217727};
  enum {position_limit=6};
  enum {anchor_limit=3};
  
  friend class _HalfCore<_Traits>;

  HalfEdge(uint incid_cell, int position, int=0) : _incid_cell(incid_cell), _position(position)
  {
  
  }
  
  /** Construtor.
  */ 
  HalfEdge() : _incid_cell(0), _position(0)
  {
  }


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
    vectorui nodes(this->getNodes(mesh));
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
  
  /** Destrutor.
  */ 
  ~HalfEdge() {}
      
protected:
  uint _incid_cell : 27;
  uint _position   : 3;
  uint _anchor     : 2;

};      







#endif