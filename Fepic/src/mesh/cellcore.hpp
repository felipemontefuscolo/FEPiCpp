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

#ifndef FEPIC_CELLCORE_HPP
#define FEPIC_CELLCORE_HPP

template<class _Traits>
class _CellCore
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<CellT*>(this)
  #define CONST_THIS static_cast<const CellT*>(this)
#endif  

public:
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::HalfT  HalfT;
  typedef typename _Traits::HalfLT HalfLT;
  typedef typename _Traits::MeshT  MeshT;

protected:
  _CellCore() {};
  _CellCore(_CellCore const&) {};
  
public:
  
  int getOrder() const
  {
    return static_cast<int>(CONST_THIS->_order);
  }

  int getNumNodes() const
  {
    return CONST_THIS->_nodes.size();
  }

  int getNumBorders() const
  {
    return CellT::n_borders;
  }

  int getNumVertices() const
  {
    return CellT::n_vertices;
  }

  uint getNodeIdx(int ith) const
  {   
#ifdef FEPIC_DEBUG_ON
    return CONST_THIS->_nodes.at(ith);
#else
    return CONST_THIS->_nodes[ith];
#endif
  }

  vectorui getBorderVertices(int ith) const
  {
    //static matrixi faces_vtx(ElementProperties<CellT, _Traits>::get_faces_vtx());
    uint vsize = CellT::borders_local_vertices[ith].size();
    vectorui vtx(vsize);
    
    for (uint i = 0; i < vsize; ++i)
      vtx[i] = CONST_THIS->_nodes[CellT::borders_local_vertices[ith][i]];
      
    return vtx;
  }

  vectorui getBorderNodes(int ith) const
  {
    matrixi const& borders_local_nodes = CellT::getBordersLocalNodes(CONST_THIS->getOrder());
    uint tam(borders_local_nodes[ith].size());
    vectorui nodes(tam);
    
    for (uint i = 0; i < tam; ++i)
      nodes[i] =  CONST_THIS->_nodes[ borders_local_nodes[ith][i] ];
      
    return nodes;
  }
  
  void setNode(int ith, uint nodeid)
  {
#ifdef FEPIC_DEBUG_ON
    THIS->_nodes.at(ith) = nodeid;
#else
    THIS->_nodes[ith] = nodeid;
#endif
  }  

  HalfT* getHalf(int ith)
  {
    FEPIC_ASSERT(ith<CellT::n_borders, "out of range");
    return &THIS->_halfs[ith];
  }

  const HalfT* getHalf(int ith) const
  {
    FEPIC_ASSERT(ith<CellT::n_borders, "out of range");
    return &CONST_THIS->_halfs[ith];
  }

  void broadcastHalf2Nodes(MeshT & mesh) const
  {
    matrixi const& borders_local_nodes(CellT::getBordersLocalNodes(CONST_THIS->getOrder()));
    for (int f = 0; f < CellT::n_borders; ++f) // loop  nas faces
      for (uint i = 0, tam=borders_local_nodes[0].size(); i < tam; ++i)
        mesh.getNode(CONST_THIS->_nodes[borders_local_nodes[f][i]])->setHalf(CONST_THIS->_halfs[f]);
  }

  void broadcastTag2Nodes(MeshT & mesh, bool force=false) const
  {
    if (force)
      for (int i = 0; i < this->node.size(); ++i)
        mesh.getNode(THIS->_nodes[i])->setTag(this->getTag());
    else
      for (int i = 0; i < this->node.size(); ++i)
        if (mesh.getNode(THIS->_nodes[i])->getTag() == 0)
          mesh.getNode(THIS->_nodes[i])->setTag(this->getTag());
  }

#undef THIS
#undef CONST_THIS  
  
};



#endif
