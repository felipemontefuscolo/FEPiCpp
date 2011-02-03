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

#ifndef FEPIC_POLY3D_HPP
#define FEPIC_POLY3D_HPP



template<class Traits>
class _Poly3d : public _Labelable
{
public:
  static const int Dim = 3;
  
  /* Se Polygon está sendo instanciado, então CellT é com certeza uma FaceT*/
  typedef typename Traits::CellT  CellT;
  typedef typename ElementProperties<CellT, Traits>::FaceT FaceT;
  typedef typename Traits::MeshT      MeshT;
  
  static const int n_borders  = ElementProperties<CellT, Traits>::n_borders;
  static const int n_vertices = ElementProperties<CellT, Traits>::n_vertices;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
    
  _Poly3d() : _Labelable()
  {
    _node.resize(n_vertices);
  }
    
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  _Poly3d(Fepic::vectorui const& nodes, int label=0) : _Labelable(label), _node(nodes)
  {
#ifdef FEPIC_DEBUG_ON
    if (_node.size() < static_cast<uint>(n_vertices))
    {
      std::cout << "error:_Poly3d constructor: número insuficiente de nós.\n";
      throw;
    }
#endif
  }

  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  _Poly3d(Fepic::vectorui&& nodes, int label=0) : _Labelable(label), _node(std::move(nodes))
  {
#ifdef FEPIC_DEBUG_ON
    if (_node.size() < n_vertices)
    {
      std::cout << "error:_Poly3d constructor: número insuficiente de nós.\n";
      throw;
    }
#endif
  }

  int getNumNodes() const
  {
    return _node.size();
  }

  /** Retorna o índice do i-ésimo nó do polígono.
  */ 
  uint getNodeIdx(int ith) const
  {   
#ifdef FEPIC_DEBUG_ON
    return _node.at(ith);
#else
    return _node[ith];
#endif
  }

  /** Defini o i-ésimo nó do polígono como nodeid
  */ 
  void setNode(int ith, uint nodeid)
  {
#ifdef FEPIC_DEBUG_ON
    _node.at(ith) = nodeid;
#else
    _node[ith] = nodeid;
#endif
  }

  /** Retorna a i-ésima half-edge.
  */ 
  HalfFace<Traits>* getHalf(int ith)
  {
#ifdef FEPIC_DEBUG_ON
    if (ith >= n_borders)
    {
      std::cout << "getHalf: out of range" << std::endl;
      throw;
    }
#endif
    return &_hf[ith];
  }
  
  /** Retorna se os vértices passados formam uma face nesse poliedro.
  * @param[in] nodes um vetor com os vértices. O tamanho do vetor deve ser exatamente o número
  *            de lados (n_borders) do elemento degenerado deste poliedro.
  * @param[out] O índice local da face.
  * @param[out] a âncora.
  * @return true se forma uma face, false caso contrário.
  * @warning A ordem da numeração dos nós é importante!
  * @note definição para ancora: dada uma face n0,n1,...,nB, se ela é a face F
  * e ancora A de um poliedro, então, o nó B-A é o primeiro nó da face F do poliedro.
  */ 
  bool isAnFace(Fepic::vectorui const& vertices, int &face, int &anchor) const
  {
    static const Fepic::matrixi faces_vtx(ElementProperties<CellT, Traits>::get_faces_vtx());
  
    // CHECK
    if (int(vertices.size()) != FaceT::n_borders)
    {
      std::cout << "isAnFace error: wrong parameter!" << std::endl;
      std::cout << "vertices.size() = " << vertices.size() << std::endl;
    }
    
    Fepic::vectorui temp(FaceT::n_borders);
    Fepic::vectorui this_vtx(FaceT::n_borders);
    
    for (int f = 0; f < n_borders; ++f) // para cada face
    {
      for (int n = 0; n < FaceT::n_borders; n++)
      {
        this_vtx[n] = _node[ faces_vtx[f][n] ];
      }                       
      for (int a = 0; a < FaceT::n_borders; ++a) // para cada ancora
      {
        rotate_copy(vertices.begin(), vertices.begin()+a, vertices.end(), temp.begin()); // temp = rot(vertices + a);
        
        if(temp == this_vtx)
        {
          face = f;
          anchor = a;
          return true;
        }
              
      }
    }
    
    return false;
  
  }
  
  /** O mesmo que a função isAnFace, mas aqui a ordem dos vértices podem ser anti-paralelos a face.
  */ 
  bool isAnParallelFace(Fepic::vectorui vertices, int &face, int &anchor)
  {
    if (this->isAnFace(vertices, face, anchor))
      return true;
            
    reverse(vertices.begin(), vertices.end());
    return (this->isAnFace(vertices, face, anchor));
  }
  
  void printSelfState(std::ostream &o) const
  {
    o << getNodeIdx(0);
    for (uint i = 1, tam=this->_node.size(); i < tam; i++)
      o << " " << getNodeIdx(i);
  }       
  
  /** Retorna um vetor com os nós da i-ésima face.
  * @param ith a i-ésima face.
  */ 
  Fepic::vectorui getBorderNodes(int ith, MeshT const& mesh) const
  {
    Fepic::vectorui nodes;
    
    for (uint i = 0, tam=mesh.border_local_nodes[ith].size(); i < tam; i++)
      nodes.push_back( _node[ mesh.border_local_nodes[ith][i] ] );
      
    return nodes;
  }
  
  /** Retorna um vetor com os vértices da i-ésima face
  */ 
  Fepic::vectorui getBorderVertices(int ith) const
  {
    static Fepic::matrixi faces_vtx(ElementProperties<CellT, Traits>::get_faces_vtx());
    uint vsize = faces_vtx[ith].size();
    Fepic::vectorui vtx(vsize);
    
    for (uint i = 0; i < vsize; ++i)
      vtx[i] = _node[faces_vtx[ith][i]];
      
    return vtx;
  }
  
  /** Atribui a cada nó deste tetraedro sua respectiva half-face.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
    for (int f = 0; f < n_borders; ++f) // loop  nas faces
      for (uint i = 0, tam=mesh.border_local_nodes[0].size(); i < tam; ++i)
        mesh.getNode(_node[mesh.border_local_nodes[f][i]])->setHalf(_hf[f]);
  }
  
  /** Atribui o o label do tetraedro em seus nós.
  * @param force Caso true, a herança do label será feita sem restrição. Caso false,
  * a herança será feita apenas se cada nó não tiver label (i.e., label = 0).
  * @warning CONFERIR!!
  */ 
  void propagateTag(MeshT & mesh, bool force=false) const
  {
    if (force)
      for (int i = 0; i < this->node.size(); i++)
        mesh.getNode(_node[i])->setTag(this->getTag());
    else
      for (int i = 0; i < this->node.size(); i++)
        if (mesh.getNode(_node[i])->getTag() == 0)
          mesh.getNode(_node[i])->setTag(this->getTag());
  }
  
  
protected:
        Fepic::vectorui                 _node;
        HalfFace<Traits>       _hf[n_borders];    
};







#endif
