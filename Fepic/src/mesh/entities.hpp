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

#ifndef ENTITIES_HPP    
#define ENTITIES_HPP

#ifndef ENTTS_DEBUG
  #define ENTTS_DEBUG
#endif

#include "Fepic/src/mesh/msh_tags.hpp"
#include "Fepic/src/util/misc.hpp"
#include "Fepic/src/custom_eigen/custom_eigen.hpp"
#include "Fepic/src/shapefunctions/parametric_pts.hpp"
#include <iostream>
#include <iomanip>
#include <deque>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <type_traits>


/* pre declarações */

template<class Traits> class iMesh;

class iLabel;

// volume
template<class Traits> class iPolyhedron;
template<class Traits> class iTetrahedron;
template<class Traits> class iHexahedron;

// surface
template<class Traits> class iPolygon;
template<class Traits> class iTriangle;
template<class Traits> class iQuadrangle;

// line
template<class Traits> class iEdge;

//iPoint
template<class Traits> class iPoint;

// Halfs
template<class Traits> class iHalfEdge;
template<class Traits> class iHalfFace;

// MHalfs
template<class Traits> class iMarkedHalfEdge;
template<class Traits> class iMarkedHalfFace;
/*--------------------------------------------------------------------*/

template<int dim>
class Simplex {};

template<int dim>
class Hypercube {};


template<class ElmType, class Traits> class ElementProperties;

template<class Traits>
class ElementProperties<Simplex<1>, Traits> {
public:
  typedef iEdge<Traits> Type;
};

template<class Traits>
class ElementProperties<Simplex<2>, Traits> {
public:
  typedef iTriangle<Traits> Type;
};

template<class Traits>
class ElementProperties<Simplex<3>, Traits> {
public:
  typedef iTetrahedron<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<1>, Traits> {
public:
  typedef iEdge<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<2>, Traits> {
public:
  typedef iQuadrangle<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<3>, Traits> {
public:
  typedef iHexahedron<Traits> Type;
};

template<class Traits>
class ElementProperties<iEdge<Traits>, Traits> {
public:
  static const int N_borders = 2;
  static const int N_vertices = 2;
};
template<class Traits>
class ElementProperties<iTriangle<Traits>, Traits> {
public:
  static const int N_borders = 3;
  static const int N_vertices = 3;
};
template<class Traits>
class ElementProperties<iTetrahedron<Traits>, Traits> {
public: static const int N_borders = 4;
  static const int N_vertices = 4;
  static Fepic::matrixi get_faces_vtx()
  {
    static const Fepic::matrixi temp = { {1,0,2}, {0,1,3}, {3,2,0}, {2,3,1} };
    return temp;
  }
  static Fepic::matrixi get_edges_vtx()
  {
    static const Fepic::matrixi temp = { {0,1}, {1,2}, {2,0}, {3,0}, {3,2}, {3,1} };
    return temp;
  }
  
  typedef iTriangle<Traits> FaceT;
};


class UndefElement {
public:
  static int getMSHTag()
  {
    return MSH_UNDEFINED_ELEM;
  }
};

/* Define type of volume: Cell::dim < 3 ? UndefVol : Cell::Volume */
template<int CellDim, class Cell>
class VolumeDef; 

template<class Cell>
class VolumeDef<1, Cell> {
public:
  typedef UndefElement VolumeT;
}; 

template<class Cell>
class VolumeDef<2, Cell> {
public:
  typedef UndefElement VolumeT;
};

template<class Cell>
class VolumeDef<3, Cell> {
public:
  typedef typename Cell::VolumeT VolumeT;
};

/* Define the type of face */
template<int Dim, class Cell>
class FaceDef;

template<class Cell>
class FaceDef<1, Cell> {
public:
  typedef UndefElement FaceT;
};

template<class Cell>
class FaceDef<2, Cell> {
public:
  typedef typename Cell::FaceT FaceT;
};

template<class Cell>
class FaceDef<3, Cell> {
public:
  typedef typename Cell::FaceT FaceT;
};

//------------------------------------

template<int Dim, class Traits>
class HalfDef;

template<class Traits>
class HalfDef<1, Traits> {
public: 
  typedef iHalfEdge<Traits> HalfT;
};

template<class Traits>
class HalfDef<2, Traits> {
public: 
  typedef iHalfEdge<Traits> HalfT;
};

template<class Traits>
class HalfDef<3, Traits> {
public: 
  typedef iHalfFace<Traits> HalfT;
};


//------------------------------------

template<int Dim, class Traits>
class MHalfDef;

template<class Traits>
class MHalfDef<1, Traits> {
public: 
  typedef iMarkedHalfEdge<Traits> MHalfT;
};

template<class Traits>
class MHalfDef<2, Traits> {
public: 
  typedef iMarkedHalfEdge<Traits> MHalfT;
};

template<class Traits>
class MHalfDef<3, Traits> {
public: 
  typedef iMarkedHalfFace<Traits> MHalfT;
};


/** @class DefaultTraits
 * 
 * Default definitions of a traits.
 * 
 * Users can do their own traits.
 * 
 * In custom traits must be defined:
 * - CellT      := type of grid cell
 * - EdgeT
 * - PointT
 * - MeshT
 * - spacedim   := dimension of the space
 * 
 * 
 */
template<int _spacedim, class CellType = Simplex<_spacedim> >
class DefaultTraits {
public:

  DefaultTraits(DefaultTraits const&) = delete; // dont copy me
  ~DefaultTraits() = delete;
  
  typedef DefaultTraits Traits;
  
  typedef typename ElementProperties<CellType, Traits>::Type CellT;
  
  typedef iPoint<Traits>  PointT;
  
  typedef iMesh<Traits>   MeshT;
  
  static const int spacedim = _spacedim;
};


// ---------------------------------------------


//======================================================================
//==========================  iLabel =============================
//======================================================================

/**
 * @class iLabel
 * 
 *  Interface para armazenar um rótulo e 3 flags. Esses atributos são
 *  compactados em um único char. Os 3 primeiros bits do char são para
 *  armazenar os flags:
 * 
 *  - char[0] := "isdead" ?;  
 *  - char[1] := "must be killed" ?;
 *  - char[2] := "other flags"
 * 
 *  O restante é para armazenar o rótulo (inteiro de 0 a 31 se o char for de 8 bits).
 */ 
class iLabel
{
public:
  /** @param Label O rótulo (inteiro de 0 a 31).
  *  @param Isdead true se está morto ou false caso contrário.
  *  @param MustBK true se deve ser morto ou false caso contrário.
  *  @param otherf Um flag qualquer. Ainda não lhe foi atribuído nenhuma função.
  */ 
  iLabel(int Label, bool Isdead=false, bool MustBK=false, bool otherf=false)
  {
    _label = (unsigned char)(Label << OffSet);
    Isdead ? set_bit(_label,0) : unset_bit(_label,0);
    MustBK ? set_bit(_label,1) : unset_bit(_label,1);
    otherf ? set_bit(_label,2) : unset_bit(_label,2);
  }
  
  iLabel()
  {
    _label = 0;
  };
  
  /** @return O rótulo.
  */ 
  int getLabel() const
  {
    return int(_label >> OffSet);
  }
  
  /** Seta o rótulo.
  *  @param c O rótulo.
  */ 
  void setLabel(unsigned char const& c)
  {
    _label &= FlagsMask;
    _label |= (c << OffSet);
  }
  
  /** Limpa os flags.
  */ 
  void resetFlags()
  {
    _label &= (~FlagsMask);
  }
  
  /** @return true se está morto ou false caso contrário.
  */ 
  bool isDead() const
  {
    return get_bit(_label,0);
  }
  
  /** @return true se deve ser morto ou false caso contrário.
  */ 
  bool mustBeKilled() const
  {
    return get_bit(_label,1);
  }
  
  /** Seta o flag "isdead".
  */ 
  void killMe()
  {
    set_bit(_label,0);
  }
  
  /** Limpa o flag "isdead".
  */ 
  void reviveMe()
  {
    unset_bit(_label,0);
  }
  
  /** Seta o flag "must be killed".
  */ 
  void killMeLater()
  {
    set_bit(_label,1);
  }
  
  /** Limpa o flag "must be killed".
  */ 
  void dontKillMe() {
    unset_bit(_label,1);
  }
  
  ~iLabel() {}
  
  
  
protected:
  unsigned char _label;
  static const unsigned char FlagsMask = 7;
  static const int  OffSet = 3;   
};



//======================================================================
//============================ iHalfEdge ===============================
//======================================================================

/** @class iHalfEdge
 *
 * Armazena o iD da célula incidente e o seu índice local em relação a esta célula.
 * 
 * O iD da iHalfEdge é feita da seguinte composição
 * 
 * <f, ith> \n\n
 * onde
 *        - f:= iD da célula (polígono)            [ 28bits ]   
 *        - ith := i-ésima aresta desta célula [  4bits ]
 * 
 * Essas informações são armazenadas em uma única variável (uint) [32bits]
 * usando manipulação de bits. Os valores de f e ith estão, respectivamente, na faixa:
 * 
 * - 0  <= f <= 268 435 456
 * - -1 <= ith <= 14
 * 
 * @warning o tipo unsigned int deve ter 32 bits.
 * 
 */ 
template<class Traits>
class iHalfEdge {

public:
  typedef typename Traits::MeshT MeshT;
  typedef typename Traits::CellT CellT;

  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;

  /** Construtor.
  *      @param cellid O iD da célula incidente.
  *  @param ith A posição da aresta nesta célula.
  */ 
  iHalfEdge(uint cellid, int ith, int=0)
  {
    _compositeID = (cellid << OffSet);
    _compositeID |= (uint)(ith+1);
    /* atenção na soma de 1 : é um truque para mais tarde poder retornar -1 */
  }
  
  /** Construtor.
  */ 
  iHalfEdge()
  {
    this->_compositeID = 0;
    this->setIDPosition(-1);
  }
  
  /** @return iD da célula incidente (triângulo, quadrângulo, ...)
  */     
  uint getIDCell() const
  {
    return _compositeID >> OffSet;
  }
  
  /** @return A posição desta iHalfEdge na célula
  */     
  int getIDPosition() const
  {
    return (_compositeID & IndexMask) - 1;
  }
  
  /** Seta o iD da célula incidente.
  *  @param cellid índice da célula.
  */ 
  void setIDCell(uint cellid)
  {
    _compositeID &= IndexMask;
    _compositeID |= (cellid << OffSet);
  }
  
  /** Seta a posição da HE na célula.
  *  @param ith A posição.
  */ 
  void setIDPosition(int ith)
  {
    _compositeID &= ~IndexMask;
    _compositeID |= (uint)(ith+1);
  }
  
  /** @param cellid O iD da célula incidente.
  *  @param ith A posição da aresta nesta célula.
  */     
  void setCompleteID(uint cellid, int ith)
  {
    _compositeID = (cellid << OffSet);
    _compositeID |= (uint)(ith+1);
  }
  
  /** Imprime a composição do iD desta iHalfEdge.
  *  @param o o stream onde se vai escrever, e.g., std::cout.
  */ 
  void printSelf(std::ostream& o) const
  {
    o << getIDCell() << " " << getIDPosition();
  }
  
  /** Retorna o comprimento da iHalfEdge, i.e., da aresta que
  *  ela representa.
  */ 
  double getLenght(MeshT& mesh) const
  {
    Fepic::vectorui nodes(std::move(this->getNodes(mesh)));
    double sum=0;
    for (int i = 0; i < nodes.size()-1; i++)
      sum += mesh.getNode(nodes[i])->getDistance(*mesh.getNode(nodes[i+1]));
    
    return sum;
  }
  
  /** Retorna um vetor com os índices dos vértices que esta iHalfEdge contém.
  *  @param mesh a malha na qual a iHalfEdge está contida.
  */ 
  Fepic::vectorui getVertices(MeshT& mesh)
  {
    CellT *cell = mesh.getCell(this->getIDCell());
    
    return cell->getBorderVertices(this->getIDPosition());
  }
  
      /** Retorna um vetor com os índices dos nós que esta iHalfEdge contém.
      *  @param mesh a malha na qual a iHalfEdge está contida.
      */ 
  Fepic::vectorui getNodes(MeshT& mesh) const
  {
    CellT *cell = mesh.getCell(this->getIDCell());
    
    return cell->getBorderNodes(this->getIDPosition(),mesh);
  }
      
  /** Verifica se esta iHalfEdge contém exatamente os vértices (índices) passados
  *  em v.
  *  @param v vetor com os índices dos vértices.
  *  @param mesh a malha na qual esta iHalfEdge está contida.
  *  @return true se os nós formam a iHalfEdge e false caso contrário.
  *  @note Os nós podem estar em qualquer orientação cíclica da original.
  */
  bool hasTheseVertices(Fepic::vectorui const& v, MeshT& mesh) const
  {
    Fepic::vectorui vtx(2);
    
    CellT *cell = mesh.getCell(this->getIDCell());
    
    vtx = cell->getBorderVertices(this->getIDPosition()) ;
    
    return arrayIsCyclicallyEqual(v, vtx);
  }
  
  /** @return uma string contendo "Half-Edge"
  */ 
  static std::string getName()
  {
    return std::string("Half-Edge");
  }
  
  /** Destrutor.
  */ 
  ~iHalfEdge() {}
      
protected:
  uint _compositeID;
  static const uint IndexMask = 15; // [00000000000000000000000000001111]
  static const int OffSet = 4;
};      

//======================================================================
//============================ iHalfFace ===============================
//======================================================================

/** @class iHalfFace
 *
 * Armazena o iD da célula incidente, o seu índice local em relação a esta célula e
 * a sua âncora (ver definição no final da descrição).
 * 
 * O iD da iHalfEdge é feita da seguinte composição
 * 
 * <v, ith, anchor> \n
 * 
 * onde
 * 
 * - v := iD da célula (poliedro)                    [ 27bits ]
 * - ith := i-ésima face desse volume        [  3bits ]
 * - anchor := k-ésimo nó da face do volume  [  2bits ]
 * 
 * Essas informações são armazenadas em uma única variável (uint) [32bits]
 * usando manipulação de bits. Os limites de v, ith e anchor são, respectivamente:
 * 
 * - 0  <= v <= 134 217 727
 * - -1 <= ith <= 6
 * - 0 <= anchor <= 4
 * 
 * <b> definição para âncora </b>: dada uma iHalfFace que contenha os vértices n0,n1,...,nB,
 * se ela é a face F e âncora A de um poliedro, então, o nó B-A é o primeiro nó da face F do poliedro. 
 */ 
template<class Traits>
class iHalfFace
{
public:
  typedef typename Traits::MeshT MeshT;
  typedef typename Traits::CellT CellT;
  
  /** @param cellid O iD da célula incidente.
  *  @param ith A posição da face nesta célula.
  *  @param anchor O índice âncora.
  */ 
  iHalfFace(uint cellid, int ith, int anchor)
  {
    /* <v, ith, anchor>
      27b   3b    2b    */
    
    _compositeID = (cellid << OffSet1);
    _compositeID |= ((uint)(ith+1) << OffSet0);
    _compositeID |= (uint)(anchor);
  }
  
  /** Construtor.
  */ 
  iHalfFace()
  {
    this->_compositeID = 0;
    this->setIDPosition(-1);
  }
  
  /** @return iD da célula incidente
  */ 
  uint getIDCell() const
  {
    return _compositeID >> OffSet1;
  }
  
  /** @return A posição desta iHalfFace na célula
  */ 
  int getIDPosition() const
  {
    return ((_compositeID & IndexMask) >> OffSet0) - 1;
  }
  
  /** @return O índice âncora
  */ 
  int getIDAnchor() const
  {
    return _compositeID & AnchorMask;
  }
  
  /** Seta o iD da célula incidente.
  *  @param cellid O índice da célula.
  */ 
  void setIDCell(uint cellid)
  {
    _compositeID &= FlagsMask;
    _compositeID |= (cellid << OffSet1);
  }
  
  /** Seta o índice local da iHalfFace na célula incidente.
  *  @param ith o índice.
  */ 
  void setIDPosition(int ith)
  {
    _compositeID &= ~IndexMask;
    _compositeID |= (uint)(ith+1) << OffSet0;
  }
  
  /** Seta o índice âncora.
  *  @param anchor O índice âncora.
  */ 
  void setIDAnchor(int anchor)
  {
    _compositeID &= ~AnchorMask;
    _compositeID |= (uint)(anchor);
  }
  
  /** @param cellid O iD da célula incidente.
  *  @param ith O índice local desta iHalfFace na célula incidente.
  *  @param anchor O índice âncora.
  */ 
  void setCompleteID(uint v, int ith, int anchor)
  {
    /* <v, ith, anchor>
        R   3b    2b    */
    
    _compositeID = (v << OffSet1);
    _compositeID |= ((uint)(ith+1) << OffSet0);
    _compositeID |= (uint)(anchor);
  }
  
  /** Imprime em um stream a composição do iD, i.e, imprime \n
  *  getIDCell() << " " << getIDPosition() << " " << getIDAnchor()
  *  @param o o stream onde se vai escrever.
  */ 
  void printSelf(std::ostream& o) const
  {
    o << getIDCell() << " " << getIDPosition() << " " << getIDAnchor();
  }
  
  /** @return uma string com o nome Half-Face.
  */ 
  static std::string getName()
  {
    return std::string("Half-Face");
  }
  
  /** Retorna um vetor com os índices dos vértices que esta iHalfFace contém.
  *  @param mesh a malha na qual a iHalfFace está contida.
  */
  Fepic::vectorui getVertices(MeshT& mesh)
  {
    CellT *cell = mesh.getCell(this->getIDCell());
    
    return cell->getBorderVertices(this->getIDPosition(), mesh);
  }

  /** Retorna um vetor com os índices dos nós que esta iHalfFace contém.
  *  @param mesh a malha na qual a iHalfFace está contida.
  */
  Fepic::vectorui getNodes(MeshT& mesh)
  {
      CellT *cell = mesh.getCell(this->getIDCell());
      
      return cell->getBorderNodes(this->getIDPosition(), mesh);
  }

  /** Verifica se esta iHalfFace contém exatamente os vértices (índices) passados
  *  em v.
  *  @param v vetor com os índices dos vértices.
  *  @param mesh a malha na qual esta iHalfFace está contida.
  *  @return true se os nós formam a iHalfFace e false caso contrário.
  *  @note Os nós podem estar em qualquer orientação cíclica da original.
  */
  bool hasTheseVertices(Fepic::vectorui const& v, MeshT& mesh)
  {
    Fepic::vectorui vtx(mesh.faces_local_nodes[0].size());
    
    CellT *cell = mesh.getCell(this->getIDCell());
    
    vtx = cell->getBorderVertices(this->getIDPosition()) ;
    
    return arrayIsCyclicallyEqual(v, vtx);
  }
        
  /** Destrutor.
  */ 
  ~iHalfFace() {}
protected:
  uint _compositeID;
  static const uint AnchorMask = 3; //  [00000000000000000000000000000011]
  static const uint IndexMask  = 28; // [00000000000000000000000000011100]
  static const uint FlagsMask  = 31; // [00000000000000000000000000011111]
  static const uint IdMask   = (~31);// [11111111111111111111111111100000]
  static const int OffSet0  = 2;
  static const int OffSet1  = 5;
};      


//======================================================================

/** @class iMarkedHalfEdge
 * 
 * Trata-se da classe iHalfEdge com herança de classe iLabel.
 */ 
template<class Traits>
class iMarkedHalfEdge : public iHalfEdge<Traits>, public iLabel
{
public:

  typedef typename Traits::MeshT MeshT;    
  
  /** Construtor.
  *  @param cellid O iD da célula incidente.
  *  @param ith A posição da aresta nesta célula.
  *  @param label o rótulo.
  */ 
  iMarkedHalfEdge(uint cellid, int ith, int label=0) : iHalfEdge<Traits>(cellid, ith),
                                                                                                          iLabel(label) {}
  iMarkedHalfEdge() : iHalfEdge<Traits>(), iLabel() {};
  
  /** Faz com que cada nó desta iMarkedHalfEdge aponte para ela.
  *  @param mesh a malha na qual a iMarkedHalfEdge está contida.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
    Fepic::vectorui v = mesh.getCell( this->getIDCell() )->getBorderNodes( this->getIDPosition(), mesh );
    
    for (uint i = 0; i < v.size(); i++)
      mesh.getNode(v[i])->setHalf(*this);
    
  }
  
  /** Atribui o rótulo desta iMarkedHalfEdge a seus nós.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem label=0;
  */ 
  void propagateLabel(MeshT & mesh, bool force=false) const
  {
    Fepic::vectori v = mesh.getCell( this->getIDCell() )->getBorderNodes( this->getIDPosition(), mesh);
    
    if (force)
      for (int i = 0; i < v.size(); ++i)
        mesh.getNode(v[i])->setLabel(this->getLabel());
    else
      for (int i = 0; i < v.size(); i++)
        if (mesh.getNode(v[i])->getLabel() == 0)
          mesh.getNode(v[i])->setLabel(this->getLabel());
  }                                                                                                                        
  
  /** Destrutor.
  */ 
  ~iMarkedHalfEdge() {}
};

/** @class iMarkedHalfFace
 * 
 * Trata-se da classe iHalfFace com herança de classe iLabel.
 */ 
template<class Traits>
class iMarkedHalfFace : public iHalfFace<Traits>, public iLabel
{
public:

  typedef typename Traits::MeshT MeshT;
  typedef typename Traits::CellT CellT;
  
  
  /** Construtor.
  *  @param cellid o iD da célula incidente.
  *  @param ith a posição da aresta nesta célula.
  *  @param anchor o índice âncora.
  *  @param label o rótulo.
  */ 
  iMarkedHalfFace(uint cellid, int ith, int anchor, int label=0) : iHalfFace<Traits> (cellid, ith, anchor),
                                                                                                                                  iLabel(label) {}
  
  iMarkedHalfFace() : iHalfFace<Traits>(), iLabel() {}
  
  /** Faz com que cada nó desta iMarkedHalfFace aponte para ela.
  *  @param mesh a malha na qual a iMarkedHalfFace está contida.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
    Fepic::vectorui v = mesh.getCell( this->getIDCell() )->getBorderNodes( this->getIDPosition(),mesh );
    
    for (uint i = 0; i < v.size(); i++)
    {
      mesh.getNode(v[i])->setHalf(*this);
    }
  }
  
  /** Atribui o rótulo desta iMarkedHalfEdge a seus nós.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem label=0;
  */ 
  void propagateLabel(MeshT & mesh, bool force=false) const
  {
    Fepic::vectori v = mesh.getCell( this->getIDCell() )->getBorderNodes( this->getIDPosition(),mesh );
    
    if (force)
      for (int i = 0; i < v.size(); ++i)
        mesh.getNode(v[i])->setLabel(this->getLabel());
    else
      for (int i = 0; i < v.size(); i++)
        if (mesh.getNode(v[i])->getLabel() == 0)
          mesh.getNode(v[i])->setLabel(this->getLabel());
  }                                                                                                                        
  
  /** Destrutor.
  */ 
  ~iMarkedHalfFace() {}
};



//======================================================================
//============================== iPoint ================================
//======================================================================

/** 
 *  Os objetos desta classe representam pontos no espaço.
 */ 
template<class Traits>
class iPoint : public iLabel {
public:
  static const int Dim = 0;
  
  typedef typename Traits::CellT                         CellT;
  typedef typename VolumeDef<CellT::Dim, CellT>::VolumeT VolumeT;      // CellT::dim < 3 ? UndefVol : CellT::Volume
  typedef typename FaceDef<CellT::Dim, CellT>::FaceT     FaceT;        // análogo
  typedef typename HalfDef<CellT::Dim, Traits>::HalfT    HalfT;
      
  typedef iEdge<Traits>          EdgeT;
  typedef typename Traits::MeshT MeshT;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
    
  /** Construtor.
  *  @param coord um vetor com Dim elementos que armazena a coordenada.
  *  @param label seu rótulo.
  */ 
  template<class T>
  iPoint(T const& coord, char label=0) : iLabel(label)
  {
    for (int i = 0; i < Traits::spacedim; ++i)
      _coord[i] = coord[i];
  }
  
  /** Construtor.
  */ 
  iPoint() : iLabel() {}
  // construtor de cópia não necessário, pois não há nenhum ponteiro.
  
  /** @return a dimensão do espaço.
  */  
  int getSpaceDim() const 
  {
    return Traits::spacedim;
  }
  
  /** Faz com que este ponto tenha a HalfT ha.
  */ 
  void setHalf(HalfT const& ha)
  {
    _half = ha;
  }
  
  /** Retorna um ponteiro para sua HalfT.
  */ 
  HalfT* getHalf()
  {
    return &_half;
  }
  
  /** Retorna um ponteiro para sua HalfT.
  */ 
  const HalfT* getHalf() const
  {
    return &_half;
  }
          
  /** Imprime as coordenadas de um ponto em um stream dado.
  *  @param o stream onde se vai imprimir.
  *  @param space Espaço entre a impressão de cada dimensão da coordenada.
  */ 
  void printSelfVTK(std::ostream &o, int space = 22) const
  {
    switch (Traits::spacedim) {
      case 1:
      {
        o << std::left << std::setw(space) << _coord[0] << " 0.0 0.0";
        break;
      }
      case 2:
      {
        o << std::left << std::setw(space) << _coord[0]
                       << std::setw(space) << _coord[1] << " 0.0";
        break;
      }
      case 3:
      {
        o << std::left << std::setw(space) << _coord[0]
                       << std::setw(space) << _coord[1]
                       << std::setw(space) << _coord[2];
        break;
      }
      default:
      {
        std::cout << "Error: invalid dimension\n";
        break;
      }
    }
  }
  
  /** Define a coordenada deste ponto.
  *  @param coord um vetor com a coordenada.
  */ 
  template<class Vec>
  void setCoord(Vec const& coord) 
  {
    for (int i = 0; i < Traits::spacedim; ++i)
      _coord[i] = coord[i];
  }
  
  /** Retorna a coordenada deste ponto em coord.
  *  @param[out] coord a coordenada.
  */ 
  template<class Vec>
  void getCoord(Vec & coord) const 
  {
    for (int i = 0; i < Traits::spacedim; ++i)
      coord[i] = _coord[i];
  }
  
  /** Retorna a i-ésima componente da coordenada
  */ 
  double getCoord(int i) const
  {
    return _coord[i];
  }
  
  /** Retorna a coordenada deste ponto.
  */
  VecT getCoord() const
  {
    return _coord;
  }
  
  /** Retorna a distância até a coordenada p.
  */     
  double getDistance(VecT const& p) const
  {
    return sqrt( (_coord-p).dot(_coord-p) );
  }
  
  /** Retorna a distância até o ponto p.
  */  
  double getDistance(iPoint const& p) const
  {
    return sqrt( (_coord-p._coord).dot(_coord-p._coord) );
  }
  
  
  /** Retorna o tag correspondente ao formato de arquivo .msh
  */ 
  static int getMSHTag()
  {
    return MSH_PNT;
  }
  
  /** NOT FOR USERS
  */ 
  static void setOrder()
  {
    // nada
  }
  
  /** Destrutor.
  */ 
  ~iPoint() {}
  
  static const int N_borders = 1;
protected:
  VecT _coord;
  HalfT _half;
};



//======================================================================
//============================== iEdge =================================
//======================================================================

/**
 *  Esta classe representa segmentos de linha orientados. As iEdge podem
 * ter 2 ou mais nós, um nó em cada extremo e o restante no seu interior
 * com espaçamento constante. Uma iEdge com N nós tem a seguinte forma: \n
 * 
 * *n0_____*n2______*n3____ ... ____*nN-1_____*n1 \n
 * 
 * Uma iEdge de N+1 nós é dita iEdge de ordem N, e sua orientação é definida
 * como a direção ao longo da iEdge de n0 até n1.
 * 
 */ 
template<class Traits>
class iEdge : public iLabel {
public:
  static const int Dim = 1;

  typedef Simplex<1>  ElmClass;
    
  typedef typename Traits::CellT                         CellT;
  typedef typename VolumeDef<CellT::Dim, CellT>::VolumeT VolumeT;
  typedef typename FaceDef<CellT::Dim, CellT>::FaceT     FaceT; 
    
  typedef typename Traits::PointT  PointT;       
  typedef typename Traits::MeshT   MeshT;        
  typedef typename Traits::PointT  BorderT;  
  typedef UndefElement             BndBorderT;

  /** Construtor.
  * @param nodes vetor com os nós que compõe a edge.
  * @note devem ser passados pelo menos dois nós.
  */ 
  iEdge(Fepic::vectorui const& nodes, int label=0) : iLabel(label), _node(nodes)
  {
#ifdef ENTTS_DEBUG
    if (_node.size() < 2)
    {
      std::cout << "erro: iEdge constructor: devem ser passados pelo menos dois nós\n";
      throw;
    }
#endif
  }
    
  /** Construtor.
  */ 
  iEdge() : iLabel(), _node({0,0})
  {
  }
  
  /** Retorna o número de nós de uma iEdge dada uma ordem order.
  */ 
  static int getNumNodesCell(int const order)
  {
    return order + 1;
  }
  
  /** Retorna o número de nós desta iEdge.
  */ 
  int getNumNodes() const
  {
    return _node.size();
  }
  
  /** Retorna o ith-ésimo nó desta iEdge.
  */ 
  uint getNodeIdx(int const ith) const
  {
    return _node[ith];
  }
  
  /** Verifica se esta iEdge está alinha com outra iEdge. O alinhamento é verificado
  * comparando-se os nós dos extremos das iEdges.
  *  @param e edge na qual se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (iEdge const& e) const
  {
    if ( (_node[0] == e._node[0]) && (_node[1] == e._node[1]) )
      return 1;
    else if ( (_node[0] == e._node[1]) && (_node[1] == e._node[0]) )
      return -1;
    else
      return 0;
  }
  
  /** Verifica se esta iEdge está alinha com dois nós dados. O alinhamento é verificado
  * comparando-se os nós dos extremos desta iEdge com os dois nós.
  *  @param nodes os vetor com os dois nós nos quais se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (Fepic::vectorui const& nodes) const
  {
    if ((_node[1] == nodes[1]) && (_node[0] == nodes[0]) )
      return 1;
    else if ( (_node[0] == nodes[1]) && (_node[1] == nodes[0]) )
      return -1;
    else
      return 0;
  }       
  
  /** Atribui o i-ésimo nó da aresta como o n-ésimo nó da malha.
  */ 
  void setNode(int const ith, uint const nth)
  {
    _node[ith] = nth;
  }
  
  /** NÃO IMPLEMENTADO
  *  Retorna o comprimento desta iEdge.
  */ 
  double lenght() const; // IMPLEMENTAR E LEMBRAR QUE PODER SER LINHA CURVA
  
  /** Imprime os nós da aresta. 
  *  @param o o stream aonde se vai escrever.
  *  @param order a ordem da aresta.
  */ 
  void printSelfVTK(std::ostream &o, int order) const
  {
  
    if (order<=1)
      o << "2 " << _node[0] << " " << _node[1];
    else
    {
      o << "2 " << _node[0] << " " << _node[2];
      for (int i = 0; i < order-2; i++)
      {   
        o << std::endl;
        o << "2 " << _node[2+i] << " " << _node[3+i];
      }
      o << std::endl;
      o << "2 " << _node[order] << " " << _node[1];
    }
          
  }
  
  /** INCOMPLETO: APENAS LINEAR
  * Retorna a tag desta iEdge definida no formato VTK.
  */ 
  static int getCellTypeVTK()
  {
    /* TRABALHO: apenas linear
    */
    return 3; // VTK_LINE (3)
  }
  
  /** Imprime esta iEdge no formato State
  */ 
  void printSelfState(std::ostream &o, int order)
  {
    o << getNodeIdx(0);
    for (int i = 0; i < order+1; ++i)
      o << " " << getNodeIdx(i);
  }
  
  /** Retorna a tag de uma iEdge no formato MSH.
  *  @param order a ordem da iEdge.
  */  
  static int getMSHTag(int order)
  {
    /* poderia fazer meta programção aqui, mas acho complicação desnecessária */
    switch (order)
    {
      case 1: return MSH_LIN_2;
      case 2: return MSH_LIN_3;
      case 3: return MSH_LIN_4;
      case 4: return MSH_LIN_5;
      case 5: return MSH_LIN_6;
      default:
      {
        std::cout << "edge order not supported yet" << std::endl;
        throw;
      }
    }
  }
  
  /** Atualiza esta iEdge para a ordem <em>order</em>.
  *  @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    _node.resize(order+1);
  }
  
  static const int N_borders = 2;
  static const int N_vertices = 2;
  
  /** Destrutor.
  */ 
  ~iEdge() {}
  
protected:
  Fepic::vectori _node;
};



//======================================================================
//============================== iPolygon ==============================
//======================================================================

/** A classe iPolygon representa os polígonos (triângulo, quadrângulos, ...), entidades
 * da malha de dimensão 2. Dependendo da ordem, os polígonos dessa classe podem ter arestas curvas.
 * 
 * @note Esta é uma classe de conceito abstrato (apesar se não ser abstrata sob o ponto
 * de vista da linguagem c++ pois não tem funções virtuais puras), logo não se deve instanciá-la
 * a não ser que se saiba exatamente o que se está fazendo.
 */ 
template<class Traits>
class iPolygon : public iLabel
{
public:
  static const int Dim = 2;
  
  /* Se Polygon está sendo instanciado, então CellT é com certeza uma FaceT*/
  typedef typename Traits::CellT    FaceT;
  typedef typename Traits::MeshT        MeshT;
  
  static const int N_borders  = ElementProperties<FaceT, Traits>::N_borders;
  static const int N_vertices = ElementProperties<FaceT, Traits>::N_vertices;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
  
  /** Construtor.
  */ 
  iPolygon() : iLabel()
  {
    _node.resize(N_vertices);
  }
  
  /** Construtor.
  *      @param nodes vetor com os nós que compõe o iPolygon.
  */
  iPolygon(Fepic::vectorui const& nodes, int label=0) : iLabel(label), _node(nodes)
  {
#ifdef ENTTS_DEBUG
    if (_node.size() < static_cast<uint>(N_vertices))
    {
      std::cout << "error:iPolygon constructor: número insuficiente de nós.\n";
      throw;
    }
#endif
  }

  /** Retorna o número de nós deste polígono.
  */ 
  int getNumNodes() const
  {
    return _node.size();
  }

  /** Retorna o índice do i-ésimo nó do polígono.
  */ 
  uint getNodeIdx(int ith) const
  {   
#ifdef ENTTS_DEBUG
    return _node.at(ith);
#else
    return _node[ith];
#endif
  }

  /** Defini o i-ésimo nó do polígono como nodeid
  */ 
  void setNode(int ith, uint nodeid)
  {
#ifdef ENTTS_DEBUG
    _node.at(ith) = nodeid;
#else
    _node[ith] = nodeid;
#endif
  }
  
  /** Retorna um ponteiro para a i-ésima iHalfEdge.
  */ 
  iHalfEdge<Traits>* getHalf(int ith)
  {
#ifdef ENTTS_DEBUG
    if (ith >= N_borders)
    {
      std::cout << "getHalf: out of range" << std::endl;
      throw;
    }
#endif
    return &_he[ith];
  }

  /** Retorna se os vertices passados formam uma aresta do polígono.
  * @param[in] vertices um vetor com exatamente 2 vertices.
  * @param[out] ith o índice local da aresta que os pontos formam.
  * @return true se e somente se forma uma aresta.
  * @warning a ordem dos nós É importante.
  */
  bool isAnEdge(Fepic::vectorui const& vertices, int &ith) const
  {
          
#ifdef ENTTS_DEBUG
    if (vertices.size() != 2)
    {
      std::cout << "isAnEdge: first argument has wrong dimension.\n" << std::endl;
      throw;
    }
#endif        
    Fepic::vectorui vtx(2);
                
    for (int i = 0; i != N_borders; ++i)
    {
      vtx[0] = _node[i];
      vtx[1] = _node[(i+1)%N_borders];
      
      if (vtx == vertices)
      {
        ith = i;
        return true;
      }
    }
    return false;
  }
        
  /** Retorna se os vertices passados formam uma aresta do polígono.
  * @param[in] vertices um vetor com exatamente 2 vertices.
  * @param[out] ith o índice local da aresta que os pontos formam.
  * @return true se e somente se forma uma aresta.
  * @warning a ordem dos nós NÃO é importante.
  */
  bool isAnParallelEdge(Fepic::vectorui&& vertices, int &ith) const
  {
#ifdef ENTTS_DEBUG
    if (vertices.size() != 2)
    {
      std::cout << "isAnEdge: first argument has wrong dimension.\n" << std::endl;
      throw;
    }
#endif        
    if (this->isAnEdge(vertices, ith))
      return true;
                        
    std::swap(vertices[0], vertices[1]);
    if (this->isAnEdge(vertices, ith))
    {
      std::swap(vertices[0], vertices[1]);
      return true;
    }
    std::swap(vertices[0], vertices[1]);
    return false;
  }

  /** OBSOLETO: refazer\n
  *  Procura por uma edge em comum com outro polígono: se encontra retorna true, caso contrário
  * retorna false.
  */ 
  bool getCommonBorder(FaceT *other, int &thisborder, int &otherborder) const
  {
#ifdef ENTTS_DEBUG
    if (other==NULL)
    {
      std::cout << "getCommonBorder: null FaceT.\n" << std::endl;
      throw;
    }
#endif
    Fepic::vectorui this_edge_vtx(2);
  
    for (int s = 0; s < N_borders; ++s)
    {
      this_edge_vtx = this->getBorderVertices(s);
      reverse(this_edge_vtx.begin(), this_edge_vtx.end());
      if (other->isAnEdge(this_edge_vtx, otherborder))
      {
        thisborder = s;
        return true;
      }
    }
  
    return false;
  }
    
  /** Atribui em cada nó deste polígono uma iHalfEdge compatível.
  * @param mesh a malha em que este polígono está contido.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
#ifdef ENTTS_DEBUG
    if (mesh.edges_local_nodes.empty())
    {
      std::cout << "edges_local_nodes must be initializated." << std::endl;
      throw;
    }
#endif
    for (uint s = 0; s < static_cast<uint>(N_borders); ++s) // loop  nas arestas
    {
      for (uint i = 0, tam=mesh.edges_local_nodes[0].size(); i < tam; ++i)
      {
        mesh.getNode(_node[mesh.edges_local_nodes[s][i]])->setHalf(_he[s]);
      }
    }
  }

  /** Atribui o rótulo deste polígono a seus nós.
  * @param mesh a malha em que este polígono está contido.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem label=0;
  * @warning CONFERIR!
  */ 
  void propagateLabel(MeshT & mesh, bool force=false) const
  {
    if (force)
      for (int i = 0; i < _node.size(); i++)
        mesh.getNode(_node[i])->setLabel(this->getLabel());
    else
      for (int i = 0; i < _node.size(); i++)
        if (mesh.getNode(_node[i])->getLabel() == 0)
          mesh.getNode(_node[i])->setLabel(this->getLabel());
  }
  
  /** Retorna um vetor com os nós da i-ésima edge.
  */  
  Fepic::vectorui getBorderNodes(int ith, MeshT& mesh) const
  {
    Fepic::vectorui nodes;
    
    for (uint i = 0; i < mesh.edges_local_nodes[ith].size(); i++)
      nodes.push_back(_node[ mesh.edges_local_nodes[ith][i] ]);
    return nodes;
  }
  
  /** Retorna um vetor com os vértices da i-ésima edge
  */ 
  Fepic::vectorui getBorderVertices(int ith) const
  {
#ifdef ENTTS_DEBUG
    if (ith >= N_borders)
    {
      std::cout << "getCommonBorder: null FaceT.\n" << std::endl;
      throw;
    }
#endif
    return {_node[ith], _node[(ith+1)%FaceT::N_borders]};
  }
    
  /** Imprime este polígono no formato State (seus nós).
  * @param o a stream onde se vai imprimir, e.g., std::cout.
  */ 
  void printSelfState(std::ostream &o) const
  {
    o << getNodeIdx(0);
    for (uint i = 1; i < _node.size(); ++i)
      o << " " << getNodeIdx(i);
  }

protected:
  Fepic::vectorui     _node;
  iHalfEdge<Traits>   _he[N_borders];
};

//======================================================================
//============================== iTriangle =============================
//======================================================================

/** A classe de triângulos. A dimensão de seu tipo na malha é 2, mas não
 * se deve confundir com a dimensão do espaço onde o triângulo está, que
 * pode ser 2 ou 3. O triângulo tem no mínimo 3 nós (um em cada vértice;
 * triângulo linear). Um triângulo de ordem <em>n</em> tem <em>
 * (n+1)(n+2)/2</em> nós.
 */ 
template<class Traits>
class iTriangle : public iPolygon<Traits>
{
public:
  static const int N_borders = 3; 
  static const int N_vertices = 3;        
        
  typedef Simplex<2>           ElmClass;
  
  typedef UndefElement           VolumeT;  // CellT::dim < 3 ? UndefVol : CellT::Volume
  typedef iTriangle<Traits>        FaceT;

  /* Se iTriangle é instanciado, então ele é a célula */
  typedef iTriangle<Traits>    CellT;
        
  typedef iEdge<Traits>             BorderT;
  typedef typename Traits::PointT   BndBorderT;
  
  typedef iEdge<Traits>             EdgeT;
  typedef typename Traits::PointT   PointT;
  typedef typename Traits::MeshT    MeshT;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
  
  /** Construtor.
  */ 
  iTriangle() : iPolygon<Traits>()
  {
  }
  
  /** Construtor.
  * @param nodes vetor com os nós que compõe o iTriangle.
  */
  iTriangle(Fepic::vectorui const& nodes, int label=0) : iPolygon<Traits>(nodes, label)
  {
  }
  
  /** Retorna a tag do formato MSH correspondente a um triângulo de ordem <em>order</em>.
  */ 
  static int getMSHTag(int order)
  {
    switch (order)
    {
      case 1: return MSH_TRI_3;
      case 2: return MSH_TRI_6;
      case 3: return MSH_TRI_10; // TRABALHO, muito TRABALHO
      case 4: return MSH_TRI_15;
      case 5: return MSH_TRI_21;
      default:
      {
        std::cout << "Triangle order not supported." << std::endl;
        throw;
      }
    }
  }
  
  /** Retorna o número de nós de um triângulo de ordem <em>order</em>
   */ 
  static int getNumNodesCell(int order)
  {
    return (order+1)*(order+2)/2;
  }
    
  /**
  *  Imprime a célula no formate VTK.
  *  @note O número de subdivisões necessários para imprimir o triângulo é order^2
  *  @warning Pula linha no final.
  */ 
  void printSelfVTK(std::ostream &o, int order) const
  {
    static iCellData  data;
    Fepic::matrixi    minimesh;
  
    minimesh = data.getMinimesh(order);
  
    o <<"3 "<<  this->getNodeIdx(minimesh[0][0]) << " " << this->getNodeIdx(minimesh[0][1]) << " " << this->getNodeIdx(minimesh[0][2]);
    for (uint i = 1, tam=minimesh.size(); i < tam; ++i)
    {
      o << std::endl;
      o <<"3 "<<  this->getNodeIdx(minimesh[i][0]) << " " << this->getNodeIdx(minimesh[i][1]) << " " << this->getNodeIdx(minimesh[i][2]);
    }
  }
  
  /** TRABALHO: apenas linear
  * retorna a tag do formato VTK correspondente a este elemento.
  */ 
  static int getCellTypeVTK()
  {
    return 5; // VTK_TRIANGLE(=5)
  }
  
  
  /** Atualiza a ordem deste elemento.
  *  @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    this->_node.resize((order+1)*(order+2)/2);
  }
  
  /**INICIALIZADOR ... NOT FOR USERS!
  */ 
  static Fepic::matrixi getEdgesLocalNodes(int order)
  {
    /* Essa função é usada para atualizar o vetor edges_local_nodes
    * que tem no imesh. Esse vetor contém a numeração local dos nós de
    * uma dada ordem.
    */ 
    
    const int E = order - 1;
    
    Fepic::matrixi en(3, Fepic::vectori(order + 1));
    
    for (int i = 0; i < 3; i++)
    {
      en[i][0] = i;
      en[i][1] = (i+1)%3;
      
      for (int j = 0; j < E; ++j)
        en[i][j+2] = 3+i*E + j;
    }
            
    return en;
  }
  
  /** NOT FO USERS
  * @param n ordem
  */ 
  static Fepic::vectori getOppELN(int n)
  {
    Fepic::vectori opp_eln(n+1);
    
    opp_eln[0]=1;
    opp_eln[n]=0;
    
    for (int i = 0; i < n-1; ++i)
    {
        opp_eln[i+1] = n-i;
    }
    
    return opp_eln;
  }
  
  /** Retorna o número de mini-células de uma mini-malha de ordem n.
  */ 
  static int getNumCellsMM(int n)
  {
    return n*n;
  }
  
  /** Retorna o número de pontos interiores a um triângulo unitário de ordem n
  */ 
  static int getNumBubbles(int n)
  {
    return (n-1)*(n-2)/2;
  }
  
  /** Retorna as coordenadas dos pontos de um triângulo unitário de ordem n
  */ 
  static std::vector<Eigen::Vector2d> getParametricPts(int n)
  {
    return genTriParametricPts(n);
  }
  
  /** NOT FOR USERS \n 
  * classe auxiliar para armazenar dados comuns a todos os objetos
  * da classe iTriangle<>
  *  
  */ 
  class iCellData
  {
  public:
    typedef Fepic::matrixi Minimesh;
    typedef int          Order;
    
    /** NOT FOR USERS \n
    *  Dado uma célula de ordem n, essa função retorna uma mini-malha que
    * subdivide essa célula. Essa função é últil para imprimir a célula em VTK.
    */ 
    Minimesh getMinimesh(Order n)
    {
      std::map<Order, Minimesh>::iterator it = table.find(n);
      if (it != table.end()) // se já existe esta tabelada esta mini-malha
        return it->second;
      else // se não, cria esta mini-malha e guarda na tabela
      {
        /* encontrar todas as mini-células que tenham o padrão:
        * - (a,b), (a+1,b), (a,b+1)
        * - (a,b), (a,b+1), (a-1,b+1)
        */ 
        Fepic::matrixi                  minimesh;       // mini-malha, inicialmente vazia
        int                           n2, n3;         // id do segundo e terceiro nó na mini-malha
        int                           a,b;            // coordenada inteira
        std::vector<Eigen::Vector2i>  coords_list = genTriParametricPtsINT(n);
        auto                          clbegin = coords_list.begin(), clend = coords_list.end();
        decltype(clbegin)             clit2, clit3;   // iteradores
                                    
        
        for (int i = 0; i < (int)coords_list.size(); ++i)
        {
          a = coords_list[i](0);
          b = coords_list[i](1);
          
          /* Primeiro padrão: (a,b), (a+1,b), (a,b+1)  */
          clit2 = find(clbegin, clend, Eigen::Vector2i(a+1,b));
          clit3 = find(clbegin, clend, Eigen::Vector2i(a,b+1));
          
          if ((clit2 != clend) && (clit3 != clend))
          {
            n2 = distance(clbegin, clit2);
            n3 = distance(clbegin, clit3);
            
            minimesh.push_back(Fepic::vectori{i,n2,n3});
          }
          
          /* Segundo padrão: (a,b), (a,b+1), (a-1,b+1)  */
          clit2 = find(clbegin, clend, Eigen::Vector2i(a,b+1));
          clit3 = find(clbegin, clend, Eigen::Vector2i(a-1,b+1));
          
          if ((clit2 != clend) && (clit3 != clend))
          {
            n2 = distance(clbegin, clit2);
            n3 = distance(clbegin, clit3);
            
            minimesh.push_back(Fepic::vectori{i,n2,n3});
          }
          
        } // end for
        
        table[n] = minimesh;
        
        return minimesh;
        
      } // end if
      
    } // end getMinimesh
    
    
    // atributos
    std::map<Order, Minimesh> table;
  }; // iCellData
  
  /** Destrutor.
  */ 
  ~iTriangle() {}
  
protected:

};


//======================================================================
//============================== iPolyhedron ===========================
//======================================================================

template<class Traits>
class iPolyhedron : public iLabel
{
public:
  static const int Dim = 3;
  
  /* Se Polygon está sendo instanciado, então CellT é com certeza uma FaceT*/
  typedef typename Traits::CellT  CellT;
  typedef typename ElementProperties<CellT, Traits>::FaceT FaceT;
  typedef typename Traits::MeshT      MeshT;
  
  static const int N_borders  = ElementProperties<CellT, Traits>::N_borders;
  static const int N_vertices = ElementProperties<CellT, Traits>::N_vertices;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
    
  iPolyhedron() : iLabel()
  {
    _node.resize(N_vertices);
  }
    
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  iPolyhedron(Fepic::vectorui const& nodes, int label=0) : iLabel(label), _node(nodes)
  {
#ifdef ENTTS_DEBUG
    if (_node.size() < static_cast<uint>(N_vertices))
    {
      std::cout << "error:iPolyhedron constructor: número insuficiente de nós.\n";
      throw;
    }
#endif
  }

  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  iPolyhedron(Fepic::vectorui&& nodes, int label=0) : iLabel(label), _node(std::move(nodes))
  {
#ifdef ENTTS_DEBUG
    if (_node.size() < N_vertices)
    {
      std::cout << "error:iPolyhedron constructor: número insuficiente de nós.\n";
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
#ifdef ENTTS_DEBUG
    return _node.at(ith);
#else
    return _node[ith];
#endif
  }

  /** Defini o i-ésimo nó do polígono como nodeid
  */ 
  void setNode(int ith, uint nodeid)
  {
#ifdef ENTTS_DEBUG
    _node.at(ith) = nodeid;
#else
    _node[ith] = nodeid;
#endif
  }

  /** Retorna a i-ésima half-edge.
  */ 
  iHalfFace<Traits>* getHalf(int ith)
  {
#ifdef ENTTS_DEBUG
    if (ith >= N_borders)
    {
      std::cout << "getHalf: out of range" << std::endl;
      throw;
    }
#endif
    return &_hf[ith];
  }
  
  /** Retorna se os vértices passados formam uma face nesse poliedro.
  * @param[in] nodes um vetor com os vértices. O tamanho do vetor deve ser exatamente o número
  *            de lados (N_borders) do elemento degenerado deste poliedro.
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
    if (int(vertices.size()) != FaceT::N_borders)
    {
      std::cout << "isAnFace error: wrong parameter!" << std::endl;
      std::cout << "vertices.size() = " << vertices.size() << std::endl;
    }
    
    Fepic::vectorui temp(FaceT::N_borders);
    Fepic::vectorui this_vtx(FaceT::N_borders);
    
    for (int f = 0; f < N_borders; ++f) // para cada face
    {
      for (int n = 0; n < FaceT::N_borders; n++)
      {
        this_vtx[n] = _node[ faces_vtx[f][n] ];
      }                       
      for (int a = 0; a < FaceT::N_borders; ++a) // para cada ancora
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
    
    for (uint i = 0, tam=mesh.faces_local_nodes[ith].size(); i < tam; i++)
      nodes.push_back( _node[ mesh.faces_local_nodes[ith][i] ] );
      
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
    for (int f = 0; f < N_borders; ++f) // loop  nas faces
      for (uint i = 0, tam=mesh.faces_local_nodes[0].size(); i < tam; ++i)
        mesh.getNode(_node[mesh.faces_local_nodes[f][i]])->setHalf(_hf[f]);
  }
  
  /** Atribui o o label do tetraedro em seus nós.
  * @param force Caso true, a herança do label será feita sem restrição. Caso false,
  * a herança será feita apenas se cada nó não tiver label (i.e., label = 0).
  * @warning CONFERIR!!
  */ 
  void propagateLabel(MeshT & mesh, bool force=false) const
  {
    if (force)
      for (int i = 0; i < this->node.size(); i++)
        mesh.getNode(_node[i])->setLabel(this->getLabel());
    else
      for (int i = 0; i < this->node.size(); i++)
        if (mesh.getNode(_node[i])->getLabel() == 0)
          mesh.getNode(_node[i])->setLabel(this->getLabel());
  }
  
  
protected:
        Fepic::vectorui                 _node;
        iHalfFace<Traits>       _hf[N_borders];    
};



//======================================================================
//============================== iTetrahedron ==========================
//======================================================================

template<class Traits>
class iTetrahedron : public iPolyhedron<Traits> {
public:
  static const int Dim = 3;
  
  typedef          Simplex<3>           ElmClass;
  typedef          iTetrahedron<Traits> VolumeT;
  typedef          iTriangle<Traits>    FaceT;
  typedef          FaceT                BorderT;
  typedef          iEdge<Traits>        BndBorderT;
  typedef          BndBorderT           EdgeT;
  typedef typename Traits::PointT       PointT;       
  typedef typename Traits::MeshT        MeshT;  
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1>  VecT;      
  
  iTetrahedron() : iPolyhedron<Traits>()
  {
  }
  
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  iTetrahedron(Fepic::vectorui const& nodes, int label=0) : iPolyhedron<Traits>(nodes, label)
  {
  }
  
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  iTetrahedron(Fepic::vectorui && nodes, int label=0) : iPolyhedron<Traits>(nodes, label)
  {
  }
  
  /** Retorna o número de nós desse poliedro.
  */ 
  static int getNumNodesCell(int order)
  {
    return (order+1)*(order+2)*(order+3)/6;;
  }
                  
  void printSelfVTK(std::ostream &o, int order) const
  {
    static iCellData data;
    Fepic::matrixi   minimesh;
  
    minimesh = data.getMinimesh(order);
    
    o <<"4 "<<  this->getNodeIdx(minimesh[0][0]) << " " << this->getNodeIdx(minimesh[0][1]) << " " << this->getNodeIdx(minimesh[0][2]) << " " << this->getNodeIdx(minimesh[0][3]);
    for (uint i = 1, tam=minimesh.size(); i < tam; ++i)
    {
      o << std::endl;
      o <<"4 "<<  this->getNodeIdx(minimesh[i][0]) << " " << this->getNodeIdx(minimesh[i][1]) << " " << this->getNodeIdx(minimesh[i][2]) << " " << this->getNodeIdx(minimesh[i][3]);
    }
          
  }
  
  /* OBS: apenas linear
  */ 
  static int getCellTypeVTK()
  {
    return 10; // VTK_TETRA (=10)
  }       
  
  static int getMSHTag(int order)
  {
    switch (order)
    {
      case 1: return MSH_TET_4;
      case 2: return MSH_TET_10;
      case 3: return MSH_TET_20;
      case 4: return MSH_TET_35;
      case 5: return MSH_TET_56;
      default:
      {
        std::cout << "invalid tetrahedron order." << std::endl;
        throw;
      }
    }
  }
  
  
  /** @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    this->_node.resize((order+1)*(order+2)*(order+3)/6);
  }
  
  
  /**INICIALIZADOR ... NÃO UTILIZAR
  */ 
  static Fepic::matrixi getEdgesLocalNodes(int order)
  {
    static Fepic::matrixi edges_vtx = ElementProperties<iTetrahedron, Traits>::get_edges_vtx();
  
    const int E = order - 1;
    
    Fepic::matrixi en(6, std::vector<int>(order + 1));
    
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 2; ++j)
        en[i][j] = edges_vtx[i][j];
      
      for (int j = 0; j < E; ++j)
        en[i][j+2] = 4 + i*E + j;
    }
            
    return en;
  }
  
  /**INICIALIZADOR ... NÃO UTILIZAR
  */ 
  static Fepic::matrixi getFacesLocalNodes(int order)
  {       
    static Fepic::matrixi faces_vtx(ElementProperties<iTetrahedron, Traits>::get_faces_vtx());
  
    const int E = order - 1;
    const int C = 4 + E*6;
    const int F = (order-1)*(order-2)/2;
    
    Fepic::matrixi fn(4, Fepic::vectori((order+1)*(order+2)/2));
    
    Fepic::matrixi ed_nds = iTetrahedron::getEdgesLocalNodes(order);
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        fn[i][j] = faces_vtx[i][j];
    
    /*
    *  face | edges         (" - ": orientação invertida) 
    * ¯¯¯¯¯¯|¯¯¯¯¯¯¯¯¯¯¯¯
    *        0  | -0 -2 -1       
    *        1  | +0 -5 +3       
    *        2  | +4 +2 -3       
    *        3  | -4 +5 +1       
    */
      
    for (int i = 0; i < E; i++)
    {
      fn[0][i+3]     = ed_nds[0][order - i];
      fn[0][i+3+E]   = ed_nds[2][order - i];
      fn[0][i+3+2*E] = ed_nds[1][order - i];
      
      fn[1][i+3]     = ed_nds[0][2+i];
      fn[1][i+3+E]   = ed_nds[5][order - i];
      fn[1][i+3+2*E] = ed_nds[3][2+i];
      
      fn[2][i+3]     = ed_nds[4][2+i];
      fn[2][i+3+E]   = ed_nds[2][2+i];
      fn[2][i+3+2*E] = ed_nds[3][order - i];
      
      fn[3][i+3]     = ed_nds[4][order - i];
      fn[3][i+3+E]   = ed_nds[5][2+i];
      fn[3][i+3+2*E] = ed_nds[1][2+i];
    }
    
    for (int f = 0; f < 4; f++)
      for (int i = 0; i < F; i++)
        fn[f][i+3*(1+E)] = C+f*F + i;
    
    return fn;
    
  }
  
  /** NOT FOR USERS
  * @param n ordem
  */ 
  static Fepic::vectori getOppELN(int n)
  {
    Fepic::vectori opp_eln(n+1);
    
    opp_eln[0]=1;
    opp_eln[n]=0;
    
    for (int i = 0; i < n-1; ++i)
      opp_eln[i+1] = n-i;
    
    return opp_eln;
  }
  
  /** NOT FOR USERS \n
  * getOppFLN(n)[anchor] = vetor com a visão da face oposta de âncora anchor.
  * @param n ordem
  */ 
  static Fepic::matrixi getOppFLN(int n)
  {
    int k;
    int anch;
    const int                           N = (n+1)*(n+2)/2;
    Fepic::matrixi                      opp_fln(3, Fepic::vectori(N));
    const std::vector<Eigen::Vector2i>  Coords = genTriParametricPtsINT(n);
    
    std::map<int, Eigen::Vector2i>      X;
    
    Eigen::Matrix2i A;
    Eigen::Vector2i b;
  
#define PAIR std::pair<int, Eigen::Vector2i>

/* transformações são feitas fazendo A*x + b */
#define DO_TRANSFORMATION_AND_APPLY \
    for (int i = 0; i != N; ++i) \
      X[i] = A*Coords[i] + b;    \
                                 \
    for (int i = 0; i != N; ++i) \
    {                            \
      k = find_if(X.begin(), X.end(), [i, &Coords](PAIR const& p)->bool {  \
        return (p.second==Coords[i]);   \
      })->first;                        \
                                        \
      opp_fln[anch][i] = k;             \
    }                                                                   
        
    /* ----------------
    *    ANCORA 0    
    * ---------------- */
    anch = 0;
    A << +1, +0,
         -1, -1;
    b << +0,
         +n;
        
    DO_TRANSFORMATION_AND_APPLY;
    
    /* ----------------
    *    ANCORA 1    
    * ---------------- */
    anch = 1;
    A << -1, -1,
         +0, +1;
    b << +n,
         +0;
        
    DO_TRANSFORMATION_AND_APPLY;
    
    /* ----------------
    *    ANCORA 0    
    * ---------------- */
    anch = 2;
    A << +0, +1,
         +1, +0;
    b << +0,
         +0;
    
        
    DO_TRANSFORMATION_AND_APPLY;
    
#undef DO_TRANSFORMATION_AND_APPLY
#undef PAIR

      return opp_fln;
    }

  /** Retorna o número de mini-células de uma mini-malha de ordem n.
  */ 
  static int getNumCellsMM(int n)
  {
    return n*n*n;
  }
  
  /** NOT FOR USERS \n 
  * classe auxiliar para armazenar dados comuns a todos os objetos
  * da classe iTriangle<>
  *  
  */ 
  class iCellData
  {
  public:
    typedef Fepic::matrixi Minimesh;
    typedef int            Order;
  
    /** NOT FOR USERS \n
    *  Dado uma célula de ordem n, essa função retorna uma mini-malha que
    * subdivide essa célula. Essa função é últil para imprimir a célula em VTK.
    */ 
    Minimesh getMinimesh(Order n)
    {
      std::map<Order, Minimesh>::iterator it = table.find(n);
      if (it != table.end()) // se já existe esta tabelada esta mini-malha
      {
        return it->second;
      }
      else // se não, cria esta mini-malha e guarda na tabela
      {
        /* encontrar todas as mini-células que tenham os padrões:
        * - (a,b,c), (a+1,b+0,c+0), (a+0,b+1,c+0), (a+0,b+0,c+1)
        * - (a,b,c), (a+1,b+0,c-1), (a+0,b+1,c-1), (a+0,b+1,c+0)
        * - (a,b,c), (a+0,b+0,c+1), (a+1,b-1,c+0), (a+1,b+0,c+0)
        * - (a,b,c), (a+1,b+0,c-1), (a+0,b+1,c+0), (a+1,b+0,c+0)
        * - (a,b,c), (a+0,b+1,c+0), (a-1,b+1,c+1), (a+0,b+0,c+1)
        * - (a,b,c), (a+1,b+0,c-1), (a+1,b+0,c+0), (a+1,b-1,c+0)
        */ 
        Fepic::matrixi                  minimesh;       // mini-malha, inicialmente vazia
        int                           n2, n3, n4;     // id do segundo e terceiro nó na mini-malha
        int                           a,b,c;          // coordenada inteira
        std::vector<Eigen::Vector3i>  coords_list = genTetParametricPtsINT(n);
        auto                          clbegin = coords_list.begin(), clend = coords_list.end();
        decltype(clbegin)             clit2, clit3, clit4;   // iteradores
                                    
        
        for (int i = 0; i < (int)coords_list.size(); ++i)
        {
            a = coords_list[i](0);
            b = coords_list[i](1);
            c = coords_list[i](2);
            
            /* Primeiro padrão: (a,b,c), (a+1,b+0,c+0), (a+0,b+1,c+0), (a+0,b+0,c+1)  */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c+0));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a+0,b+1,c+0));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+0,b+0,c+1));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            /* Segundo padrão: (a,b,c), (a+1,b+0,c-1), (a+0,b+1,c-1), (a+0,b+1,c+0)  */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c-1));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a+0,b+1,c-1));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+0,b+1,c+0));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            /* Terceiro padrão: (a,b,c), (a+0,b+0,c+1), (a+1,b-1,c+0), (a+1,b+0,c+0)  */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+0,b+0,c+1));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a+1,b-1,c+0));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c+0));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            /* Quarto padrão: (a,b,c), (a+1,b+0,c-1), (a+0,b+1,c+0), (a+1,b+0,c+0)  */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c-1));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a+0,b+1,c+0));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c+0));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            /* Quinto padrão: (a,b,c), (a+0,b+1,c+0), (a-1,b+1,c+1), (a+0,b+0,c+1)  */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+0,b+1,c+0));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a-1,b+1,c+1));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+0,b+0,c+1));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            /* Sexto padrão: (a,b,c), (a+1,b+0,c-1), (a+1,b+0,c+0), (a+1,b-1,c+0)   */
            clit2 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c-1));
            clit3 = find(clbegin, clend, Eigen::Vector3i(a+1,b+0,c+0));
            clit4 = find(clbegin, clend, Eigen::Vector3i(a+1,b-1,c+0));
            
            if ((clit2 != clend) && (clit3 != clend) && (clit4 != clend))
            {
              n2 = distance(clbegin, clit2);
              n3 = distance(clbegin, clit3);
              n4 = distance(clbegin, clit4);
              
              minimesh.push_back(Fepic::vectori{i,n2,n3,n4});
            }
            
            
          } // end for
          
          table[n] = minimesh;
          return minimesh;
          
      } // end if
      
    } // end getMinimesh
  
    /* iCellData members */
    std::map<Order, Minimesh> table;
    
  }; // class iCellData
      
  static std::string name()
  {
    return "Tetrahedron";
  }
    
protected:

        
};





//======================================================================
//----------------------------------------------------------------------
//------------------------------ TOOLS ---------------------------------
//----------------------------------------------------------------------
//======================================================================




/** versão estática.\n
 *  Faz um mapeamento de pontos na célula unitário para a célula real
 *  @param list_pts uma lista com pontos no triângulo unitário nos quais se deseja fazer o mapeamento
 *  @param intp_pts a lista de pontos de interpolação da célula; ex, para interpolação linear, passar os vértices
 *  @param Phi funções de interpolação que correspondem aos pontos de interpolação
 *  @return uma lista com as coordenadas dos pontos passados em list_pts na célula real
 *  @warning as funções Phi DEVEM corresponder aos pontos de interpolação
 */
template<class Traits, class ShapeFun, 
         int   sdim = Traits::spacedim,
         int   cdim = Traits::CellT::Dim,
         class VecT = Eigen::Matrix<double, sdim, 1>,   // vetor no espaço da célula real
         class VecU = Eigen::Matrix<double, cdim, 1> >  // vetor na espaço da célula unitária
std::vector<VecT> map2RealCell(std::vector<VecU> const& list_pts,
                               std::vector<VecT> const& intp_pts,
                               ShapeFun          const& Phi)
{
  int tam = static_cast<int>( list_pts.size() );
  std::vector<VecT> ret(tam, VecT::Zero());
  
  for (int i = 0; i < tam; ++i)
    for (int k = 0, size=static_cast<int>( intp_pts.size() ); k < size; ++k)
      ret[i] += Phi(list_pts[i], k)*intp_pts[k];
  
  return ret;
}










#undef ENTTS_DEBUG



#endif





