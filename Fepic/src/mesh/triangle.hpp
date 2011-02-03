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

#ifndef FEPIC_TRIANGLE_HPP
#define FEPIC_TRIANGLE_HPP




/** A classe de triângulos. A dimensão de seu tipo na malha é 2, mas não
 * se deve confundir com a dimensão do espaço onde o triângulo está, que
 * pode ser 2 ou 3. O triângulo tem no mínimo 3 nós (um em cada vértice;
 * triângulo linear). Um triângulo de ordem <em>n</em> tem <em>
 * (n+1)(n+2)/2</em> nós.
 */ 
template<class Traits>
class Triangle : public _Poly2d<Traits>
{
public:
  static const int n_borders = 3; 
  static const int n_vertices = 3;        
        
  typedef Simplex<2>           ElmClass;
  
  typedef UndefElement           VolumeT;  // CellT::dim < 3 ? UndefVol : CellT::Volume
  typedef Triangle<Traits>        FaceT;

  /* Se Triangle é instanciado, então ele é a célula */
  typedef Triangle<Traits>    CellT;
        
  typedef Edge<Traits>             BorderT;
  typedef typename Traits::PointT   BndBorderT;
  
  typedef Edge<Traits>             EdgeT;
  typedef typename Traits::PointT   PointT;
  typedef typename Traits::MeshT    MeshT;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
  
  /** Construtor.
  */ 
  Triangle() : _Poly2d<Traits>()
  {
  }
  
  /** Construtor.
  * @param nodes vetor com os nós que compõe o Triangle.
  */
  Triangle(Fepic::vectorui const& nodes, int label=0) : _Poly2d<Traits>(nodes, label)
  {
  }
  
  /** Retorna a tag do formato Msh correspondente a um triângulo de ordem <em>order</em>.
  */ 
  static int getMshTag(int order)
  {
    switch (order)
    {
      case 1: return Msh_TRI_3;
      case 2: return Msh_TRI_6;
      case 3: return Msh_TRI_10; // TRABALHO, muito TRABALHO
      case 4: return Msh_TRI_15;
      case 5: return Msh_TRI_21;
      default:
      {
        std::cout << "Triangle order not supported." << std::endl;
        throw;
      }
    }
  }
  
  /** Retorna o número de nós de um triângulo de ordem <em>order</em>
   */ 
  static int getNumCellNodes(int order)
  {
    return (order+1)*(order+2)/2;
  }
    
  /**
  *  Imprime a célula no formate Vtk.
  *  @note O número de subdivisões necessários para imprimir o triângulo é order^2
  *  @warning Pula linha no final.
  */ 
  void printSelfVtk(std::ostream &o, int order) const
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
  * retorna a tag do formato Vtk correspondente a este elemento.
  */ 
  static int getCellTypeVtk()
  {
    return 5; // Vtk_TRIANGLE(=5)
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
  * da classe Triangle<>
  *  
  */ 
  class iCellData
  {
  public:
    typedef Fepic::matrixi Minimesh;
    typedef int          Order;
    
    /** NOT FOR USERS \n
    *  Dado uma célula de ordem n, essa função retorna uma mini-malha que
    * subdivide essa célula. Essa função é últil para imprimir a célula em Vtk.
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
  ~Triangle() {}
  
protected:

};





#endif
