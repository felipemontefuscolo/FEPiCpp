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

#ifndef FEPIC_TETRAHEDRON_HPP
#define FEPIC_TETRAHEDRON_HPP



template<class Traits>
class Tetrahedron : public _Poly3d<Traits> {
public:
  static const int Dim = 3;
  
  typedef          Simplex<3>           ElmClass;
  typedef          Tetrahedron<Traits> VolumeT;
  typedef          Triangle<Traits>    FaceT;
  typedef          FaceT                BorderT;
  typedef          Edge<Traits>        BndBorderT;
  typedef          BndBorderT           EdgeT;
  typedef typename Traits::PointT       PointT;       
  typedef typename Traits::MeshT        MeshT;  
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1>  VecT;      
  
  Tetrahedron() : _Poly3d<Traits>()
  {
  }
  
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  Tetrahedron(Fepic::vectorui const& nodes, int label=0) : _Poly3d<Traits>(nodes, label)
  {
  }
  
  /** @param nodes vetor com os nós que compõe o polígono.
  */ 
  Tetrahedron(Fepic::vectorui && nodes, int label=0) : _Poly3d<Traits>(nodes, label)
  {
  }
  
  /** Retorna o número de nós desse poliedro.
  */ 
  static int getNumCellNodes(int order)
  {
    return (order+1)*(order+2)*(order+3)/6;;
  }
                  
  void printSelfVtk(std::ostream &o, int order) const
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
  static int getCellTypeVtk()
  {
    return 10; // Vtk_TETRA (=10)
  }       
  
  static int getMshTag(int order)
  {
    switch (order)
    {
      case 1: return Msh_TET_4;
      case 2: return Msh_TET_10;
      case 3: return Msh_TET_20;
      case 4: return Msh_TET_35;
      case 5: return Msh_TET_56;
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
    static Fepic::matrixi edges_vtx = ElementProperties<Tetrahedron, Traits>::get_edges_vtx();
  
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
    static Fepic::matrixi faces_vtx(ElementProperties<Tetrahedron, Traits>::get_faces_vtx());
  
    const int E = order - 1;
    const int C = 4 + E*6;
    const int F = (order-1)*(order-2)/2;
    
    Fepic::matrixi fn(4, Fepic::vectori((order+1)*(order+2)/2));
    
    Fepic::matrixi ed_nds = Tetrahedron::getEdgesLocalNodes(order);
    
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
  * da classe Triangle<>
  *  
  */ 
  class iCellData
  {
  public:
    typedef Fepic::matrixi Minimesh;
    typedef int            Order;
  
    /** NOT FOR USERS \n
    *  Dado uma célula de ordem n, essa função retorna uma mini-malha que
    * subdivide essa célula. Essa função é últil para imprimir a célula em Vtk.
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








#endif
