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
template<class _Traits>
class Triangle : public _Labelable, public _Poly2d<_Traits>
{
public:

  typedef          Simplex<2>      PolytopeT;
  typedef typename _Traits::PointT PointT;
  typedef typename _Traits::HalfT  HalfT;      
  typedef typename _Traits::MeshT  MeshT;  
  
  friend class _Poly2d<_Traits>;
  friend class _CellCore<_Traits>;
  
  enum { dim=2,
         n_vertices=3,
         n_borders=3,
         n_vertices_per_border=2};
  
  template<class... LabeableArgs>
  Triangle(Eigen::VectorXi const& nodes, int order, LabeableArgs... args) :
                          _Labelable(args...), _nodes(nodes), _order(order)
  {
    FEPIC_CHECK(nodes.size()==numNodes<Simplex<2>>(order), "", std::invalid_argument);
  }
  
  template<class... LabeableArgs>
  Triangle(Eigen::VectorXi && nodes, int order, LabeableArgs... args) :
                      _Labelable(args...), _nodes(nodes), _order(order)
  {
    FEPIC_CHECK(nodes.size()==numNodes<Simplex<2>>(order), "", std::invalid_argument);
  }  
  
  Triangle() : _nodes(Eigen::Vector3i(-1,-1,-1)), _order(1) {};
  Triangle(Triangle const&) = default;
  ~Triangle() = default;

    /** Atualiza a ordem deste elemento.
  *  @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    this->_nodes.resize((order+1)*(order+2)/2);
    _order = order;
  }
  
  /** NOT FO USERS
  * @param n ordem
  */ 
  static Eigen::VectorXi getOppELN(int n)
  {
    FEPIC_CHECK(n>0, "", std::invalid_argument);
    Eigen::VectorXi opp_eln(n+1);
    
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
  static int getNumSubdivisions(int n)
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

  static Eigen::Matrix3Xi const& getBordersLocalNodes(int order)
  {
    FEPIC_CHECK(order<=MAX_CELLS_ORDER,  "order not supported yet", std::out_of_range);
    return Triangle::table1[order];
  }

  static Eigen::MatrixX3i const& getMinimesh(int order)
  {
    FEPIC_CHECK(order<=MAX_CELLS_ORDER, "order not supported yet", std::out_of_range);
    return Triangle::table2[order];
  }

protected:

  static Eigen::Matrix3Xi _getBordersLocalNodes(int order)
  {
    const int E = order - 1;
    
    Eigen::Matrix3Xi en(3, order+1);
    
    for (int i = 0; i < 3; i++)
    {
      en(i,0) = i;
      en(i,1) = (i+1)%3;
      
      for (int j = 0; j < E; ++j)
        en(i,j+2) = 3+i*E + j;
    }
            
    return en;    
  }

  static Eigen::MatrixX3i _getMinimesh(int order)
  {
    /* encontrar todas as mini-células que tenham o padrão:
    * - (a,b), (a+1,b), (a,b+1)
    * - (a,b), (a,b+1), (a-1,b+1)
    */
    Eigen::MatrixX3i minimesh(order*order, 3);
    
    int n2, n3; // id do segundo e terceiro nó na mini-malha
    int a,b;    // coordenada inteira (a,b)
    std::vector<Eigen::Vector2i> coords_list = genTriParametricPtsINT(order);
    std::vector<Eigen::Vector2i>::iterator clbegin = coords_list.begin();
    std::vector<Eigen::Vector2i>::iterator clend  = coords_list.end();
    std::vector<Eigen::Vector2i>::iterator clit2, clit3;
                                
    int counter(0);
    for (int i = 0; i < static_cast<int>(coords_list.size()); ++i)
    {
      a = coords_list[i](0);
      b = coords_list[i](1);
      
      /* Primeiro padrão: (a,b), (a+1,b), (a,b+1) */
      clit2 = find(clbegin, clend, Eigen::Vector2i(a+1,b));
      clit3 = find(clbegin, clend, Eigen::Vector2i(a,b+1));
      
      if ((clit2 != clend) && (clit3 != clend))
      {
        n2 = distance(clbegin, clit2);
        n3 = distance(clbegin, clit3);
        
        minimesh.row(counter++) = Eigen::Vector3i(i,n2,n3);
      }
      
      /* Segundo padrão: (a,b), (a,b+1), (a-1,b+1) */
      clit2 = find(clbegin, clend, Eigen::Vector2i(a,b+1));
      clit3 = find(clbegin, clend, Eigen::Vector2i(a-1,b+1));
      
      if ((clit2 != clend) && (clit3 != clend))
      {
        n2 = distance(clbegin, clit2);
        n3 = distance(clbegin, clit3);
        
        minimesh.row(counter++) = Eigen::Vector3i(i,n2,n3);
      }
      
    } // end for
    
    return minimesh;
  } // end getMinimesh 

  static std::vector<Eigen::Matrix3Xi> _initTable1()
  {
    // the table1 at index n store the local nodes of borders of a cell
    // with order n.
    // WARNING: must be constructed after _initTable0 be called
    std::vector<Eigen::Matrix3Xi> table;
    
    table.reserve(MAX_CELLS_ORDER+1);
    table.push_back(Triangle::_getBordersLocalNodes(1)); // no order 0,plz
    for (int i = 1; i <= MAX_CELLS_ORDER; ++i)
    {
      table.push_back(Triangle::_getBordersLocalNodes(i));
    }
    
    return table;    
  }

  static std::vector<Eigen::MatrixX3i> _initTable2()
  {
    // the table2 at index n store the minimesh of a cell with order n.
    std::vector<Eigen::MatrixX3i> table;
    
    table.reserve(MAX_CELLS_ORDER+1);
    table.push_back(Triangle::_getMinimesh(1)); // no order 0,plz
    for (int i = 1; i <= MAX_CELLS_ORDER; ++i)
    {
      table.push_back(Triangle::_getMinimesh(i));
    }
    
    return table;    
  }
  
public:  
  static const Eigen::Matrix32i borders_local_vertices; // edges

  static const std::vector<Eigen::Matrix3Xi> table1; // order / borders local nodes
  static const std::vector<Eigen::MatrixX3i> table2; // order / minimeshs

protected:
  unsigned char   _order;
  HalfT           _halfs[n_borders];
  Eigen::VectorXi _nodes;
  
};

template<class _Traits>
const Eigen::Matrix32i Triangle<_Traits>::borders_local_vertices = ((Eigen::Matrix32i()<< 0,1,  1,2,  2,0).finished());

template<class _Traits>
const std::vector<Eigen::Matrix3Xi> Triangle<_Traits>::table1 = Triangle::_initTable1();

template<class _Traits>
const std::vector<Eigen::MatrixX3i> Triangle<_Traits>::table2 = Triangle::_initTable2();




#endif
