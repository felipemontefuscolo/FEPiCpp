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

#ifndef SHAPEFUNCTIONS_HPP
#define SHAPEFUNCTIONS_HPP

#include "Fepic/src/shapefunctions/parametric_pts.hpp"
#include <algorithm>
#include "Fepic/src/util/misc.hpp"

/* Computa:  product(n*x-k, k, 0, Q-1)
 * É a definição que está na documentação.
 */
/** NOT FOR USERS
 */
inline
double _omega(int n, int Q, double x)
{
  double prod=1;
  for (int k = 0; k < Q; k++)
  {
    prod *= n*x-k;
  }
  return prod;
}

/* Computa:  ddx( product(n*x-k, k, 0, Q-1)  )
 * É a definição que está na documentação.
 */
/** NOT FOR USERS
 */
inline
double _ddxomega(int n, int Q, double x)
{
  if (Q==0)
    return 0.;
  return (n*x - Q + 1.)*_ddxomega(n, Q-1, x) + n*_omega(n, Q-1, x);
}

/*
 * Computa: product(n-3*k,k,0,Q-1)
 */
/** NOT FOR USERS
 */
inline
long int _prod3(int n, int Q)
{
  long int prod=1;
  for (int k = 0; k < Q; k++)
    prod *= n-3*k;
  return prod;
}

/*
 * Computa: product(n-4*k,k,0,Q-1)
 */
/** NOT FOR USERS
 */
inline
long int _prod4(int n, int Q)
{
  long int prod=1;
  for (int k = 0; k < Q; k++)
    prod *= n-4*k;
  return prod;
}


template<class T> class ShapeFunction;


/*
 *
 *
 *    SEGMENT
 *
 *  FUNCIONANDO
 */
template<>
class ShapeFunction<Simplex<1> >
{
public:
  ShapeFunction(int degree=1) : _degree(degree)
  {
  };

  void setDegree(int degree)
  {
    _degree=(degree);
  }


  /** Retorna a função avaliada nas coordenadas do segmento unitário [-1,1]
   */
  double operator() (double qsi, int ith) const
  {
    if (_degree==0)
    {
      return 1.;
    }
    if (_degree==1)
    {
      switch (ith)
      {
        case 0: return 1.-qsi;
        case 1: return qsi;
      }
    }

    double prod=1;
        double xk;
        double xi = (2.*ith-_degree)/_degree;
    for (int k = 0; k < _degree+1; ++k)
      if (k!=ith)
      {
        xk = (2.*k-_degree)/_degree;
        prod *= (qsi - xk)/(xi - xk);
      }

    return prod;
  }

  /** component = 0 => componentada com respeito a L1
   *  component = 1 => componentada com respeito a L2
   *  outro valor: error
   */
  double gradL(double qsi, int ith) const
  {

    if (_degree == 0)
      return 0;

    if (_degree == 1)
    {
      switch (ith)
      {
        case 0:
          return -1./2;
        case 1:
          return 1./2;

        default:
          std::cout << "gradL: invalid ith" << std::endl;
          throw;
      }
    }

    double denom = 1;
    double prod;
    double sum=0;
    double xl,xk,xi = (2.*ith-_degree)/_degree;

    for (int k = 0; k < _degree+1; ++k)
      if (k!=ith)
      {
        xk = (2.*k-_degree)/_degree;
        prod=1;
        for (int l = 0; l < _degree+1; ++l)
        {
          if ((l!=ith) && (l!=k))
          {
            xl = (2.*l-_degree)/_degree;
            prod *= (qsi - xl);
          }
        }
        sum += prod;
        denom *= (xi-xk);
      }

    return sum/denom;
  }

  int getDegree() const
  {
    return _degree;
  }

  int getNumDof() const
  {
    return _degree+1;
  }

  /* member functions */
  int _degree;
};

/*
 *    (0,1)
 *      |`\
 *      |  `\
 *      |    `\
 *      |      `\
 *      |        `\
 *      |__________`\
 *    (0,0)        (1,0)
 *
 */
/**
 *  Cria funções de interpolação de Lagrange no triângulo unitário, podendo
 *  ser enriquecidos com função bolha (hierárquica ou não).
 */
template<>
class ShapeFunction<Simplex<2> >
{
public:
  /**
   *  @param degree o grau de intepolação
   *  @param flags 8 bits: \n
   *                 bit0 se com bolha \n
   *           bit1 se com bolha hierárquica
   *
   *  @note se setar o bit1 automaticamente é setado o bit0
   * */
  ShapeFunction(int degree=1, char flags=0) : _degree(degree), _flags(flags)
  {
    if (get_bit(_flags,1))
      set_bit(_flags,0);
    _ndof = (degree+1)*(degree+2)/2 + get_bit(_flags,0);
    this->_bar_nodes = genTriParametricPtsINT(degree);
    this->calculateDenomList();
  };

  void setDegree(int degree)
  {
    if (get_bit(_flags,1))
      set_bit(_flags,0);
    _ndof = (degree+1)*(degree+2)/2 + get_bit(_flags,0);
    this->_bar_nodes = genTriParametricPtsINT(degree);
    this->calculateDenomList();
		_degree = degree;
  }

  /** Defini o grau de função de lagrange.
   *  @param degree o grau
   *  @param flags os flags indicando se é uma função enriquecida, etc ... \n
   *           são 8 bits: \n
   *                 bit0 se com bolha \n
   *           bit1 se com bolha hierárquica
   *
   *  @note se setar o bit1 automaticamente é setado o bit0
   */
  void setDegree(int degree, char flags)
  {
    _flags = flags;
    this->setDegree(degree);
  }

  /*  calcula os "coeficientes (demominadores)" das funções de lagrange ... ver doc ..
   *  calcular a cada nova atualização de _degree;
   */
  /** NOT FOR USERS
   */
  void calculateDenomList()
  {
    this->_denom.resize(_bar_nodes.size());

    /* long int para evitar overflow */
    long int  denom_partF, denom_partG, denom_partH, Q0, Q1, Q2;

    // para cada ponto, calcular o denominador
    for (uint i = 0; i < _bar_nodes.size(); i++)
    {
      denom_partF = 1;
      denom_partG = 1;
      denom_partH = 1;

      Q1 = this->_bar_nodes[i][0];
      Q2 = this->_bar_nodes[i][1];
      Q0 = this->_degree-Q1-Q2;

      for (uint k = 0; k < Fepic::max(Q0, Q1, Q2); ++k)
      {
        if (k<Q0) denom_partF *= Q0 - k;
        if (k<Q1) denom_partG *= Q1 - k;
        if (k<Q2) denom_partH *= Q2 - k;
      }

      this->_denom[i] = denom_partF*denom_partG*denom_partH;

    }
  }

  /** Retorna a função avaliada nas coordenadas (L1, L2) do triângulo unitário
   */
  double operator() (double L1, double L2, int ith) const
  {
    if ((ith<0) || (ith>=_ndof))
    {
      std::cout << "error:ShapeFunction::operator(): invalid index." << std::endl;
      throw;
    }

    if (_degree==0)
      return 1.;

    /* função bolha */
    if (get_bit(_flags,0) && (ith==_ndof-1))
      return 27.*(1. - L1 - L2)*L1*L2;

    double L0 = 1. - L1 - L2;
    int Q1 = this->_bar_nodes[ith][0],
        Q2 = this->_bar_nodes[ith][1],
        Q0 = this->_degree-Q1-Q2;

    double result = _omega(_degree, Q0, L0)*
                    _omega(_degree, Q1, L1)*
                    _omega(_degree, Q2, L2);

    /* correção da bolha */
    if (get_bit(_flags,0) && (!get_bit(_flags,1)))
    {
      result -= pow(3,3-_degree)*_prod3(_degree, Q0)*
      _prod3(_degree, Q1)*
      _prod3(_degree, Q2)* (L0*L1*L2);
    }

    return result / _denom[ith];
  }

  /** Retorna a função avaliada nas coordenadas (L1, L2) do triângulo unitário
   */
  double operator() (Eigen::Vector2d const& L, int ith) const
  {
    return this->operator()(L[0], L[1], ith);
  }

  /**
   *  Retorna o gradiente no triangulo unitário, i.e.
   *
   */
  Eigen::Vector2d gradL(double L1, double L2, int ith) const
  {
    if ((ith<0) || (ith>=_ndof))
    {
      std::cout << "error:ShapeFunction::gradL(): invalid index." << std::endl;
      throw;
    }

    if (_degree == 0)
      return Eigen::Vector2d(0,0);

    double L0 = 1. - L1 - L2;
    int Q1 = this->_bar_nodes[ith][0],
        Q2 = this->_bar_nodes[ith][1],
        Q0 = this->_degree-Q1-Q2,
        n  = _degree;

    Eigen::Vector2d V;

    /* função bolha */
    if (get_bit(_flags,0) && (ith==_ndof-1))
    {
      V[0] = 27.*L2*(L0-L1);
      V[1] = 27.*L1*(L0-L2);
      return V;
    }

    /* lembrete: ddL1 = -ddL0 */
    V[0] = -_ddxomega(n, Q0, L0)*_omega(n, Q1, L1)*_omega(n, Q2, L2) +
            _ddxomega(n, Q1, L1)*_omega(n, Q2, L2)*_omega(n, Q0, L0);

    /* lembrete: ddL2 = -ddL0 */
    V[1] = -_ddxomega(n, Q0, L0)*_omega(n, Q1, L1)*_omega(n, Q2, L2) +
            _ddxomega(n, Q2, L2)*_omega(n, Q1, L1)*_omega(n, Q0, L0);

    /* correção da bolha */
    if (get_bit(_flags,0) && (!get_bit(_flags,1)))
    {
      V[0] -= pow(3,3-_degree)*_prod3(_degree, Q0)*
                               _prod3(_degree, Q1)*
                               _prod3(_degree, Q2)* L2*(L0-L1);

      V[1] -= pow(3,3-_degree)*_prod3(_degree, Q0)*
                               _prod3(_degree, Q1)*
                               _prod3(_degree, Q2)* L1*(L0-L2);
    }

    return V/_denom[ith];
  }

  Eigen::Vector2d gradL(Eigen::Vector2d const& L, int ith)
  {
    return this->gradL(L[0], L[1], ith);
  }

  /** Retorna o grau do polinômio.
   */
  int getDegree() const
  {
    return _degree;
  }

  /** Retorna o número de graus de liberdade.
   */
  int getNumDof() const
  {
    return (_degree+1)*(_degree+2)/2. + get_bit(_flags,0);
  }

  std::vector<Eigen::Vector2i> _bar_nodes; // integer baricentric nodes
  std::vector<long int>    _denom;   // ver doc
  int  _degree;
  char _flags;                             // _flags[0] : enriched_bubble;  _flags[1] : hierarc_bubble
  int  _ndof;
};

/*
 *                     (0,0,1)
 *                      /|`\
 *                     / |  `\
 *                    /  |    `\
 *                   /   |      `\
 *                  /    |        `\
 *                 /     |__________`\ (0,1,0)
 *                /     / (0,0,0)   ,'
 *               /    /         , '
 *              /   /       , '
 *             /  /     , '
 *            / /   , '
 *           // , '
 *          /,'
 *        (1,0,0)
 *
 */
/**
 *  Cria funções de interpolação de Lagrange no tetraedro unitário, podendo
 *  ser enriquecidos com função bolha (hierárquica ou não).
 */
template<>
class ShapeFunction<Simplex<3> >
{
public:
  /**
   *  @param degree o grau de intepolação
   *  @param flags 8 bits: \n
   *                 bit0 se com bolha \n
   *           bit1 se com bolha hierárquica
   *
   *  @note se setar o bit1 automaticamente é setado o bit0
   * */
  ShapeFunction(int degree=1, char flags=0) : _degree(degree), _flags(flags)
  {
    if (get_bit(_flags,1))
      set_bit(_flags,0);
    _ndof = (degree+1)*(degree+2)*(degree+3)/6 + get_bit(_flags,0);
    this->_bar_nodes = genTetParametricPtsINT(degree);
    this->calculateDenomList();
  };

  void setDegree(int degree)
  {
    if (get_bit(_flags,1))
      set_bit(_flags,0);
    _ndof = (degree+1)*(degree+2)*(degree+3)/6 + get_bit(_flags,0);
    this->_bar_nodes = genTetParametricPtsINT(degree);
    this->calculateDenomList();
        _degree = degree;
  }

  /** Defini o grau de função de lagrange.
   *  @param degree o grau
   *  @param flags os flags indicando se é uma função enriquecida, etc ... \n
   *           são 8 bits: \n
   *                 bit0 se com bolha \n
   *           bit1 se com bolha hierárquica
   *
   *  @note se setar o bit1 automaticamente é setado o bit0
   */
  void setDegree(int degree, char flags)
  {
    _flags = flags;
    this->setDegree(degree);
  }

  /*  calcula os "coeficientes" das funções de lagrange ... ver doc ..
   *  calcular a cada nova atualização de _degree;
   */
  /** NOT FOR USERS
   */
  void calculateDenomList()
  {
    this->_denom.resize(_bar_nodes.size());

    long int  denom_partF, denom_partG, denom_partH, denom_partI, Q0, Q1, Q2, Q3;

    // para cada ponto, calcular o denominador
    for (uint i = 0; i < _bar_nodes.size(); i++)
    {
      denom_partF = 1;
      denom_partG = 1;
      denom_partH = 1;
      denom_partI = 1;

      Q1 = this->_bar_nodes[i][0];
      Q2 = this->_bar_nodes[i][1];
      Q3 = this->_bar_nodes[i][2];
      Q0 = this->_degree-Q1-Q2-Q3;

      for (uint k = 0; k < Fepic::max(Q0,Q1,Q2,Q3); ++k)
      {
        if (k<Q0) denom_partF *= Q0 - k;
        if (k<Q1) denom_partG *= Q1 - k;
        if (k<Q2) denom_partH *= Q2 - k;
        if (k<Q3) denom_partI *= Q3 - k;
      }

      this->_denom[i] = denom_partF*denom_partG*denom_partH*denom_partI;

    }
  }

  /** Retorna a função avaliada nas coordenadas (L1, L2) do triângulo unitário
   */
  double operator() (double L1, double L2, double L3, int ith) const
  {
    if ((ith<0) || (ith>=_ndof))
    {
      std::cout << "error:ShapeFunction::operator(): invalid index." << std::endl;
      throw;
    }

    if (_degree==0)
    {
      return 1.;
    }

    /* função bolha */
    if (get_bit(_flags,0) && (ith==_ndof-1))
    {
      return 256.*(1.-L1-L2-L3)*L1*L2*L3;
    }

    double L0 = 1.-L1-L2-L3;
    int Q1 = this->_bar_nodes[ith][0],
        Q2 = this->_bar_nodes[ith][1],
        Q3 = this->_bar_nodes[ith][2],
        Q0 = this->_degree-Q1-Q2-Q3;

    double result = _omega(_degree, Q0, L0)*
                    _omega(_degree, Q1, L1)*
                    _omega(_degree, Q2, L2)*
                    _omega(_degree, Q3, L3);

    /* correção da bolha */
    if (get_bit(_flags,0) && (!get_bit(_flags,1)))
    {
      result -= pow(4,4-_degree)*_prod4(_degree, Q0)*
                                 _prod4(_degree, Q1)*
                                 _prod4(_degree, Q2)*
                                 _prod4(_degree, Q3)* (L0*L1*L2*L3);
    }

    return result / _denom[ith];
  }

  /** Retorna a função avaliada nas coordenadas (L1, L2) do triângulo unitário
   */
  double operator() (Eigen::Vector3d const& L, int ith) const
  {
    return this->operator()(L[0], L[1], L[2], ith);
  }

  /**
   *  Retorna o gradiente no triangulo unitário, i.e.
   *
   */
  Eigen::Vector3d gradL(double L1, double L2, double L3, int ith) const
  {
    if ((ith<0) || (ith>=_ndof))
    {
      std::cout << "error:ShapeFunction::gradL(): invalid index." << std::endl;
      throw;
    }

    if (_degree == 0)
      return Eigen::Vector3d(0,0,0);

    double L0 = 1.-L1-L2-L3;
    int Q1 = this->_bar_nodes[ith][0],
        Q2 = this->_bar_nodes[ith][1],
        Q3 = this->_bar_nodes[ith][2],
        Q0 = this->_degree-Q1-Q2-Q3,
        n  = _degree;

    Eigen::Vector3d V;

    /* função bolha */
    if (get_bit(_flags,0) && (ith==_ndof-1))
    {
      V[0] = 256.*L2*L3*(L0-L1);
      V[1] = 256.*L1*L3*(L0-L2);
      V[2] = 256.*L1*L2*(L0-L3);
      return V;
    }

    /* lembrete: ddL1 = -ddL0 */
    V[0] = -_ddxomega(n, Q0, L0)*_omega(n, Q1, L1)*_omega(n, Q2, L2)*_omega(n, Q3, L3) +
            _ddxomega(n, Q1, L1)*_omega(n, Q0, L0)*_omega(n, Q2, L2)*_omega(n, Q3, L3);

    /* lembrete: ddL2 = -ddL0 */
    V[1] = -_ddxomega(n, Q0, L0)*_omega(n, Q1, L1)*_omega(n, Q2, L2)*_omega(n, Q3, L3) +
            _ddxomega(n, Q2, L2)*_omega(n, Q1, L1)*_omega(n, Q0, L0)*_omega(n, Q3, L3);

    /* lembrete: ddL3 = -ddL0 */
    V[2] = -_ddxomega(n, Q0, L0)*_omega(n, Q1, L1)*_omega(n, Q2, L2)*_omega(n, Q3, L3) +
            _ddxomega(n, Q3, L3)*_omega(n, Q1, L1)*_omega(n, Q2, L2)*_omega(n, Q0, L0);

    /* correção da bolha */
    if (get_bit(_flags,0) && (!get_bit(_flags,1)))
    {
      V[0] -= pow(4,4-_degree)*_prod4(_degree, Q0)*
                               _prod4(_degree, Q1)*
                               _prod4(_degree, Q2)*
                               _prod4(_degree, Q3)* L2*L3*(L0-L1);

      V[1] -= pow(4,4-_degree)*_prod4(_degree, Q0)*
                               _prod4(_degree, Q1)*
                               _prod4(_degree, Q2)*
                               _prod4(_degree, Q3)* L1*L3*(L0-L2);

      V[2] -= pow(4,4-_degree)*_prod4(_degree, Q0)*
                               _prod4(_degree, Q1)*
                               _prod4(_degree, Q2)*
                               _prod4(_degree, Q3)* L1*L2*(L0-L3);
    }

    return V/_denom[ith];
  }

  Eigen::Vector3d gradL(Eigen::Vector3d const& L, int ith)
  {
    return this->gradL(L[0], L[1], L[2], ith);
  }

  /** Retorna o grau do polinômio.
   */
  int getDegree() const
  {
    return _degree;
  }

  /** Retorna o número de graus de liberdade.
   */
  int getNumDof() const
  {
    return (_degree+1)*(_degree+2)*(_degree+3)/6 + get_bit(_flags,0);
  }

  std::vector<Eigen::Vector3i> _bar_nodes; // integer baricentric nodes
  std::vector<long int>    _denom;   // ver doc
  int  _degree;
  char _flags;                             // _flags[0] : enriched_bubble;  _flags[1] : hierarc_bubble
  int  _ndof;
};

double gradL(double const& L, ShapeFunction<Simplex<1> > PHI, int ith)
{
  return PHI.gradL(L,ith);
}

Eigen::Vector2d gradL(Eigen::Vector2d const& L, ShapeFunction<Simplex<2> > PHI, int ith)
{
  return PHI.gradL(L,ith);
}

Eigen::Vector3d gradL(Eigen::Vector3d const& L, ShapeFunction<Simplex<3> > PHI, int ith)
{
  return PHI.gradL(L,ith);
}








#endif
