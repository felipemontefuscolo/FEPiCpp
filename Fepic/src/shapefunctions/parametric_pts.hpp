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

#ifndef PARAMETRIC_PTS_HPP
#define PARAMETRIC_PTS_HPP

#include <iostream>
#include "Fepic/src/custom_eigen/custom_eigen.hpp"
#include <vector>
#include <algorithm>

/* protótipos */
std::vector<double> genLineParametricPts(int n);

std::vector<Eigen::Vector2d> genTriParametricPts(int n);
std::vector<Eigen::Vector2i> genTriParametricPtsINT(int n, int a=0, int b=-1);

std::vector<Eigen::Vector3d> genTetParametricPts(int n);
std::vector<Eigen::Vector3i> genTetParametricPtsINT(int n, int a=0, int b=-1);





/*
     0----------1 --> u      0-----2----1     0----2----3----1
   |          |
   -1         +1

*/
/** Retorna um vetor com as coordenadas paramétricas dos pontos de uma
 *  reta de ordem n.
 */
std::vector<double> genLineParametricPts(int n)
{
  std::vector<double> list;
  if(n==0)
  {
    list.push_back(0);
    return list;
  }
  if (n >= 1)
  {
    list.push_back(-1);
    list.push_back(+1);

    for (int i = 1; i <= n-1; i++)
    {
      list.push_back(-1 + i*2./n);
    }

    return list;
  }
  return list;
}

/*
     v
     ^
     |
     2 (0,1)
     |`\
     |  `\
     |    `\
     |      `\
     |        `\
     0----------1 --> u
   (0,0)       (1,0)

*/
/** Retorna um vetor com as coordenadas paramétricas dos pontos de um
 *  triângulo de ordem n.
 */
std::vector<Eigen::Vector2d> genTriParametricPts(int n)
{
  std::vector<Eigen::Vector2d> list;

  if(n==0)
  {
    list.push_back(Eigen::Vector2d(1./3., 1./3.));
    return list;
  }
  if(n>0)
  {
    /* é usado o gerador de pontos de coordenadas inteiras para reduzir
     * o erro */
    std::vector<Eigen::Vector2i> temp(genTriParametricPtsINT(n));
    std::vector<Eigen::Vector2i>::iterator tit = temp.begin();

    for (; tit != temp.end() ; ++tit)
    {
      list.push_back((*tit).cast<double>()/n);
    }
    return list;
  }

  return list; //vazio
}

/* Retorna um vetor com as coordenadas (em inteiros) dos pontos de um
 *  triângulo de ordem n.
 *
 * As coordenadas são retornadas em inteiros para a função retornar a coordenada
 * exata.
 * Dessa forma, se o triângulo é de ordem n, seu lado tem comprimento n.
 *
 * @warning APENAS USAR COM O PRIMEIRO ARGUMENTO ... os outros argumentos
 * são usados para funcionalidades internas.
 */
/** NOT FOR USERS
 */
std::vector<Eigen::Vector2i> genTriParametricPtsINT(int n, int a, int b)
{
  if (b==-1)
  {
    b=n;
  }

  std::vector<Eigen::Vector2i> list;

  int l = b-a;

  if(n==0)
  {
    list.push_back(Eigen::Vector2i(a+l/3, a+l/3));
    return list;
  }
  if(n>0)
  {
    list.push_back(Eigen::Vector2i(a, a));
    list.push_back(Eigen::Vector2i(b, a));
    list.push_back(Eigen::Vector2i(a, b));

    // aresta 0 (0-1)
    for (int i = 1; i <= n-1; ++i)
    {
      list.push_back(Eigen::Vector2i(a + i, a));
    }

    // aresta 1 (1-2)
    for (int i = n-1; i > 0; --i)
    {
      list.push_back(Eigen::Vector2i(a + i, b - i));
    }

    // aresta 2 (2-0)
    for (int i = 1; i <= n-1; ++i)
    {
      list.push_back(Eigen::Vector2i(a, b - i));
    }

    // pontos interiores
    if (n>2)
    {
      std::vector<Eigen::Vector2i> temp;
      int aa = a + 1;
      int bb = b - 2;
      temp = genTriParametricPtsINT(n-3, aa, bb);
      list.insert(list.end(), temp.begin(), temp.end());
      return list;
    }
    else
    {
      return list;
    }

  }

  return list; // vazio
}

/*
                    v
                      .
                    ,/
                   /
                2
              ,/|`\
            ,/  |  `\
          ,/    '.   `\
        ,/       |     `\
      ,/         |       `\
     0-----------'.--------1 --> u
      `\.         |      ,/
         `\.      |    ,/
            `\.   '. ,/
               `\. |/
                  `3
                     `\.
                        ` w
*/
/** Retorna um vetor com as coordenadas paramétricas dos pontos de um
 *  triângulo de ordem n.
 */
std::vector<Eigen::Vector3d> genTetParametricPts(int n)
{
  std::vector<Eigen::Vector3d> list;

  if(n==0)
  {
    list.push_back(Eigen::Vector3d(1./4., 1./4., 1./4.));
    return list;
  }
  if(n>0)
  {
    /* é usado o gerador de pontos de coordenadas inteiras para reduzir
     * o erro */
    std::vector<Eigen::Vector3i> temp(genTetParametricPtsINT(n));
    std::vector<Eigen::Vector3i>::iterator tit = temp.begin();

    for (; tit != temp.end() ; ++tit)
    {
      list.push_back((*tit).cast<double>()/n);
    }
    return list;
  }

  return list; //vazio
}

/* Retorna um vetor com as coordenadas (em inteiros) dos pontos de um
 *  tetraedro de ordem n.
 *
 * As coordenadas são retornadas em inteiros para a função retornar a coordenada
 * exata.
 * Dessa forma, se o tetraedro é de ordem n, seu lado tem comprimento n.
 *
 * @warning APENAS USAR COM O PRIMEIRO ARGUMENTO ... os outros argumentos
 * são usados para funcionalidades internas.
 */
/** NOT FOR USERS
 */
std::vector<Eigen::Vector3i> genTetParametricPtsINT(int n, int a, int b)
{
  if (b==-1)
    b=n;

  std::vector<Eigen::Vector3i> list;

  int l = b-a;

  if(n==0)
  {
    list.push_back(Eigen::Vector3i(a+l/4, a+l/4, a+l/4));
    return list;
  }
  if(n>0)
  {
    list.push_back(Eigen::Vector3i(a, a, a));
    list.push_back(Eigen::Vector3i(b, a, a));
    list.push_back(Eigen::Vector3i(a, b, a));
    list.push_back(Eigen::Vector3i(a, a, b));

    // aresta 0 (0-1)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(a + i, a, a));

    // aresta 1 (1-2)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(b - i, a + i, a));

    // aresta 2 (2-0)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(a, b - i, a));

    // aresta 3 (3-0)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(a, a, b-i));

    // aresta 4 (3-2)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(a, a+i, b-i));

    // aresta 5 (3-1)
    for (int i = 1; i <= n-1; ++i)
      list.push_back(Eigen::Vector3i(a+i, a, b-i));

    // faces
    if(n>2)
    {
      std::vector<Eigen::Vector2i> face(genTriParametricPtsINT(n-3));
      Eigen::MatrixXi M(3,2);
      Eigen::Vector3i vi, vj, vk; // vértices da face reduzida {i,j,k}

      // face 0 {1,0,2}
      vi << b-2, a+1, a;
      vj << a+1, a+1, a;
      vk << a+1, b-2, a;
      M << (vi-vk),(vj-vk);
      for (uint i = 0; i < face.size(); i++)
      {
        int N = n==3 ? 1 : (n-3);
        list.push_back(M*face[i]/N + vk);
        /* M*face[i]/(n-3) + vk : é uma transformação.
           Divide-se por N p/ normalização de face[i] */
      }

      // face 1 {0,1,3}
      vi << a+1, a, a+1;
      vj << b-2, a, a+1;
      vk << a+1, a, b-2;
      M << (vi-vk),(vj-vk);
      for (uint i = 0; i < face.size(); i++)
      {
        int N = n==3 ? 1 : (n-3);
        list.push_back(M*face[i]/N + vk);
      }

      // face 2 {3,2,0}
      vi << a, a+1, b-2;
      vj << a, a+1, a+1;
      vk << a, b-2, a+1;
      M << (vi-vk),(vj-vk);
      for (uint i = 0; i < face.size(); i++)
      {
        int N = n==3 ? 1 : (n-3);
        list.push_back(M*face[i]/N + vk);
      }

      // face 3 {2,3,1}
      vi << a+1, b-2, a+1;
      vj << a+1, a+1, b-2;
      vk << b-2, a+1, a+1;
      M << (vi-vk),(vj-vk);
      for (uint i = 0; i < face.size(); i++)
      {
        int N = n==3 ? 1 : (n-3);
        list.push_back(M*face[i]/N + vk);
      }
    }

    // pontos interiores
    if (n>3)
    {
      std::vector<Eigen::Vector3i> temp;
      int aa = a + 1;
      int bb = b - 3;
      temp = genTetParametricPtsINT(n-4, aa, bb);
      list.insert(list.end(), temp.begin(), temp.end());
      return list;
    }
    else
      return list;

  }

  return list;
}


#endif











