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

#ifndef FEPIC_METAMESH_HPP
#define FEPIC_METAMESH_HPP



/*
  ==================================================================
  ========================== _MeshMethods ===========================
  ==================================================================
*/

template<class Mesh, int CellDim>
class _MeshMethods;

template<class Mesh>
class _MeshMethods<Mesh, 2>
{
public:
  static void buildAdjacency(Mesh & mesh)
  {
    mesh.buildAdjacency4face();
  }

  static void remodelCellsNodes(Mesh & mesh, int order)
  {
    mesh.remodelCellsNodes4face(order);
  }

};

template<class Mesh>
class _MeshMethods<Mesh, 3>
{
public:
  static void buildAdjacency(Mesh & mesh)
  {
    mesh.buildAdjacency4volume();
  }

  static void remodelCellsNodes(Mesh & mesh, int order)
  {
    mesh.remodelCellsNodes4volume(order);
  }

};



/*
  ==================================================================
  ================== _MeshReadMarkedElementsMsh =====================
  ==================================================================
*/

template<class Mesh, int CellDim>
class _MeshReadMarkedElementsMsh;

template<class Mesh>
class _MeshReadMarkedElementsMsh<Mesh, 1>
{
public:
  static void buildAdjacency(std::ifstream &File)
  {
    Mesh::readMarkedElementsMsh4edge(File);
  }
};

template<class Mesh>
class _MeshReadMarkedElementsMsh<Mesh, 2>
{
public:
  static void buildAdjacency(std::ifstream &File)
  {
    Mesh::readMarkedElementsMsh4face(File);
  }
};


template<class Mesh>
class _MeshReadMarkedElementsMsh<Mesh, 3>
{
public:
  static void buildAdjacency(std::ifstream &File)
  {
    Mesh::readMarkedElementsMsh4volume(File);
  }
};


#endif // METAMESH_HPP

