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

/* linear Triangle */
class Triangle3 :  public _CellCore<Triangle3>
{
public:

  friend class _CellCore<Triangle3>;

  typedef Simplex<2>  PolytopeT;
  typedef Edge2       Derived;

  enum { dim=2,
         n_vertices=3,
         n_nodes=3,
         n_facets=3,
         n_corners=3,
         n_vertices_per_facet=2,
         n_vertices_per_corner=1,
         n_nodes_per_facet=2,
         n_nodes_per_corner=1};

  static const bool has_edge_nodes = false;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = TRIANGLE3;
  static const EMshTag msh_tag = MSH_TRI_3;
  
  
  virtual ~Triangle3() {}

protected:
  char _icells_pos[3]; // positions of icells
  int  _nodes[3];      // nodes id
  int  _facets[3];     // facets id
  int  _icells[3];     // incident cells id
  int  _corners[3];

  /* não utilizados */
  int _icells_anchors[0];

public:
  
  static const int table_fC_x_vC[3][2];
  static const int table_fC_x_nC[3][2];
  static const int table_vC_x_fC[3][2];
  static const int table_fC_x_bC[3][2];
  static const int table_bC_x_vC[3][1];
  static const int table_bC_x_nC[3][1];
  static const int table_bC_x_fC[3][2];
  
};



/* quadratic Triangle */
class Triangle6 :  public _CellCore<Triangle6>
{
public:

  friend class _CellCore<Triangle6>;

  typedef  Simplex<2> PolytopeT;
  typedef  Edge3      Derived;

  enum { dim=2,
         n_vertices=3,
         n_nodes=6,
         n_facets=3,
         n_corners=3,
         n_vertices_per_facet=2,
         n_vertices_per_corner=1,
         n_nodes_per_facet=3,
         n_nodes_per_corner=1};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = TRIANGLE6;
  static const EMshTag msh_tag = MSH_TRI_6;
  
  virtual ~Triangle6() {};
  
protected:
  char _icells_pos[3]; // positions of icells
  int  _facets[3];     // facets id
  int  _icells[3];     // incident cells id
  int  _nodes[6];      // nodes id
  int  _corners[3];

  /* não utilizados */
  int _icells_anchors[0];
  

public:

  static const int table_fC_x_vC[3][2];
  static const int table_fC_x_nC[3][3];
  static const int table_vC_x_fC[3][2];
  static const int table_fC_x_bC[3][2];
  static const int table_bC_x_vC[3][1];
  static const int table_bC_x_nC[3][1];
  static const int table_bC_x_fC[3][2];

};




#endif



