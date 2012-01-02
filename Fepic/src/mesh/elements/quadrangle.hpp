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

#ifndef FEPIC_QUADRANGLE_HPP
#define FEPIC_QUADRANGLE_HPP

/* linear Quadrangle */
class Quadrangle4 :  public _CellCore<Quadrangle4>
{
public:
  friend class _CellCore<Quadrangle4>;
  
  typedef Hypercube<2> PolytopeT;
  typedef Edge2        Derived;

  enum { dim=2,
         n_vertices=4,
         n_nodes=4,
         n_facets=4,
         n_corners=4,
         n_vertices_per_facet=2,
         n_vertices_per_corner=1,
         n_nodes_per_facet=2,
         n_nodes_per_corner=1};

  static const bool has_edge_nodes = false;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = QUADRANGLE4;
  static const EMshTag msh_tag = MSH_QUA_4;

  virtual ~Quadrangle4() {}

protected:
  char _icells_pos[4];  // positions of icells
  int _facets[4];     // facets id
  int _icells[4];     // incident cells id
  int _nodes[4];      // nodes id

  /* não utilizados */
  int _icells_anchors[0];
  int _corners[0];



public:
  
  static const int table_fC_x_vC[4][2];
  static const int table_fC_x_nC[4][2];
  static const int table_vC_x_fC[4][2];
  static const int table_fC_x_bC[4][2];
  static const int table_bC_x_vC[4][1];
  static const int table_bC_x_nC[4][1];
  static const int table_bC_x_fC[4][2];
    
};



/* quadratic Quadrangle (serendipity) */
class Quadrangle8 :  public _CellCore<Quadrangle8>
{
public:
  friend class _CellCore<Quadrangle8>;

  typedef Hypercube<2> PolytopeT;
  typedef Edge3        Derived;

  enum { dim=2,
         n_vertices=4,
         n_nodes=8,
         n_facets=4,
         n_corners=4,
         n_vertices_per_facet=2,
         n_vertices_per_corner=1,
         n_nodes_per_facet=3,
         n_nodes_per_corner=1};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = QUADRANGLE8;
  static const EMshTag msh_tag = MSH_QUA_8;

  virtual ~Quadrangle8() {}

protected:
  char _icells_pos[4];  // positions of icells
  int _facets[4];     // facets id
  int _icells[4];     // incident cells id
  int _nodes[8];      // nodes id

  /* não utilizados */
  int _icells_anchors[0];
  int _corners[0];
  
public:
  
  static const int table_fC_x_vC[4][2];
  static const int table_fC_x_nC[4][3];
  static const int table_vC_x_fC[4][2];
  static const int table_fC_x_bC[4][2];
  static const int table_bC_x_vC[4][1];
  static const int table_bC_x_nC[4][1];
  static const int table_bC_x_fC[4][2];

};



/* quadratic Quadrangle */
class Quadrangle9 :  public _CellCore<Quadrangle9>
{
public:
  friend class _CellCore<Quadrangle9>;

  typedef Hypercube<2> PolytopeT;
  typedef Edge3        Derived;

  enum { dim=2,
         n_vertices=4,
         n_nodes=9,
         n_facets=4,
         n_corners=4,
         n_vertices_per_facet=2,
         n_vertices_per_corner=1,
         n_nodes_per_facet=3,
         n_nodes_per_corner=1};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = true;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = QUADRANGLE9;
  static const EMshTag msh_tag = MSH_QUA_9;

  virtual ~Quadrangle9() {}

protected:
  char _icells_pos[4];  // positions of icells
  int _facets[4];     // facets id
  int _icells[4];     // incident cells id
  int _nodes[9];      // nodes id

  /* não utilizados */
  int _icells_anchors[0];
  int _corners[0];
  
public:
  
  static const int table_fC_x_vC[4][2];
  static const int table_fC_x_nC[4][3];
  static const int table_vC_x_fC[4][2];
  static const int table_fC_x_bC[4][2];
  static const int table_bC_x_vC[4][1];
  static const int table_bC_x_nC[4][1];
  static const int table_bC_x_fC[4][2];
    
};



#endif



