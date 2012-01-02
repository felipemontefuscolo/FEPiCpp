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

#ifndef FEPIC_HEXAHEDRON_HPP
#define FEPIC_HEXAHEDRON_HPP

template<class CT> class _CellCore;

/* linear hexahedron */
class Hexahedron8 :  public _CellCore<Hexahedron8>
{
public:
  friend class _CellCore<Hexahedron8>;

  typedef Hypercube<3> PolytopeT;
  typedef Quadrangle4  Derived;

  enum { dim=3,
         n_vertices=8,
         n_nodes=8,
         n_facets=6,
         n_corners=12,
         n_vertices_per_facet=4,
         n_vertices_per_corner=2,
         n_nodes_per_facet=4,
         n_nodes_per_corner=2};

  static const bool has_edge_nodes = false;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = HEXAHEDRON8;
  static const EMshTag msh_tag = MSH_HEX_8;

  
protected:
  char _icells_pos[6];     // positions of icells
  char _icells_anchors[6]; // anchors of icells  
  int _facets[6];          // facets id
  int _icells[6];          // incident cells id
  int _nodes[8];           // nodes id
  int _corners[12];        // edges id

public:
  static const int table_vC_x_fC[8][6];
  static const int table_fC_x_bC[6][4];
  static const int table_fC_x_vC[6][8];
  static const int table_fC_x_nC[6][4];
  static const int table_bC_x_vC[12][2];
  static const int table_bC_x_nC[12][2];
  static const int table_bC_x_fC[12][4];
    
};



/* quadratic hexahedron (serendipity) */
class Hexahedron20 :  public _CellCore<Hexahedron20>
{
public:
  friend class _CellCore<Hexahedron20>;

  typedef Hypercube<3> PolytopeT;
  typedef Quadrangle8  Derived;

  enum { dim=3,
         n_vertices=8,
         n_nodes=20,
         n_facets=6,
         n_corners=12,
         n_vertices_per_facet=4,
         n_vertices_per_corner=2,
         n_nodes_per_facet=8,
         n_nodes_per_corner=3};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = HEXAHEDRON20;
  static const EMshTag msh_tag = MSH_HEX_20;

protected:
  char _icells_pos[6];     // positions of icells
  char _icells_anchors[6]; // anchors of icells  
  int _facets[6];          // facets id
  int _icells[6];          // incident cells id
  int _nodes[20];           // nodes id
  int _corners[12];        // edges id

public:
  static const int table_vC_x_fC[8][6];
  static const int table_fC_x_bC[6][4];
  static const int table_fC_x_vC[6][8];
  static const int table_fC_x_nC[6][8];
  static const int table_bC_x_vC[12][2];
  static const int table_bC_x_nC[12][3];
  static const int table_bC_x_fC[12][4];
    
};



/* quadratic hexahedron (serendipity) */
class Hexahedron27 :  public _CellCore<Hexahedron27>
{
public:
  friend class _CellCore<Hexahedron27>;

  typedef Hypercube<3> PolytopeT;
  typedef Quadrangle9  Derived;

  enum { dim=3,
         n_vertices=8,
         n_nodes=27,
         n_facets=6,
         n_corners=12,
         n_vertices_per_facet=4,
         n_vertices_per_corner=2,
         n_nodes_per_facet=9,
         n_nodes_per_corner=3};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = true;
  static const bool has_volume_nodes = true;
  static const ECellType fep_tag = HEXAHEDRON27;
  static const EMshTag msh_tag = MSH_HEX_27;

protected:
  char _icells_pos[6];     // positions of icells
  char _icells_anchors[6]; // anchors of icells  
  int _facets[6];          // facets id
  int _icells[6];          // incident cells id
  int _nodes[27];           // nodes id
  int _corners[12];        // edges id

public:
  static const int table_vC_x_fC[8][6];
  static const int table_fC_x_bC[6][4];
  static const int table_fC_x_vC[6][8];
  static const int table_fC_x_nC[6][9];
  static const int table_bC_x_vC[12][2];
  static const int table_bC_x_nC[12][3];
  static const int table_bC_x_fC[12][4];
    
};



#endif
