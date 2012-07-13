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

/* linear tetrahedron */
class Tetrahedron4 :  public _CellCore<Tetrahedron4>
{
public:

  friend class _CellCore<Tetrahedron4>;

  typedef Simplex<3> PolytopeT;
  typedef Triangle3  Derived;

  enum { dim=3,
         n_vertices=4,
         n_nodes=4,
         n_facets=4,
         n_corners=6,
         n_vertices_per_facet=3,
         n_vertices_per_corner=2,
         n_nodes_per_facet=3,
         n_nodes_per_corner=2};

  static const bool has_edge_nodes = false;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = TETRAHEDRON4;
  static const EMshTag msh_tag = MSH_TET_4;

  virtual ~Tetrahedron4() {}
  
protected:
  char _icells_pos[4];      // positions of icells
  char _icells_anchors[4];  // anchors of icells  
  int  _nodes[4];      // nodes id
  int  _facets[4];     // facets id
  int  _icells[4];     // incident cells id
  int  _corners[6];    // edges id

public:
  static const int table_vC_x_fC[4][6];
  static const int table_fC_x_bC[4][3];
  static const int table_fC_x_vC[4][6];
  static const int table_fC_x_nC[4][3];
  static const int table_bC_x_vC[6][2];
  static const int table_bC_x_nC[6][2];
  static const int table_bC_x_fC[6][4];
    
};


/* quadratic tetrahedron */
class Tetrahedron10 : public _CellCore<Tetrahedron10>
{
public:

  friend class _CellCore<Tetrahedron10>;

  typedef Simplex<3> PolytopeT;
  typedef Triangle6  Derived;

  enum { dim=3,
         n_vertices=4,
         n_nodes=10,
         n_facets=4,
         n_corners=6,
         n_vertices_per_facet=3,
         n_vertices_per_corner=2,
         n_nodes_per_facet=6,
         n_nodes_per_corner=3};

  static const bool has_edge_nodes = true;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = TETRAHEDRON10;
  static const EMshTag msh_tag = MSH_TET_10;

  virtual ~Tetrahedron10() {}

protected:
  char _icells_pos[4];      // positions of icells
  char _icells_anchors[4];  // anchors of icells  
  int  _facets[4];     // facets id
  int  _icells[4];     // incident cells id
  int  _corners[6];    // edges id
  int  _nodes[10];      // nodes id
  int  _conn_comp_id;  // connected component to which the cell belongs.

public:
  static const int table_vC_x_fC[4][6];
  static const int table_fC_x_bC[4][3];
  static const int table_fC_x_vC[4][6];
  static const int table_fC_x_nC[4][6];
  static const int table_bC_x_vC[6][2];
  static const int table_bC_x_nC[6][3];
  static const int table_bC_x_fC[6][4];
    
};

// only change table_fC_x_nC and table_bC_x_nC compared to Tetrahedron4




#endif
