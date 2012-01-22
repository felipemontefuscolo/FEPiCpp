#ifndef FEPIC_ENUMS_HPP
#define FEPIC_ENUMS_HPP


enum EMshTag {

  MSH_LIN_2  = 1 ,   // 2-node line.
  MSH_TRI_3  = 2 ,   // 3-node triangle.
  MSH_QUA_4  = 3 ,   // 4-node quadrangle.
  MSH_TET_4  = 4 ,   // 4-node tetrahedron.
  MSH_HEX_8  = 5 ,   // 8-node hexahedron.
  MSH_PRI_6  = 6 ,   // 6-node prism.
  MSH_PYR_5  = 7 ,   // 5-node pyramid.
  MSH_LIN_3  = 8 ,   // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
  MSH_TRI_6  = 9 ,   // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
  MSH_QUA_9  = 10,   // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
  MSH_TET_10 = 11,   // 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
  MSH_HEX_27 = 12,   // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
  MSH_PRI_18 = 13,   // 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
  MSH_PYR_14 = 14,   // 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
  MSH_PNT    = 15,   // 1-node point.
  MSH_QUA_8  = 16,   // 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
  MSH_HEX_20 = 17,   // 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
  MSH_PRI_15 = 18,   // 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
  MSH_PYR_13 = 19,   // 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
  MSH_TRI_9  = 20,   // 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
  MSH_TRI_10 = 21,   // 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
  MSH_TRI_12 = 22,   // 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
  MSH_TRI_15 = 23,   // 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
  MSH_TRI_15I= 24,   // 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
  MSH_TRI_21 = 25,   // 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
  MSH_LIN_4  = 26,   // 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
  MSH_LIN_5  = 27,   // 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
  MSH_LIN_6  = 28,   // 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
  MSH_TET_20 = 29,   // 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
  MSH_TET_35 = 30,   // 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
  MSH_TET_56 = 31,   // 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
  MSH_TET_34 = 32,
  MSH_UNDEFINED_ELEM = 33,
  MSH_MAX_ELEMENT_NODES = 56,
  MSH_MAX_INDEX= 32  // must be equal to the higher valid index
};

// ------------------------------------------------
// WARNING: enums must be power of two
// ------------------------------------------------
enum ECellType
{
  UNDEFINED_CELLT  = 0x0,
  POINT1           = 0x1,
  EDGE2            = 0x2,
  EDGE3            = 0x4,
  TRIANGLE3        = 0x8,
  TRIANGLE6        = 0x10,
  QUADRANGLE4      = 0x20,
  QUADRANGLE8      = 0x40,
  QUADRANGLE9      = 0x80,
  TETRAHEDRON4     = 0x100,
  TETRAHEDRON10    = 0x200,
  HEXAHEDRON8      = 0x400,
  HEXAHEDRON20     = 0x800,
  HEXAHEDRON27     = 0x1000,
  N_CELL_TYPES     = 13    // must be equal to the number of valid cells

};

enum ECellClass
{
  UNDEFINED_CELLC  = 0x0,
  POINT            = 0x1,
  EDGE             = 0x2,
  TRIANGLE         = 0x4,
  QUADRANGLE       = 0x8,
  TETRAHEDRON      = 0x10,
  HEXAHEDRON       = 0x20,
  N_CELL_CLASSES   = 6
  
};

enum ECellFamily
{
  UNDEFINED_CELLF  = 0x0,
  SIMPLEX          = 0x1,
  HCUBE            = 0x2,
  N_CELL_FAMILIES  = 2
  
};



enum ETag
{
  UNDEFINED_TAG  = 0x0
};

#endif


