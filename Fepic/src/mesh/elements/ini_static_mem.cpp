#include "cellcore.hpp"
#include "elements.hpp"

/*       NOTATION
 *  V  = global vertex id
 *  C  = global cell id
 *  N  = global node id
 *  F  = global facet id
 *  B  = global corner id
 * 
 *  v or vC  = vertex id within the cell C
 *  n or nC  = node id within the cell C
 *  f or fC  = facet id within the cell C
 *  b or bC  = corner id within the cell C
 * 
 *  vf = v-th vertex of the facet f (within a cell)
 *  nf = n-th node of the facet f (within a cell)
 *  bf = b-th corner of the facet f (within a cell)
 *  
 *  fv = f-th incident facet of the vertex v (within a cell)
 *  fn = f-th incident facet of the node v (within a cell)
 * 
 *  and so on ...
 */ 


// ========================  EDGE =====================================

#define FEPIC_INSTANTIATE_MEMBERS_FUN(Order)                                \
template<> const int Edge<Order>::table_fC_x_vC[2][1] = {{0},{1}};          \
template<> const int Edge<Order>::table_fC_x_nC[2][1] = {{0},{1}};          \
                                                                            \
template<> const int Edge<Order>::table_bC_x_vC[0][0] = {};                 \
template<> const int Edge<Order>::table_bC_x_nC[0][0] = {};                 \
template<> const int Edge<Order>::table_fC_x_bC[0][0] = {};                 \
template<> const int Edge<Order>::table_bC_x_fC[0][0] = {};                 \
template<> const int Edge<Order>::table_vC_x_fC[2][1] = {{0},{1}};

FEPIC_INSTANTIATE_MEMBERS_FUN(1)
FEPIC_INSTANTIATE_MEMBERS_FUN(2)

#undef FEPIC_INSTANTIATE_MEMBERS_FUN

// ========================  TRIANGLE =====================================


const int Triangle3::table_fC_x_vC[3][2] = {{0,1}, {1,2}, {2,0}};
const int Triangle3::table_fC_x_nC[3][2] = {{0,1}, {1,2}, {2,0}};
const int Triangle3::table_vC_x_fC[3][2] = {{0,2}, {1,0}, {2,1}};
const int Triangle3::table_fC_x_bC[3][2] = {{0,1}, {1,2}, {2,0}};
const int Triangle3::table_bC_x_vC[3][1] = {{0},{1},{2}};
const int Triangle3::table_bC_x_nC[3][1] = {{0},{1},{2}};
const int Triangle3::table_bC_x_fC[3][2] = {{0,2}, {1,0}, {2,1}};



const int Triangle6::table_fC_x_vC[3][2] = {{0,1}, {1,2}, {2,0}};
const int Triangle6::table_fC_x_nC[3][3] = {{0,1,3}, {1,2,4}, {2,0,5}};
const int Triangle6::table_vC_x_fC[3][2] = {{0,2}, {1,0}, {2,1}};
const int Triangle6::table_fC_x_bC[3][2] = {{0,1}, {1,2}, {2,0}};
const int Triangle6::table_bC_x_vC[3][1] = {{0},{1},{2}};
const int Triangle6::table_bC_x_nC[3][1] = {{0},{1},{2}};
const int Triangle6::table_bC_x_fC[3][2] = {{0,2}, {1,0}, {2,1}};



// ========================  QUADRANGLE =====================================


const int Quadrangle4::table_fC_x_vC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle4::table_fC_x_nC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle4::table_vC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};
const int Quadrangle4::table_fC_x_bC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle4::table_bC_x_vC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle4::table_bC_x_nC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle4::table_bC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};


const int Quadrangle8::table_fC_x_vC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle8::table_fC_x_nC[4][3] = {{0,1,4}, {1,2,5}, {2,3,6}, {3,0,7}};
const int Quadrangle8::table_vC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};
const int Quadrangle8::table_fC_x_bC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle8::table_bC_x_vC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle8::table_bC_x_nC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle8::table_bC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};


const int Quadrangle9::table_fC_x_vC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle9::table_fC_x_nC[4][3] = {{0,1,4}, {1,2,5}, {2,3,6}, {3,0,7}};
const int Quadrangle9::table_vC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};
const int Quadrangle9::table_fC_x_bC[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
const int Quadrangle9::table_bC_x_vC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle9::table_bC_x_nC[4][1] = {{0}, {1}, {2}, {3}};
const int Quadrangle9::table_bC_x_fC[4][2] = {{0,3}, {1,0}, {2,1}, {3,2}};


// ========================  TETRAHEDRON =====================================

                                                  // fC        vf
const int Tetrahedron4::table_vC_x_fC[4][6]  = {{0, 1, 2,   1, 0, 2},
                                                {0, 1, 3,   0, 1, 2},
                                                {0, 2, 3,   2, 1, 0},
                                                {1, 2, 3,   2, 0, 1}};
const int Tetrahedron4::table_fC_x_bC[4][3]   = {{0, 2, 1},
                                                 {0, 5, 3},
                                                 {4, 2, 3},
                                                 {4, 5, 1}};
                                                 
                                                 // vC   //fv
const int Tetrahedron4::table_fC_x_vC[4][6]  = {{1,0,2,   0,0,0},
                                                {0,1,3,   1,1,0},
                                                {3,2,0,   1,1,2},
                                                {2,3,1,   2,2,2}};
// same as above
const int Tetrahedron4::table_fC_x_nC[4][3] = {{1,0,2},{0,1,3},{3,2,0},{2,3,1}};

const int Tetrahedron4::table_bC_x_vC[6][2] = {{0,1},{1,2},{2,0},{3,0},{3,2},{3,1}};

const int Tetrahedron4::table_bC_x_nC[6][2] = {{0,1},{1,2},{2,0},{3,0},{3,2},{3,1}};

const int Tetrahedron4::table_bC_x_fC[6][4] = {{0, 1, 0, 0},
                                               {0, 3, 2, 2},
                                               {0, 2, 1, 1},
                                               {2, 1, 2, 2},
                                               {3, 2, 0, 0},
                                               {1, 3, 1, 1}};

const int Tetrahedron10::table_vC_x_fC[4][6]  = {{0, 1, 2, 1, 0, 2},
                                                 {0, 1, 3, 0, 1, 2},
                                                 {0, 2, 3, 2, 1, 0},
                                                 {1, 2, 3, 2, 0, 1}};
const int Tetrahedron10::table_fC_x_bC[4][3]   = {{0, 2, 1},
                                                  {0, 5, 3},
                                                  {4, 2, 3},
                                                  {4, 5, 1}};
                                                  
                                                 // vC   //fv
const int Tetrahedron10::table_fC_x_vC[4][6]  = {{1,0,2,   0,0,0},
                                                 {0,1,3,   1,1,0},
                                                 {3,2,0,   1,1,2},
                                                 {2,3,1,   2,2,2}};

const int Tetrahedron10::table_fC_x_nC[4][6]     = {{1,0,2,4,6,5},
                                                    {0,1,3,4,9,7},
                                                    {3,2,0,8,6,7},
                                                    {2,3,1,8,9,5}};

const int Tetrahedron10::table_bC_x_vC[6][2] = {{0,1},{1,2},{2,0},{3,0},{3,2},{3,1}};

const int Tetrahedron10::table_bC_x_nC[6][3]    = {{0,1,4},
                                                   {1,2,5},
                                                   {2,0,6},
                                                   {3,0,7},
                                                   {3,2,8},
                                                   {3,1,9}};

const int Tetrahedron10::table_bC_x_fC[6][4]   = {{0, 1, 0, 0},
                                                  {0, 3, 2, 2},
                                                  {0, 2, 1, 1},
                                                  {2, 1, 2, 2},
                                                  {3, 2, 0, 0},
                                                  {1, 3, 1, 1}};


// ========================  HEXAHEDRON =====================================



const int Hexahedron8::table_vC_x_fC[8][6]  = {{0, 1, 2, 0, 0, 0},
                                               {0, 1, 3, 3, 1, 0},
                                               {0, 3, 4, 2, 1, 0},
                                               {0, 2, 4, 1, 3, 1},
                                               {1, 2, 5, 3, 1, 0},
                                               {1, 3, 5, 2, 3, 1},
                                               {3, 4, 5, 2, 3, 2},
                                               {2, 4, 5, 2, 2, 3}};
const int Hexahedron8::table_fC_x_bC[6][4]   = {{1, 5 , 3 , 0},
                                                {0, 4 , 8 , 2},
                                                {2, 9 , 7 , 1},
                                                {3, 6 , 10, 4},
                                                {5, 7 , 11, 6},
                                                {8, 10, 11, 9}};
                                                
                                                // vC         fv
const int Hexahedron8::table_fC_x_vC[6][8]  = {{0,3,2,1,   0,0,0,0},
                                               {0,1,5,4,   1,1,0,0},
                                               {0,4,7,3,   2,1,0,1},
                                               {1,2,6,5,   2,1,0,1},
                                               {2,3,7,6,   2,2,1,1},
                                               {4,5,6,7,   2,2,2,2}};
// same as above
const int Hexahedron8::table_fC_x_nC[6][4]     = {{0,3,2,1},{0,1,5,4},{0,4,7,3},{1,2,6,5},{2,3,7,6},{4,5,6,7}};

const int Hexahedron8::table_bC_x_vC[12][2]= {{0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},{5,6},{6,7}};

const int Hexahedron8::table_bC_x_nC[12][2]   = {{0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},{5,6},{6,7}};

const int Hexahedron8::table_bC_x_fC[12][4]  = {{0, 1, 3, 0},
                                                {2, 0, 3, 0},
                                                {1, 2, 3, 0},
                                                {0, 3, 2, 0},
                                                {3, 1, 3, 1},
                                                {0, 4, 1, 0},
                                                {4, 3, 3, 1},
                                                {2, 4, 2, 1},
                                                {1, 5, 2, 0},
                                                {5, 2, 3, 1},
                                                {3, 5, 2, 1},
                                                {4, 5, 2, 2}};









// only change table_fC_x_nC and table_bC_x_nC compared to Hexaheron8

const int Hexahedron20::table_vC_x_fC[8][6]  = {{0, 1, 2, 0, 0, 0},
                                                {0, 1, 3, 3, 1, 0},
                                                {0, 3, 4, 2, 1, 0},
                                                {0, 2, 4, 1, 3, 1},
                                                {1, 2, 5, 3, 1, 0},
                                                {1, 3, 5, 2, 3, 1},
                                                {3, 4, 5, 2, 3, 2},
                                                {2, 4, 5, 2, 2, 3}};
                                          
const int Hexahedron20::table_fC_x_bC[6][4] = {{1, 5 , 3 , 0},
                                               {0, 4 , 8 , 2},
                                               {2, 9 , 7 , 1},
                                               {3, 6 , 10, 4},
                                               {5, 7 , 11, 6},
                                               {8, 10, 11, 9}};
                                          
                                                // vC         fv
const int Hexahedron20::table_fC_x_vC[6][8]  = {{0,3,2,1,   0,0,0,0},
                                                {0,1,5,4,   1,1,0,0},
                                                {0,4,7,3,   2,1,0,1},
                                                {1,2,6,5,   2,1,0,1},
                                                {2,3,7,6,   2,2,1,1},
                                                {4,5,6,7,   2,2,2,2}};

const int Hexahedron20::table_fC_x_nC[6][8] = {{0, 3, 2, 1, 9,  13, 11, 8 },
                                               {0, 1, 5, 4, 8,  12, 16, 10},
                                               {0, 4, 7, 3, 10, 17, 15, 9 },
                                               {1, 2, 6, 5, 11, 14, 18, 12},
                                               {2, 3, 7, 6, 13, 15, 19, 14},
                                               {4, 5, 6, 7, 16, 18, 19, 17}};

const int Hexahedron20::table_bC_x_vC[12][2]= {{0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},{5,6},{6,7}};

const int Hexahedron20::table_bC_x_nC[12][3] = {{0, 1, 8 }, 
                                                {0, 3, 9 }, 
                                                {0, 4, 10},
                                                {1, 2, 11},
                                                {1, 5, 12},
                                                {2, 3, 13},
                                                {2, 6, 14},
                                                {3, 7, 15},
                                                {4, 5, 16},
                                                {4, 7, 17},
                                                {5, 6, 18},
                                                {6, 7, 19}};

const int Hexahedron20::table_bC_x_fC[12][4] = {{0, 1, 3, 0},
                                                {2, 0, 3, 0},
                                                {1, 2, 3, 0},
                                                {0, 3, 2, 0},
                                                {3, 1, 3, 1},
                                                {0, 4, 1, 0},
                                                {4, 3, 3, 1},
                                                {2, 4, 2, 1},
                                                {1, 5, 2, 0},
                                                {5, 2, 3, 1},
                                                {3, 5, 2, 1},
                                                {4, 5, 2, 2}};













// only change table_fC_x_nC and table_bC_x_nC compared to Hexaheron8

const int Hexahedron27::table_vC_x_fC[8][6]  = {{0, 1, 2, 0, 0, 0},
                                                {0, 1, 3, 3, 1, 0},
                                                {0, 3, 4, 2, 1, 0},
                                                {0, 2, 4, 1, 3, 1},
                                                {1, 2, 5, 3, 1, 0},
                                                {1, 3, 5, 2, 3, 1},
                                                {3, 4, 5, 2, 3, 2},
                                                {2, 4, 5, 2, 2, 3}};
                                          
const int Hexahedron27::table_fC_x_bC[6][4] = {{1, 5 , 3 , 0},
                                               {0, 4 , 8 , 2},
                                               {2, 9 , 7 , 1},
                                               {3, 6 , 10, 4},
                                               {5, 7 , 11, 6},
                                               {8, 10, 11, 9}};
                                          
                                                // vC         fv
const int Hexahedron27::table_fC_x_vC[6][8]  = {{0,3,2,1,   0,0,0,0},
                                               {0,1,5,4,   1,1,0,0},
                                               {0,4,7,3,   2,1,0,1},
                                               {1,2,6,5,   2,1,0,1},
                                               {2,3,7,6,   2,2,1,1},
                                               {4,5,6,7,   2,2,2,2}};
// same as above
const int Hexahedron27::table_fC_x_nC[6][9] = {{0, 3, 2, 1, 9,  13, 11, 8,  20},
                                               {0, 1, 5, 4, 8,  12, 16, 10, 21},
                                               {0, 4, 7, 3, 10, 17, 15, 9,  22},
                                               {1, 2, 6, 5, 11, 14, 18, 12, 23},
                                               {2, 3, 7, 6, 13, 15, 19, 14, 24},
                                               {4, 5, 6, 7, 16, 18, 19, 17, 25}};

const int Hexahedron27::table_bC_x_vC[12][2]= {{0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},{5,6},{6,7}};

const int Hexahedron27::table_bC_x_nC[12][3] = {{0, 1, 8 }, 
                                                {0, 3, 9 }, 
                                                {0, 4, 10},
                                                {1, 2, 11},
                                                {1, 5, 12},
                                                {2, 3, 13},
                                                {2, 6, 14},
                                                {3, 7, 15},
                                                {4, 5, 16},
                                                {4, 7, 17},
                                                {5, 6, 18},
                                                {6, 7, 19}};

const int Hexahedron27::table_bC_x_fC[12][4] = {{0, 1, 3, 0},
                                                {2, 0, 3, 0},
                                                {1, 2, 3, 0},
                                                {0, 3, 2, 0},
                                                {3, 1, 3, 1},
                                                {0, 4, 1, 0},
                                                {4, 3, 3, 1},
                                                {2, 4, 2, 1},
                                                {1, 5, 2, 0},
                                                {5, 2, 3, 1},
                                                {3, 5, 2, 1},
                                                {4, 5, 2, 2}};
