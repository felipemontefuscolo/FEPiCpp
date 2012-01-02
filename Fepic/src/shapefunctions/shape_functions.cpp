#include "shape_functions.hpp"
#include "../util/assert.hpp"
#include "../mesh/fepic_tags.hpp"

EShapeType ShapeFunction::defaultEShapeType(ECellType ct)
{
  switch (ct)
  {
  case EDGE2:
    return P1;
    break;
  case EDGE3:
    return P2;
    break;
  case TRIANGLE3:
    return P1;
    break;
  case TRIANGLE6:
    return P2;
    break;
  case QUADRANGLE4:
    return Q1;
    break;
  case QUADRANGLE8:
    return Q2ser;
    break;
  case QUADRANGLE9:
    return Q2;
    break;
  case TETRAHEDRON4:
    return P1;
    break;
  case TETRAHEDRON10:
    return P2;
    break;
  case HEXAHEDRON8:
    return Q1;
    break;
  case HEXAHEDRON20:
    return Q2ser;
    break;
  case HEXAHEDRON27:
    return Q2;
    break;
  default:
    return UNDEFINED_SHAPE;
  }

}


ShapeFunction* ShapeFunction::create(ECellType ct, EShapeType sf)
{
  if (sf == UNDEFINED_SHAPE)
    sf = defaultEShapeType(ct);

  if (sf & (P0 | Q0))
    return new Shape_CONST;

  if (ct & (EDGE2 | EDGE3))
  {
    if (sf & (P1 | Q1 | Pm1))
      return new Shape_EDGE_P1;
    if (sf & (P2 | Q2))
      return new Shape_EDGE_P2;
  }
  else
  if (ct & (TRIANGLE3 | TRIANGLE6))
  {
    if (sf == P1)
      return new Shape_TRI_P1;
    if (sf == P2)
      return new Shape_TRI_P2;
    if (sf == P1ph)
      return new Shape_TRI_P1ph;
    if (sf == P2ph)
      return new Shape_TRI_P2ph;
    if (sf == Pm1)
      return new Shape_TRI_Pm1;
    if (sf == BUBBLE)
      return new Shape_TRI_BUBBLE;
  }
  else
  if (ct & (QUADRANGLE4 | QUADRANGLE8 | QUADRANGLE9))
  {
    if (sf == Q1)
      return new Shape_QUAD_Q1;
    if (sf == Q2)
      return new Shape_QUAD_Q2;
    if (sf == Q1ph)
      return new Shape_QUAD_Q1ph;
    if (sf == Pm1)
      return new Shape_QUAD_Pm1;
    if (sf == Q2ser)
      return new Shape_QUAD_Q2ser;
    if (sf == BUBBLE)
      return new Shape_QUAD_BUBBLE;
  }
  else
  if (ct & (TETRAHEDRON4 | TETRAHEDRON10))
  {
    if (sf == P1)
      return new Shape_TET_P1;
    if (sf == P2)
      return new Shape_TET_P2;
    if (sf == P1ph)
      return new Shape_TET_P1ph;
    if (sf == P2ph)
      return new Shape_TET_P2ph;
    if (sf == Pm1)
      return new Shape_TET_Pm1;
    if (sf == BUBBLE)
      return new Shape_TET_BUBBLE;
  }
  else
  if (ct & (HEXAHEDRON8 | HEXAHEDRON20 | HEXAHEDRON27))
  {
    if (sf == Q1)
      return new Shape_HEX_Q1;
    if (sf == Q2)
      return new Shape_HEX_Q2;
    if (sf == Q1ph)
      return new Shape_HEX_Q1ph;
    if (sf == Pm1)
      return new Shape_HEX_Pm1;
    if (sf == Q2ser)
      return new Shape_HEX_Q2ser;
    if (sf == BUBBLE)
      return new Shape_HEX_BUBBLE;
  }

  printf("cell type request: %d\n", (int)ct);
  FEPIC_CHECK(false, "invalid shape functions request", std::invalid_argument);

  return NULL;
}


//. ========================================================================
//.                                                                        .
//.                                                                        .
//.                                                                        .
//.   _ ____    ____  _                                                    .
//.  / |  _ \  / ___|| |__   __ _ _ __   ___  ___                          .
//.  | | | | | \___ \| '_ \ / _` | '_ \ / _ \/ __|                         .
//.  | | |_| |  ___) | | | | (_| | |_) |  __/\__ \                         .
//.  |_|____/  |____/|_| |_|\__,_| .__/ \___||___/                         .
//.                              |_|                                       .
//.                                                                        .
//.                                                                        .
//.                                                                        .
//. ========================================================================

// =========================== Shape_EDGE_P1 ================================

/*

        0 _____________ 1

     -1.0 ------------- 1.0
*/

Real Shape_EDGE_P1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0];

  if (ith == 0)
    return (1. - x)/2;
  else
    return (1. + x)/2;
}

Real Shape_EDGE_P1::gradL(Real const*, int ith, int) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof,
              "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith == 0)
    return -0.5;
  else
    return +0.5;

}



// =========================== Shape_EDGE_P2 ================================

/*

        0 ______2______ 1

     -1.0 ------------- 1.0
*/

Real Shape_EDGE_P2::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0];

  if (ith == 0)
    return x*(x-1.)/2;
  else
  if (ith == 1)
    return x*(x+1.)/2;
  else
    return 1.-x*x;
}

Real Shape_EDGE_P2::gradL(Real const*X, int ith, int) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof,
              "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0];

  if (ith == 0)
    return x-0.5;
  else
  if (ith == 1)
    return x+0.5;
  else
    return -2.*x;
}



//. ========================================================================
//.                                                                        .
//.                                                                        .
//.  ____  ____    ____  _                                                 .
//. |___ \|  _ \  / ___|| |__   __ _ _ __   ___  ___                       .
//.   __) | | | | \___ \| '_ \ / _` | '_ \ / _ \/ __|                      .
//.  / __/| |_| |  ___) | | | | (_| | |_) |  __/\__ \                      .
//. |_____|____/  |____/|_| |_|\__,_| .__/ \___||___/                      .
//.                                 |_|                                    .
//.                                                                        .
//.                                                                        .
//.                                                                        .
//. ========================================================================



// =========================== Shape_TRI_P1 ================================

/*

     X[1]
     ^
     |
1.0  2
|    |`\
|    |  `\
|    |    `\
|    |      `\
|    |        `\
0.0  0----------1 --> X[0]

   0.0----------1.0
*/

Real Shape_TRI_P1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith == 0)
    return 1.- X[0] - X[1];
  else
    return X[ith-1];
}



Real Shape_TRI_P1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith == 0)
    return -1;
  else
    return ith-1 == c ? 1 : 0;
}







// =========================== Shape_TRI_P2 ================================

/*

  2
  |`\
  |  `\
  5    `4
  |      `\
  |        `\
  0-----3----1


*/


Real Shape_TRI_P2::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0], y = X[1];
  Real L = 1.-x-y;

  if(ith == 0) return L*(2*L-1);
  if(ith == 1) return x*(2*x-1);
  if(ith == 2) return y*(2*y-1);
  if(ith == 3) return 4.*x*L;
  if(ith == 4) return 4.*x*y;
  else         return 4.*y*L;

}


Real Shape_TRI_P2::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x = X[0], y = X[1];
  Real L = 1.-x-y;

  if (c == 0)
  {
    if(ith == 0) return 1 - 4*L;
    if(ith == 1) return 4*x-1;
    if(ith == 2) return 0;
    if(ith == 3) return 4*(L - x);
    if(ith == 4) return 4*y;
    else         return -4*y;
  }
  else
  {
    if(ith == 0) return 1 - 4*L;
    if(ith == 1) return 0;
    if(ith == 2) return 4*y-1;
    if(ith == 3) return -4*x;
    if(ith == 4) return 4*x;
    else         return 4*(L - y);
  }


}



















// =========================== Shape_TRI_BUBBLE ================================

/*


  |`\
  |  `\
  |    `\
  |  0   `\
  |        `\
  ------------


*/

Real Shape_TRI_BUBBLE::operator() (Real const*X, int) const
{
  return 27.*X[0]*X[1]*(1. - X[0] - X[1]);
};

Real Shape_TRI_BUBBLE::gradL(Real const*X, int, int c) const
{
  if (c==0)
    return 27.*X[1]*(1. - 2*X[0] - X[1]);
  else
    return 27.*X[0]*(1. - 2*X[1] - X[0]);
}






// =========================== Shape_TRI_P1ph ================================

/*

     y
     ^
     |
     2
     |`\
     |  `\
     |    `\
     |  3   `\
     |        `\
     0----------1 --> x


*/

Real Shape_TRI_P1ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith < Shape_TRI_P1::n_dof)
    return Shape_TRI_P1().operator()(X, ith);
  else
    return Shape_TRI_BUBBLE().operator()(X, 0);
}

Real Shape_TRI_P1ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith < Shape_TRI_P1::n_dof)
    return Shape_TRI_P1().gradL(X, ith, c);
  else
    return Shape_TRI_BUBBLE().gradL(X, 0, c);
}






// =========================== Shape_TRI_P2ph ================================

/*


  2
  |`\
  |  `\
  5    `4
  |  6   `\
  |        `\
  0-----3----1



*/

Real Shape_TRI_P2ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith < Shape_TRI_P2::n_dof)
    return Shape_TRI_P2().operator()(X, ith);
  else
    return Shape_TRI_BUBBLE().operator()(X, 0);
}

Real Shape_TRI_P2ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith < Shape_TRI_P2::n_dof)
    return Shape_TRI_P2().gradL(X, ith, c);
  else
    return Shape_TRI_BUBBLE().gradL(X, 0, c);
}




// =========================== Shape_TRI_Pm1 ================================

Real Shape_TRI_Pm1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith == 0)
    return 1.;
  else
    return X[ith - 1] - 0.333333333333333L;

}


Real Shape_TRI_Pm1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith == 0)
    return 0.;
  else
    return ith-1 == c ? 1. : 0;

}









// =========================== Shape_QUAD_Q1 ================================

/*
                 y
                 ^
                 |
+1.0       3-----------2
 |         |     |     |
 |         |     |     |
 |         |     +---- | --> x
 |         |           |
 |         |           |
-1.0       0-----------1

        -1.0-----------+1.0
*/

const int Shape_QUAD_Q1::_nds_pos[4][2] = {{-1,-1}, {+1,-1}, {+1,+1}, {-1,+1}};

Real Shape_QUAD_Q1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x=X[0], y=X[1];

  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  return (1 + x*xi)*(1 + y*yi)/4;

}

Real Shape_QUAD_Q1::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x=X[0], y=X[1];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  if (c == 0)
    return xi*(1 + y*yi)/4;
  else // (c == 1)
    return yi*(1 + x*xi)/4;

}





// =========================== Shape_QUAD_Q2 ================================

/*
                 y
                 ^
                 |
+1.0       3-----6-----2
 |         |     |     |
 |         |     |     |
 |         7     8---- 5 --> x
 |         |           |
 |         |           |
-1.0       0-----4-----1

        -1.0-----------+1.0
*/

const int Shape_QUAD_Q2::_nds_pos[9][2] = {{-1,-1}, {+1,-1}, {+1,+1}, {-1,+1},
                                           { 0,-1}, {+1, 0}, { 0,+1}, {-1, 0}, {0,0} };

Real Shape_QUAD_Q2::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x=X[0], y=X[1];

  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  if (ith < 4)
    return x*y*(x + xi)*(y + yi)/4;
  else
  if (ith == 8)
    return (1 - x*x)*(1 - y*y);
  else
  if (xi == 0)
    return (1 - x*x)*(y*y + y*yi)/2;
  else
  if (yi == 0)
    return (1 - y*y)*(x*x + x*xi)/2;

  throw;
  return 0.;
}

Real Shape_QUAD_Q2::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x=X[0], y=X[1];

  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  if (c == 0)
  {
    if (ith < 4)
      return ((xi+2*x)*y*(yi+y))/4;
    else
    if (ith == 8)
      return 2*x*(y-1)*(y+1);
    else
    if (xi == 0)
      return -x*y*(yi+y);
    else
    if (yi == 0)
      return -((xi+2*x)*(y-1)*(y+1))/2;
  }
  else
  if (c == 1)
  {
    if (ith < 4)
      return (x*(xi+x)*(yi+2*y))/4;
    else
    if (ith == 8)
      return 2*(x-1)*(x+1)*y;
    else
    if (xi == 0)
      return -((x-1)*(x+1)*(yi+2*y))/2;
    else
    if (yi == 0)
      return -x*(xi+x)*y;
  }


  throw;
  return 0.;
}











// =========================== Shape_QUAD_BUBBLE ================================

Real Shape_QUAD_BUBBLE::operator() (Real const*X, int) const
{
  return Shape_QUAD_Q2().operator() (X, 8);
}
Real Shape_QUAD_BUBBLE::gradL(Real const*X, int, int c) const
{
  return Shape_QUAD_Q2().gradL(X, 8, c);
}






// =========================== Shape_QUAD_Q1ph ================================


Real Shape_QUAD_Q1ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);
  if (static_cast<unsigned>(ith)<Shape_QUAD_Q1::n_dof)
    return Shape_QUAD_Q1().operator()(X, ith);
  else
    return Shape_QUAD_BUBBLE().operator()(X, 0);
}

Real Shape_QUAD_Q1ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (static_cast<unsigned>(ith)<Shape_QUAD_Q1::n_dof)
    return Shape_QUAD_Q1().gradL(X, ith, c);
  else
    return Shape_QUAD_BUBBLE().gradL(X, 0, c);
}




// =========================== Shape_QUAD_Pm1 ================================

Real Shape_QUAD_Pm1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith == 0)
    return 1.;
  else
    return X[ith - 1];

}


Real Shape_QUAD_Pm1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith == 0)
    return 0.;
  else
    return ith-1 == c ? 1. : 0;

}







// =========================== Shape_QUAD_Q2ser ================================

/*
                 y
                 ^
                 |
+1.0       3-----6-----2
 |         |     |     |
 |         |     |     |
 |         7     +---- 5 --> x
 |         |           |
 |         |           |
-1.0       0-----4-----1

        -1.0-----------+1.0
*/

const int Shape_QUAD_Q2ser::_nds_pos[8][2] = {{-1,-1},{+1,-1},{+1,+1},{-1,+1},
                                              { 0,-1},{+1, 0},{ 0,+1},{-1, 0}};

Real Shape_QUAD_Q2ser::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x=X[0], y=X[1];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  if (ith < 4)
    return (1 + x*xi)*(1 + y*yi)*(x*xi + y*yi - 1)/4;
  else
  if (xi == 0)
    return (1 - x*x)*(1 + y*yi)/2;
  else
  if (yi == 0)
    return (1 - y*y)*(1 + x*xi)/2;

  throw;
  return 0;

}


Real Shape_QUAD_Q2ser::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x=X[0], y=X[1];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1];

  if (c == 0)
  {
    if (ith < 4)
      return (xi*(y*yi+1)*(y*yi+2*x*xi))/4;
    else
    if (xi == 0)
      return -x*(y*yi+1);
    else
    if (yi == 0)
      return -(xi*(y-1)*(y+1))/2;
  }
  else // c == 1
  {
    if (ith < 4)
      return ((x*xi+1)*yi*(2*y*yi+x*xi))/4;
    else
    if (xi == 0)
      return -((x-1)*(x+1)*yi)/2;
    else
    if (yi == 0)
      return -(x*xi+1)*y;
  }


  throw;
  return 0;
}









// =========================================================================
//.                                                                        .
//.                                                                        .
//.                                                                        .
//.  _____ ____        ____  _                                             .
//. |___ /|  _ \      / ___|| |__   __ _ _ __   ___  ___                   .
//.   |_ \| | | |     \___ \| '_ \ / _` | '_ \ / _ \/ __|                  .
//.  ___) | |_| |      ___) | | | | (_| | |_) |  __/\__ \                  .
//. |____/|____/      |____/|_| |_|\__,_| .__/ \___||___/                  .
//.                                     |_|                                .
//.                                                                        .
//.                                                                        .
//.                                                                        .
// =========================================================================







// =========================== Shape_TET_P1 ================================

/*

 *                     (0,0,1)
 *                       3
 *                      /|`\
 *                     / |  `\
 *                    /  |    `\
 *                   /   |      `\
 *                  /    |        `\
 *                 /     0__________`\  (0,1,0)
 *                /     / (0,0,0)   ,'2
 *               /    /         , '
 *              /   /       , '
 *             /  /     , '
 *            / /   , '
 *           // , '
 *          /,'
 *          1  (1,0,0)
 *
//*/


Real Shape_TET_P1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0], y = X[1], z = X[2];

  if (ith == 0)
    return 1-x-y-z;
  else
    return X[ith - 1];
}

Real Shape_TET_P1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith == 0)
    return -1;
  else
    return ith-1 == c ? 1 : 0;

}





// =========================== Shape_TET_P2 ================================

/*

 *                     (0,0,1)
 *                       3
 *                      /|`\
 *                     / |  `\
 *                    /  7    `8
 *                   /   |      `\
 *                  /    |        `\
 *                 9     0_____6____`\  (0,1,0)
 *                /     / (0,0,0)   ,'2
 *               /    /         , '
 *              /   4       5 '
 *             /  /     , '
 *            / /   , '
 *           // , '
 *          /,'
 *          1  (1,0,0)
 *
*/













Real Shape_TET_P2::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x = X[0], y = X[1], z = X[2];
  Real L = 1 - x - y - z;

  if (ith == 0) return L*(2*L-1);
  if (ith == 1) return x*(2*x-1);
  if (ith == 2) return y*(2*y-1);
  if (ith == 3) return z*(2*z-1);
  if (ith == 4) return 4.*x*L;
  if (ith == 5) return 4.*x*y;
  if (ith == 6) return 4.*y*L;
  if (ith == 7) return 4.*L*z;
  if (ith == 8) return 4.*y*z;
  if (ith == 9) return 4.*x*z;

  throw; return 0;
}

Real Shape_TET_P2::gradL(Real const*X, int ith, int c) const
{

  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x = X[0], y = X[1], z = X[2];

  if (c == 0)
  {
    if (ith == 0) return 4*z+4*y+4*x-3;
    if (ith == 1) return 4*x-1;
    if (ith == 2) return 0;
    if (ith == 3) return 0;
    if (ith == 4) return -4*(z+y+2*x-1);
    if (ith == 5) return 4*y;
    if (ith == 6) return -4*y;
    if (ith == 7) return -4*z;
    if (ith == 8) return 0;
    if (ith == 9) return 4*z;
  }
  if (c == 1)
  {
    if (ith == 0) return 4*z+4*y+4*x-3;
    if (ith == 1) return 0;
    if (ith == 2) return 4*y-1;
    if (ith == 3) return 0;
    if (ith == 4) return -4*x;
    if (ith == 5) return 4*x;
    if (ith == 6) return -4*(z+2*y+x-1);
    if (ith == 7) return -4*z;
    if (ith == 8) return 4*z;
    if (ith == 9) return 0;
  }
  else
  {
    if (ith == 0) return 4*z+4*y+4*x-3;
    if (ith == 1) return 0;
    if (ith == 2) return 0;
    if (ith == 3) return 4*z-1;
    if (ith == 4) return -4*x;
    if (ith == 5) return 0;
    if (ith == 6) return -4*y;
    if (ith == 7) return -4*(2*z+y+x-1);
    if (ith == 8) return 4*y;
    if (ith == 9) return 4*x;
  }
  throw;
  return 0;
}




// =========================== Shape_TET_BUBBLE ================================


Real Shape_TET_BUBBLE::operator() (Real const*X, int) const
{
  return 256.*X[0]*X[1]*X[2]*(1. - X[0] - X[1] - X[2]);
};

Real Shape_TET_BUBBLE::gradL(Real const*X, int, int c) const
{
  Real x = X[0], y = X[1], z = X[2], L = 1-x-y-z;
  if (c==0)
    return 256*y*z*(L-x);
  if (c==1)
    return 256*x*z*(L-y);
  else
    return 256*x*y*(L-z);
}






// =========================== Shape_TET_P1ph ================================

/*

     y
     ^
     |
     2
     |`\
     |  `\
     |    `\
     |  3   `\
     |        `\
     0----------1 --> x


*/

Real Shape_TET_P1ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);
  if (static_cast<unsigned>(ith)<Shape_TET_P1::n_dof)
    return Shape_TET_P1().operator()(X, ith);
  else
    return Shape_TET_BUBBLE().operator()(X, 0);
}

Real Shape_TET_P1ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (static_cast<unsigned>(ith)<Shape_TET_P1::n_dof)
    return Shape_TET_P1().gradL(X, ith, c);
  else
    return Shape_TET_BUBBLE().gradL(X, 0, c);
}






// =========================== Shape_TET_P2ph ================================

/*


  2
  |`\
  |  `\
  5    `4
  |  6   `\
  |        `\
  0-----3----1



*/

Real Shape_TET_P2ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);
  if (static_cast<unsigned>(ith)<Shape_TET_P2::n_dof)
    return Shape_TET_P2().operator()(X, ith);
  else
    return Shape_TET_BUBBLE().operator()(X, 0);
}

Real Shape_TET_P2ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith < Shape_TET_P2::n_dof)
    return Shape_TET_P2().gradL(X, ith, c);
  else
    return Shape_TET_BUBBLE().gradL(X, 0, c);
}




// =========================== Shape_TET_Pm1 ================================

Real Shape_TET_Pm1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith==0)
    return 1.;
  else
    return X[ith - 1] - 0.25;
}

Real Shape_TET_Pm1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith==0)
    return 0.;
  else
    return ith-1 == c ? 1. : 0.;
}








// =========================== Shape_HEX_Q1 ================================

/*
            z
     4----------7                                                         .
     |\     ^   |\                                                        .
     | \    |   | \                                                       .
     |  \   |   |  \                                                      .
     |   5------+---6                                                     .
     |   |  +-- |-- | -> y                                                .
     0---+---\--3   |                                                     .
      \  |    \  \  |                                                     .
       \ |     \  \ |                                                     .
        \|      x  \|                                                     .
         1----------2                                                     .
*/

const int Shape_HEX_Q1::_nds_pos[8][3] = {{-1,-1,-1},{+1,-1,-1},{+1,+1,-1},{-1,+1,-1},
                                          {-1,-1,+1},{+1,-1,+1},{+1,+1,+1},{-1,+1,+1}};

Real Shape_HEX_Q1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x=X[0], y=X[1], z=X[2];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1], zi = _nds_pos[ith][2];

  return (1 + x*xi)*(1 + y*yi)*(1 + z*zi)/8;

}

Real Shape_HEX_Q1::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x=X[0], y=X[1], z=X[2];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1], zi = _nds_pos[ith][2];

  if (c == 0)
  {
    return xi*(1 + y*yi)*(1 + z*zi)/8;
  }
  if (c == 1)
  {
    return yi*(1 + x*xi)*(1 + z*zi)/8;
  }
  else
    return zi*(1 + x*xi)*(1 + y*yi)/8;
}









// =========================== Shape_HEX_Q2 ================================

/*

     4----17----7                                                         .
     |\         |\                                                        .
     | 16   25  | 19                                                      .
     10 \ 22    15 \      center point: 26                                                .
     |   5----18+---6                                                     .
     |21 |      | 24|
     0---+-9----3   |                                                     .
      \  12   23 \  14                                                     .
       8 |  20    13|                                                     .
        \|         \|                                                     .
         1----11----2                                                     .
*/

const int Shape_HEX_Q2::_nds_idx_lag[27][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                               {0,0,1},{1,0,1},{1,1,1},{0,1,1},
                                               {2,0,0},{0,2,0},{0,0,2},{1,2,0},
                                               {1,0,2},{2,1,0},{1,1,2},{0,1,2},
                                               {2,0,1},{0,2,1},{1,2,1},{2,1,1},
                                               {2,2,0},{2,0,2},{0,2,2},{1,2,2},
                                               {2,1,2},{2,2,1},{2,2,2}};

Real Shape_HEX_Q2::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  return lag(X+0, _nds_idx_lag[ith][0])*
         lag(X+1, _nds_idx_lag[ith][1])*
         lag(X+2, _nds_idx_lag[ith][2]);
}

Real Shape_HEX_Q2::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);


  if (c==0)
    return lag.gradL(X+0, _nds_idx_lag[ith][0],0)*
           lag      (X+1, _nds_idx_lag[ith][1])*
           lag      (X+2, _nds_idx_lag[ith][2]);
  else
  if(c==1)
    return lag.gradL(X+1, _nds_idx_lag[ith][1],0)*
           lag      (X+2, _nds_idx_lag[ith][2])*
           lag      (X+0, _nds_idx_lag[ith][0]);
  else
    return lag.gradL(X+2, _nds_idx_lag[ith][2],0)*
           lag      (X+0, _nds_idx_lag[ith][0])*
           lag      (X+1, _nds_idx_lag[ith][1]);
}





// =========================== Shape_HEX_BUBBLE ================================

Real Shape_HEX_BUBBLE::operator() (Real const*X, int) const
{
  return (1. - X[0]*X[0])*(1. - X[1]*X[1])*(1. - X[2]*X[2]);
}
Real Shape_HEX_BUBBLE::gradL(Real const*X, int, int c) const
{
  Real x=X[0], y=X[1], z=X[2];
  if (c==0)
    return -2*x*(y*y-1)*(z*z-1);
  else
  if (c==1)
    return -2*y*(x*x-1)*(z*z-1);
  else
    return -2*z*(x*x-1)*(y*y-1);
}






// =========================== Shape_HEX_Q1ph ================================


Real Shape_HEX_Q1ph::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<int>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);
  if (static_cast<unsigned>(ith)<Shape_HEX_Q1::n_dof)
    return Shape_HEX_Q1().operator()(X, ith);
  else
    return Shape_HEX_BUBBLE().operator()(X, 0);
}

Real Shape_HEX_Q1ph::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (static_cast<unsigned>(ith)<Shape_HEX_Q1::n_dof)
    return Shape_HEX_Q1().gradL(X, ith, c);
  else
    return Shape_HEX_BUBBLE().gradL(X, 0, c);
}




// =========================== Shape_HEX_Pm1 ================================

Real Shape_HEX_Pm1::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK((unsigned)ith<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  if (ith==0)
    return 1.;
  else
    return X[ith - 1];
}

Real Shape_HEX_Pm1::gradL(Real const*, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  if (ith==0)
    return 0.;
  else
    return ith-1 == c ? 1. : 0.;
}







// =========================== Shape_HEX_Q2ser ================================

/*

     4----17----7                                                         .
     |\         |\                                                        .
     | 16       | 19                                                      .
     10 \       15 \                                                       .
     |   5----18+---6                                                     .
     |   |      |   |
     0---+-9----3   |                                                     .
      \  12      \  14                                                     .
       8 |        13|                                                     .
        \|         \|                                                     .
         1----11----2                                                     .
*/

const int Shape_HEX_Q2ser::_nds_pos[20][3] = {{-1,-1,-1},{+1,-1,-1},{+1,+1,-1},{-1,+1,-1},
                                              {-1,-1,+1},{+1,-1,+1},{+1,+1,+1},{-1,+1,+1},
                                              { 0,-1,-1},{-1, 0,-1},{-1,-1, 0},{+1, 0,-1},
                                              {+1,-1, 0},{ 0,+1,-1},{+1,+1, 0},{-1,+1, 0},
                                              { 0,-1,+1},{-1, 0,+1},{+1, 0,+1},{ 0,+1,+1}};

Real Shape_HEX_Q2ser::operator() (Real const*X, int ith) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof, "invalid index: ith="+itoa(ith), std::out_of_range);

  Real x=X[0], y=X[1], z=X[2];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1], zi = _nds_pos[ith][2];

  if (ith < 8)
    return (1 + x*xi)*(1 + y*yi)*(1 + z*zi)*(x*xi + y*yi + z*zi - 2)/8;
  else
  if (xi == 0)
    return (1 - x*x)*(1 + y*yi)*(1 + z*zi)/4;
  else
  if (yi == 0)
    return (1 - y*y)*(1 + x*xi)*(1 + z*zi)/4;
  else
  if (zi == 0)
    return (1 - z*z)*(1 + y*yi)*(1 + x*xi)/4;

  throw;
  return 0;

}

Real Shape_HEX_Q2ser::gradL(Real const*X, int ith, int c) const
{
  FEPIC_CHECK(static_cast<unsigned>(ith)<n_dof && static_cast<unsigned>(c)<dim,
              "invalid index: ith="+itoa(ith)+", c= "+itoa(c), std::out_of_range);

  Real x=X[0], y=X[1], z=X[2];
  int xi = _nds_pos[ith][0], yi = _nds_pos[ith][1], zi = _nds_pos[ith][2];

  if (ith < 8)
  {
    if (c==0)
      return (xi*(y*yi+1)*(z*zi+1)*(z*zi+y*yi+2*x*xi-1))/8;
    else
    if (c==1)
      return ((x*xi+1)*yi*(z*zi+1)*(z*zi+2*y*yi+x*xi-1))/8;
    else
      return ((x*xi+1)*(y*yi+1)*zi*(2*z*zi+y*yi+x*xi-1))/8;
  }
  else
  if (xi == 0)
  {
    if (c==0)
      return -(x*(y*yi+1)*(z*zi+1))/2;
    else
    if (c==1)
      return -((x-1)*(x+1)*yi*(z*zi+1))/4;
    else
      return -((x-1)*(x+1)*(y*yi+1)*zi)/4;
  }
  else
  if (yi == 0)
  {
    if (c==0)
      return -(xi*(y-1)*(y+1)*(z*zi+1))/4;
    else
    if (c==1)
      return -((x*xi+1)*y*(z*zi+1))/2;
    else
      return -((x*xi+1)*(y-1)*(y+1)*zi)/4;
  }
  else
  if (zi == 0)
  {
    if (c==0)
      return -(xi*(y*yi+1)*(z-1)*(z+1))/4;
    else
    if (c==1)
      return -((x*xi+1)*yi*(z-1)*(z+1))/4;
    else
      return -((x*xi+1)*(y*yi+1)*z)/2;
  }

  throw;
  return 0;
}



