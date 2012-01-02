#ifndef FEPIC_ENUMS_SF_HPP
#define FEPIC_ENUMS_SF_HPP

enum EShapeType {

        UNDEFINED_SHAPE      = 0x0  ,
/* 1*/  P0                   = 0x1  ,
/* 2*/  P1                   = 0x2  ,
/* 3*/  P2                   = 0x4  ,
/* 4*/  P1ph                 = 0x8  ,     // P1+ (hierarchical)
/* 5*/  P2ph                 = 0x10 ,     // P2+ (hierarchical)
/* 6*/  Pm1                  = 0x20 ,     // P1 discontinuous
/* 7*/  Q0                   = 0x40 ,     
/* 8*/  Q1                   = 0x80 ,     
/* 9*/  Q2                   = 0x100,     
/*10*/  Q1ph                 = 0x200,     // Q1+ (hierarchical)
/*11*/  Q2ser                = 0x400,     // serendipity
/*12*/  BUBBLE               = 0x800,
        N_SHAPE_TYPES        = 12

};

inline
EShapeType facetof(EShapeType sf)
{
  if (sf == P1ph)
    return P1;
  if (sf == P2ph)
    return P2;
  if (sf == Q1ph)
    return Q1;
  
  return sf;
}

#endif

