#include <tr1/array>
#include "fepic_tags.hpp"
#include "../util/misc2.hpp"
#include "../util/assert.hpp"

namespace _FepicTagsInitializers
{
  std::tr1::array<int, N_CELL_TYPES> ctypeNumNodes_ini()
  {
    std::tr1::array<int, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = 2;
    tab[log2_i32(EDGE2        ) ] = 2;
    tab[log2_i32(EDGE3        ) ] = 3;
    tab[log2_i32(TRIANGLE3    ) ] = 3;
    tab[log2_i32(TRIANGLE6    ) ] = 6;
    tab[log2_i32(QUADRANGLE4  ) ] = 4;
    tab[log2_i32(QUADRANGLE8  ) ] = 8;
    tab[log2_i32(QUADRANGLE9  ) ] = 9;
    tab[log2_i32(TETRAHEDRON4 ) ] = 4;
    tab[log2_i32(TETRAHEDRON10) ] = 10;
    tab[log2_i32(HEXAHEDRON8  ) ] = 8;
    tab[log2_i32(HEXAHEDRON20 ) ] = 20;
    tab[log2_i32(HEXAHEDRON27 ) ] = 27;

    return tab;
  }

  std::tr1::array<int, N_CELL_TYPES> ctypeDegree_ini()
  {
    std::tr1::array<int, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = 1;
    tab[log2_i32(EDGE2        ) ] = 1;
    tab[log2_i32(EDGE3        ) ] = 2;
    tab[log2_i32(TRIANGLE3    ) ] = 1;
    tab[log2_i32(TRIANGLE6    ) ] = 2;
    tab[log2_i32(QUADRANGLE4  ) ] = 1;
    tab[log2_i32(QUADRANGLE8  ) ] = 2;
    tab[log2_i32(QUADRANGLE9  ) ] = 2;
    tab[log2_i32(TETRAHEDRON4 ) ] = 1;
    tab[log2_i32(TETRAHEDRON10) ] = 2;
    tab[log2_i32(HEXAHEDRON8  ) ] = 1;
    tab[log2_i32(HEXAHEDRON20 ) ] = 2;
    tab[log2_i32(HEXAHEDRON27 ) ] = 2;

    return tab;
  }

  std::tr1::array<int, N_CELL_TYPES> ctypeDim_ini()
  {
    std::tr1::array<int, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = 0;
    tab[log2_i32(EDGE2        ) ] = 1;
    tab[log2_i32(EDGE3        ) ] = 1;
    tab[log2_i32(TRIANGLE3    ) ] = 2;
    tab[log2_i32(TRIANGLE6    ) ] = 2;
    tab[log2_i32(QUADRANGLE4  ) ] = 2;
    tab[log2_i32(QUADRANGLE8  ) ] = 2;
    tab[log2_i32(QUADRANGLE9  ) ] = 2;
    tab[log2_i32(TETRAHEDRON4 ) ] = 3;
    tab[log2_i32(TETRAHEDRON10) ] = 3;
    tab[log2_i32(HEXAHEDRON8  ) ] = 3;
    tab[log2_i32(HEXAHEDRON20 ) ] = 3;
    tab[log2_i32(HEXAHEDRON27 ) ] = 3;

    return tab;
  }

  std::tr1::array<const char*, N_CELL_TYPES> ctypeName_ini()
  {
    std::tr1::array<const char*, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = "Point";
    tab[log2_i32(EDGE2        ) ] = "Linear edge";
    tab[log2_i32(EDGE3        ) ] = "Quadratic edge";
    tab[log2_i32(TRIANGLE3    ) ] = "Linear triangle";
    tab[log2_i32(TRIANGLE6    ) ] = "Quadratic triangle";
    tab[log2_i32(QUADRANGLE4  ) ] = "Bilinear Quadrangle";
    tab[log2_i32(QUADRANGLE8  ) ] = "Quadratic serendipity quadrangle";
    tab[log2_i32(QUADRANGLE9  ) ] = "Quadratic quadrangle";
    tab[log2_i32(TETRAHEDRON4 ) ] = "Linear tetrahedron";
    tab[log2_i32(TETRAHEDRON10) ] = "Quadratic tetrahedron";
    tab[log2_i32(HEXAHEDRON8  ) ] = "Trilinear hexahedron";
    tab[log2_i32(HEXAHEDRON20 ) ] = "Quadratic serendipity hexahedron";
    tab[log2_i32(HEXAHEDRON27 ) ] = "Quadratic hexahedron";

    return tab;
  }

  std::tr1::array<ECellType, MSH_MAX_INDEX+1> mshTag2ctype_ini()
  {
    std::tr1::array<ECellType, MSH_MAX_INDEX+1> tab;

    tab[MSH_PNT   ] = POINT1;
    tab[MSH_LIN_2 ] = EDGE2;
    tab[MSH_TRI_3 ] = TRIANGLE3;
    tab[MSH_QUA_4 ] = QUADRANGLE4;
    tab[MSH_TET_4 ] = TETRAHEDRON4;
    tab[MSH_HEX_8 ] = HEXAHEDRON8;
    tab[MSH_LIN_3 ] = EDGE3;
    tab[MSH_TRI_6 ] = TRIANGLE6;
    tab[MSH_QUA_9 ] = QUADRANGLE9;
    tab[MSH_TET_10] = TETRAHEDRON10;
    tab[MSH_HEX_27] = HEXAHEDRON27;
    tab[MSH_QUA_8 ] = QUADRANGLE8;
    tab[MSH_HEX_20] = HEXAHEDRON20;

    return tab;
  }

  std::tr1::array<EMshTag, N_CELL_TYPES> ctype2mshTag_ini()
  {
    std::tr1::array<EMshTag, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = MSH_PNT;
    tab[log2_i32(EDGE2        ) ] = MSH_LIN_2 ;
    tab[log2_i32(EDGE3        ) ] = MSH_LIN_3 ;
    tab[log2_i32(TRIANGLE3    ) ] = MSH_TRI_3 ;
    tab[log2_i32(TRIANGLE6    ) ] = MSH_TRI_6 ;
    tab[log2_i32(QUADRANGLE4  ) ] = MSH_QUA_4 ;
    tab[log2_i32(QUADRANGLE8  ) ] = MSH_QUA_8 ;
    tab[log2_i32(QUADRANGLE9  ) ] = MSH_QUA_9 ;
    tab[log2_i32(TETRAHEDRON4 ) ] = MSH_TET_4 ;
    tab[log2_i32(TETRAHEDRON10) ] = MSH_TET_10;
    tab[log2_i32(HEXAHEDRON8  ) ] = MSH_HEX_8 ;
    tab[log2_i32(HEXAHEDRON20 ) ] = MSH_HEX_20;
    tab[log2_i32(HEXAHEDRON27 ) ] = MSH_HEX_27;

    return tab;
  }

  std::tr1::array<ECellClass, N_CELL_TYPES> ctype2cclass_ini()
  {
    std::tr1::array<ECellClass, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = POINT ;
    tab[log2_i32(EDGE2        ) ] = EDGE ;
    tab[log2_i32(EDGE3        ) ] = EDGE ;
    tab[log2_i32(TRIANGLE3    ) ] = TRIANGLE;
    tab[log2_i32(TRIANGLE6    ) ] = TRIANGLE;
    tab[log2_i32(QUADRANGLE4  ) ] = QUADRANGLE;
    tab[log2_i32(QUADRANGLE8  ) ] = QUADRANGLE;
    tab[log2_i32(QUADRANGLE9  ) ] = QUADRANGLE;
    tab[log2_i32(TETRAHEDRON4 ) ] = TETRAHEDRON;
    tab[log2_i32(TETRAHEDRON10) ] = TETRAHEDRON;
    tab[log2_i32(HEXAHEDRON8  ) ] = HEXAHEDRON;
    tab[log2_i32(HEXAHEDRON20 ) ] = HEXAHEDRON;
    tab[log2_i32(HEXAHEDRON27 ) ] = HEXAHEDRON;

    return tab;
  }

  std::tr1::array<ECellFamily, N_CELL_TYPES> ctype2cfamily_ini()
  {
    std::tr1::array<ECellFamily, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = static_cast<ECellFamily>(SIMPLEX | HCUBE);
    tab[log2_i32(EDGE2        ) ] = static_cast<ECellFamily>(SIMPLEX | HCUBE);
    tab[log2_i32(EDGE3        ) ] = static_cast<ECellFamily>(SIMPLEX | HCUBE);
    tab[log2_i32(TRIANGLE3    ) ] = SIMPLEX;
    tab[log2_i32(TRIANGLE6    ) ] = SIMPLEX;
    tab[log2_i32(QUADRANGLE4  ) ] = HCUBE;
    tab[log2_i32(QUADRANGLE8  ) ] = HCUBE;
    tab[log2_i32(QUADRANGLE9  ) ] = HCUBE;
    tab[log2_i32(TETRAHEDRON4 ) ] = SIMPLEX;
    tab[log2_i32(TETRAHEDRON10) ] = SIMPLEX;
    tab[log2_i32(HEXAHEDRON8  ) ] = HCUBE;
    tab[log2_i32(HEXAHEDRON20 ) ] = HCUBE;
    tab[log2_i32(HEXAHEDRON27 ) ] = HCUBE;

    return tab;
  }

  std::tr1::array<ECellFamily, N_CELL_CLASSES> cclass2cfamily_ini()
  {
    std::tr1::array<ECellFamily, N_CELL_CLASSES> tab;

    tab[log2_i32(POINT      ) ] = static_cast<ECellFamily>(SIMPLEX | HCUBE);
    tab[log2_i32(EDGE       ) ] = static_cast<ECellFamily>(SIMPLEX | HCUBE);
    tab[log2_i32(TRIANGLE   ) ] = SIMPLEX;
    tab[log2_i32(QUADRANGLE ) ] = HCUBE;
    tab[log2_i32(TETRAHEDRON) ] = SIMPLEX;
    tab[log2_i32(HEXAHEDRON ) ] = HCUBE;

    return tab;
  }

  std::tr1::array<ECellType, N_CELL_TYPES> facetof_ini()
  {
    std::tr1::array<ECellType, N_CELL_TYPES> tab;

    tab[log2_i32(POINT1       ) ] = UNDEFINED_CELLT;
    tab[log2_i32(EDGE2        ) ] = POINT1;
    tab[log2_i32(EDGE3        ) ] = POINT1;
    tab[log2_i32(TRIANGLE3    ) ] = EDGE2;
    tab[log2_i32(TRIANGLE6    ) ] = EDGE3;
    tab[log2_i32(QUADRANGLE4  ) ] = EDGE2;
    tab[log2_i32(QUADRANGLE8  ) ] = EDGE3;
    tab[log2_i32(QUADRANGLE9  ) ] = EDGE3;
    tab[log2_i32(TETRAHEDRON4 ) ] = TRIANGLE3;
    tab[log2_i32(TETRAHEDRON10) ] = TRIANGLE6;
    tab[log2_i32(HEXAHEDRON8  ) ] = QUADRANGLE4;
    tab[log2_i32(HEXAHEDRON20 ) ] = QUADRANGLE8;
    tab[log2_i32(HEXAHEDRON27 ) ] = QUADRANGLE9;

    return tab;
  }

}

int ctypeNumNodes(ECellType type)
{
  static const
  std::tr1::array<int, N_CELL_TYPES> n_nds = _FepicTagsInitializers::ctypeNumNodes_ini();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return n_nds[idx];

}

int ctypeDegree(ECellType type)
{

  static const
  std::tr1::array<int, N_CELL_TYPES> orders = _FepicTagsInitializers::ctypeDegree_ini();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return orders[idx];
}

int ctypeDim(ECellType type)
{

  static const
  std::tr1::array<int, N_CELL_TYPES> dims = _FepicTagsInitializers::ctypeDim_ini();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return dims[idx];
}

const char* ctypeName(ECellType type)
{
  static const
  std::tr1::array<const char*, N_CELL_TYPES> names = _FepicTagsInitializers::ctypeName_ini();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return NULL;
  }

  return names[idx];
}

ECellType mshTag2ctype(EMshTag type)
{
  static const
  std::tr1::array<ECellType, MSH_MAX_INDEX+1> tab = _FepicTagsInitializers::mshTag2ctype_ini();

  unsigned idx = static_cast<unsigned>(type);
  if (idx-1 >= MSH_MAX_INDEX)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return UNDEFINED_CELLT;
  }

  return tab[idx];

}

EMshTag ctype2mshTag(ECellType type)
{
  static const
  std::tr1::array<EMshTag, N_CELL_TYPES> tab = _FepicTagsInitializers::ctype2mshTag_ini();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return MSH_UNDEFINED_ELEM;
  }

  return tab[idx];
}


ECellClass ctype2cclass(ECellType type)
{
  std::tr1::array<ECellClass, N_CELL_TYPES> tab = _FepicTagsInitializers::ctype2cclass_ini();
  
  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return UNDEFINED_CELLC;
  }
  return tab[idx];
}

ECellFamily ctype2cfamily(ECellType type)
{
  std::tr1::array<ECellFamily, N_CELL_TYPES> tab = _FepicTagsInitializers::ctype2cfamily_ini();
  
  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return UNDEFINED_CELLF;
  }
  return tab[idx];
}

ECellFamily cclass2cfamily(ECellClass type)
{
  std::tr1::array<ECellFamily, N_CELL_CLASSES> tab = _FepicTagsInitializers::cclass2cfamily_ini();
  
  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_CLASSES)
  {
    FEPIC_CHECK(false, "invalid cell class", std::invalid_argument);
    return UNDEFINED_CELLF;
  }
  return tab[idx];
}

ECellType facetof(ECellType type)
{
  std::tr1::array<ECellType, N_CELL_TYPES> tab = _FepicTagsInitializers::facetof_ini();
  
  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return UNDEFINED_CELLT;
  }
  return tab[idx];
}

