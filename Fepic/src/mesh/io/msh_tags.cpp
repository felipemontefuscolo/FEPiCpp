#include "msh_tags.hpp"
#include "../../util/assert.hpp"
#include "../../util/misc2.hpp"
#include <tr1/array>

namespace _MshTagsInitializers
{
  std::tr1::array<int, MSH_MAX_INDEX+1> n_nodes_4_elm_t()
  {
    std::tr1::array<int, MSH_MAX_INDEX+1> tab;
    
    tab[MSH_LIN_2  ] = 2;
    tab[MSH_TRI_3  ] = 3;
    tab[MSH_QUA_4  ] = 4;
    tab[MSH_TET_4  ] = 4;
    tab[MSH_HEX_8  ] = 8;
    tab[MSH_PRI_6  ] = 6;
    tab[MSH_PYR_5  ] = 5;
    tab[MSH_LIN_3  ] = 3;
    tab[MSH_TRI_6  ] = 6;
    tab[MSH_QUA_9  ] = 9;
    tab[MSH_TET_10 ] = 10;
    tab[MSH_HEX_27 ] = 27;
    tab[MSH_PRI_18 ] = 18;
    tab[MSH_PYR_14 ] = 14;
    tab[MSH_PNT    ] = 1 ;
    tab[MSH_QUA_8  ] = 8 ;
    tab[MSH_HEX_20 ] = 20;
    tab[MSH_PRI_15 ] = 15;
    tab[MSH_PYR_13 ] = 13;
    tab[MSH_TRI_9  ] = 9 ;
    tab[MSH_TRI_10 ] = 10;
    tab[MSH_TRI_12 ] = 12;
    tab[MSH_TRI_15 ] = 15;
    tab[MSH_TRI_15I] = 12;
    tab[MSH_TRI_21 ] = 21;
    tab[MSH_LIN_4  ] = 4 ;
    tab[MSH_LIN_5  ] = 5 ;
    tab[MSH_LIN_6  ] = 6 ;
    tab[MSH_TET_20 ] = 20;
    tab[MSH_TET_35 ] = 35;
    tab[MSH_TET_56 ] = 56;
    tab[MSH_TET_34 ] = 34;    
    
    return tab;
  }
 
  std::tr1::array<int, MSH_MAX_INDEX+1> order_4_elm_t()
  {
    std::tr1::array<int, MSH_MAX_INDEX+1> tab;
    
    tab[MSH_LIN_2  ] = 1;
    tab[MSH_TRI_3  ] = 1;
    tab[MSH_QUA_4  ] = 1;
    tab[MSH_TET_4  ] = 1;
    tab[MSH_HEX_8  ] = 1;
    tab[MSH_PRI_6  ] = 1;
    tab[MSH_PYR_5  ] = 1;
    tab[MSH_LIN_3  ] = 2;
    tab[MSH_TRI_6  ] = 2;
    tab[MSH_QUA_9  ] = 2;
    tab[MSH_TET_10 ] = 2;
    tab[MSH_HEX_27 ] = 2;
    tab[MSH_PRI_18 ] = 2;
    tab[MSH_PYR_14 ] = 2;
    tab[MSH_PNT    ] = 0;
    tab[MSH_QUA_8  ] = 2;
    tab[MSH_HEX_20 ] = 2;
    tab[MSH_PRI_15 ] = 4;
    tab[MSH_PYR_13 ] = 0;
    tab[MSH_TRI_9  ] = 3;
    tab[MSH_TRI_10 ] = 3;
    tab[MSH_TRI_12 ] = 4;
    tab[MSH_TRI_15 ] = 4;
    tab[MSH_TRI_15I] = 4;
    tab[MSH_TRI_21 ] = 5;
    tab[MSH_LIN_4  ] = 3;
    tab[MSH_LIN_5  ] = 4;
    tab[MSH_LIN_6  ] = 5;
    tab[MSH_TET_20 ] = 3;
    tab[MSH_TET_35 ] = 4;
    tab[MSH_TET_56 ] = 5;
    tab[MSH_TET_34 ] = 4;
    
    return tab;
  }
 
  std::tr1::array<int, MSH_MAX_INDEX+1> dim_4_elm_t()
  {
    std::tr1::array<int, MSH_MAX_INDEX+1> tab;
    
    tab[MSH_LIN_2  ] = 1;
    tab[MSH_TRI_3  ] = 2;
    tab[MSH_QUA_4  ] = 2;
    tab[MSH_TET_4  ] = 3;
    tab[MSH_HEX_8  ] = 3;
    tab[MSH_PRI_6  ] = 3;
    tab[MSH_PYR_5  ] = 3;
    tab[MSH_LIN_3  ] = 1;
    tab[MSH_TRI_6  ] = 2;
    tab[MSH_QUA_9  ] = 2;
    tab[MSH_TET_10 ] = 3;
    tab[MSH_HEX_27 ] = 3;
    tab[MSH_PRI_18 ] = 3;
    tab[MSH_PYR_14 ] = 3;
    tab[MSH_PNT    ] = 0;
    tab[MSH_QUA_8  ] = 2;
    tab[MSH_HEX_20 ] = 3;
    tab[MSH_PRI_15 ] = 3;
    tab[MSH_PYR_13 ] = 3;
    tab[MSH_TRI_9  ] = 2;
    tab[MSH_TRI_10 ] = 2;
    tab[MSH_TRI_12 ] = 2;
    tab[MSH_TRI_15 ] = 2;
    tab[MSH_TRI_15I] = 2;
    tab[MSH_TRI_21 ] = 2;
    tab[MSH_LIN_4  ] = 1;
    tab[MSH_LIN_5  ] = 1;
    tab[MSH_LIN_6  ] = 1;
    tab[MSH_TET_20 ] = 3;
    tab[MSH_TET_35 ] = 3;
    tab[MSH_TET_56 ] = 3;
    tab[MSH_TET_34 ] = 3;
    
    return tab;
  }
  
};


int numNodeForMshTag(EMshTag type)
{
  static const
  std::tr1::array<int, MSH_MAX_INDEX+1> nds = _MshTagsInitializers::n_nodes_4_elm_t();

  unsigned idx = static_cast<unsigned>(type);
  if (idx-1 >= MSH_MAX_INDEX)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return nds[idx];
}


int orderForMshTag(EMshTag type)
{
  static const
  std::tr1::array<int, MSH_MAX_INDEX+1> orders = _MshTagsInitializers::order_4_elm_t();
  
  unsigned idx = static_cast<unsigned>(type);
  if (idx-1 >= MSH_MAX_INDEX)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return orders[idx];
}


int dimForMshTag(EMshTag type)
{
  static const
  std::tr1::array<int, MSH_MAX_INDEX+1> dims = _MshTagsInitializers::dim_4_elm_t();

  unsigned idx = static_cast<unsigned>(type);
  if (idx-1 >= MSH_MAX_INDEX)
  {
    FEPIC_CHECK(false, "invalid cell type", std::invalid_argument);
    return -1;
  }

  return dims[idx];
}



