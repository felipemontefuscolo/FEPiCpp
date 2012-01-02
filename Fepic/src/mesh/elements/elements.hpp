#ifndef FEPIC_ELEMENTS_HPP
#define FEPIC_ELEMENTS_HPP

template<typename CellT>
class _CellCore;

#include "../../util/common.hpp" // needed by point.hpp
#include "../../util/forward_declarations.hpp"
#include "../../util/typedefs.hpp"
#include "../io/msh_tags.hpp"
#include "../fepic_tags.hpp"

#include "../labelable.hpp"

#include "facet.hpp"
#include "corner.hpp"
#include "edge.hpp"
#include "hexahedron.hpp"
#include "point.hpp"
#include "quadrangle.hpp"
#include "tetrahedron.hpp"
#include "triangle.hpp"

#endif


