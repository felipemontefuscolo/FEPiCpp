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

#ifndef FEPIC_USERCONFIG_MACROS_HPP
#define FEPIC_USERCONFIG_MACROS_HPP

/* user config */
#ifdef _OPENMP
#  define FEP_HAS_OPENMP 1
#endif

#define FEP_HAS_OMPTL 0
#define FEP_HAS_METIS 0
#define FEP_HAS_BOOST 0

#define FEPIC_DEBUG_ON

#define FEPIC_SCALAR_TYPE double

/* Print format of consecutive number. Eg, file00001.vtk */
#define FEPIC_FILE_FILL 5


/* limiters */
#define FEPIC_MAX_N_NODES_P_CELL 32
#define FEPIC_MAX_ICELLS 128

//#define FEPIC_   _INTEGER_TYPE __uint128_t

#endif
