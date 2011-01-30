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

#ifndef MSH_TAGS_HPP
#define MSH_TAGS_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>

// Element types in .msh file format
#define MSH_LIN_2  1    // 2-node line.
#define MSH_TRI_3  2    // 3-node triangle.
#define MSH_QUA_4  3    // 4-node quadrangle.
#define MSH_TET_4  4    // 4-node tetrahedron.
#define MSH_HEX_8  5    // 8-node hexahedron.
#define MSH_PRI_6  6    // 6-node prism.
#define MSH_PYR_5  7    // 5-node pyramid.
#define MSH_LIN_3  8    // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
#define MSH_TRI_6  9    // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
#define MSH_QUA_9  10   // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
#define MSH_TET_10 11   // 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
#define MSH_HEX_27 12   // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
#define MSH_PRI_18 13   // 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
#define MSH_PYR_14 14   // 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
#define MSH_PNT    15   // 1-node point.
#define MSH_QUA_8  16   // 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
#define MSH_HEX_20 17   // 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
#define MSH_PRI_15 18   // 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
#define MSH_PYR_13 19   // 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
#define MSH_TRI_9  20   // 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
#define MSH_TRI_10 21   // 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
#define MSH_TRI_12 22   // 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
#define MSH_TRI_15 23   // 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
#define MSH_TRI_15I 24  // 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
#define MSH_TRI_21 25   // 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
#define MSH_LIN_4  26   // 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
#define MSH_LIN_5  27   // 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
#define MSH_LIN_6  28   // 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
#define MSH_TET_20 29   // 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
#define MSH_TET_35 30   // 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
#define MSH_TET_56 31   // 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
#define MSH_TET_34 32
#define MSH_UNDEFINED_ELEM 33
#define MSH_MAX_ELEMENT_NODES 56



static
int getNumVerticesForElementTypeMSH(int type)
{
	switch (type) {
	case MSH_PNT    : return 1; break;
	case MSH_LIN_2  : return 2; break;
	case MSH_LIN_3  : return 3; break;
	case MSH_LIN_4  : return 4; break;
	case MSH_LIN_5  : return 5; break;
	case MSH_LIN_6  : return 6; break;
	case MSH_TRI_3  : return 3; break;
	case MSH_TRI_6  : return 6; break;
	case MSH_TRI_9  : return 9; break;
	case MSH_TRI_10  : return 10; break;
	case MSH_TRI_12  : return 12; break;
	case MSH_TRI_15  : return 15; break;
	case MSH_TRI_15I : return 15; break;
	case MSH_TRI_21 : return 21; break;
	case MSH_QUA_4  : return 4; break;
	case MSH_QUA_8  : return 8; break;
	case MSH_QUA_9  : return 9; break;
	case MSH_HEX_8  : return 8; break;
	case MSH_HEX_20 : return 20; break;
	case MSH_HEX_27 : return 27; break;
	case MSH_PRI_6  : return 6; break;
	case MSH_PRI_15 : return 15; break;
	case MSH_PRI_18 : return 18; break;
	case MSH_PYR_5  : return 5; break;
	case MSH_PYR_13 : return 13; break;
	case MSH_PYR_14 : return 14; break;
	case MSH_TET_4  : return 4; break;
	case MSH_TET_10 : return 10; break;
	case MSH_TET_20 : return 20; break;
	case MSH_TET_35 : return 35; break;
	case MSH_TET_34 : return 34; break;
	case MSH_TET_56 : return 56; break;
	default:
		printf("Error (getMSHTags.h): Unknown type of element %d\n", type);
		throw;
	}
	return 0;
}

static
int getOrderForElementTypeMSH(int type)
{
	switch (type) {
	case MSH_LIN_2  : return 1;  
	case MSH_LIN_3  : return 2;  
	case MSH_LIN_4  : return 3;  
	case MSH_LIN_5  : return 4;  
	case MSH_LIN_6  : return 5;
	case MSH_TRI_3  : return 1;  
	case MSH_TRI_6  : return 2;  
	case MSH_TRI_10  :return 3;   
	case MSH_TRI_15  :return 4;   
	case MSH_TRI_21 : return 5;  
	case MSH_QUA_4  : return 1;  
	case MSH_QUA_9  : return 2;
	case MSH_TET_4  : return 1;  
	case MSH_TET_10 : return 2;  
	case MSH_HEX_8  : return 1;  
	case MSH_HEX_27 : return 2;  
	case MSH_TET_20 : return 3;  
	case MSH_TET_35 : return 4;  
	case MSH_TET_56 : return 5;
	default:
		std::cout << "MSH Error:Invalid order request: "<<type<<"\n";
	}
	return -1;	
}

static
int getDimForElementTypeMSH(int type)
{
	switch (type) {
	case MSH_PNT    : return 0;
	case MSH_LIN_2  :
	case MSH_LIN_3  :
	case MSH_LIN_4  :
	case MSH_LIN_5  :
	case MSH_LIN_6  : return 1;
	case MSH_TRI_3  :
	case MSH_TRI_6  :
	case MSH_TRI_9  :
	case MSH_TRI_10  :
	case MSH_TRI_12  :
	case MSH_TRI_15  :
	case MSH_TRI_15I :
	case MSH_TRI_21 :
	case MSH_QUA_4  :
	case MSH_QUA_8  :
	case MSH_QUA_9  : return 2;
	case MSH_TET_4  :
	case MSH_TET_10 :
	case MSH_HEX_8  :
	case MSH_HEX_20 :
	case MSH_HEX_27 :
	case MSH_PRI_6  :
	case MSH_PRI_15 :
	case MSH_PRI_18 :
	case MSH_PYR_5  :
	case MSH_PYR_13 :
	case MSH_PYR_14 :
	case MSH_TET_20 :
	case MSH_TET_35 :
	case MSH_TET_34 :
	case MSH_TET_56 : return 3;
	default:
		printf("Error (getMSHTags.h): Unknown type of element %d\n", type);
	}
	return -1;
}

static
std::string getElementNameMSH(int type)
{
	std::string name;
	switch (type) {
	case MSH_PNT     : name = "Point"; break;
	case MSH_LIN_2   : name = "Linear edge"; break;
	case MSH_LIN_3   : name = "Quadratic edge"; break;
	case MSH_LIN_4   : name = "Cubic edge"; break;
	case MSH_LIN_5   : name = "Quartic edge"; break;
	case MSH_LIN_6   : name = "Pentic edge"; break;
	case MSH_TRI_3   : name = "Linear triangle"; break;
	case MSH_TRI_6   : name = "Quadratic triangle"; break;
	case MSH_TRI_9   : name = "Cubic serendipity triangle"; break;
	case MSH_TRI_10  : name = "Cubic triangle"; break;
	case MSH_TRI_12  : name = "Quartic serendipity triangle"; break;
	case MSH_TRI_15  : name = "Quartic triangle"; break;
	case MSH_TRI_15I : name = "Pentic serendipity triangle"; break;
	case MSH_TRI_21  : name = "Pentic triangle"; break;
	case MSH_QUA_4   : name = "Bilinear Quadrangle"; break;
	case MSH_QUA_8   : name = "Quadratic serendipity quadrangle"; break;
	case MSH_QUA_9   : name = "Quadratic quadrangle"; break;
	case MSH_TET_4   : name = "Linear tetrahedron"; break;
	case MSH_TET_10  : name = "Quadratic tetrahedron"; break;
	case MSH_HEX_8   : name = "Trilinear hexahedron"; break;
	case MSH_HEX_20  : name = "Quadratic edge serendipity hexahedron"; break;
	case MSH_HEX_27  : name = "Quadratic serendipity hexahedron"; break;
	case MSH_PRI_6   : name = "Linear prism"; break;
	case MSH_PRI_15  : name = "Quadratic edge serendipity prism"; break;
	case MSH_PRI_18  : name = "Quadratic serendipity prism"; break;
	case MSH_PYR_5   : name = "Linear pyramid"; break;
	case MSH_PYR_13  : name = "Quadratic edge serendipity pyramid"; break;
	case MSH_PYR_14  : name = "Quadratic serendipty pyramid"; break;
	case MSH_TET_20  : name = "Cubic tetrahedron";break;
	case MSH_TET_35  : name = "Quartic tetrahedron";break;
	case MSH_TET_34  : name = "Quartic serendipity tetrahedron";break;
	case MSH_TET_56  : name = "Pentic tetrahedron";break;
	case MSH_UNDEFINED_ELEM : name = "Undefined Element"; break;
	default:
	  std::stringstream ss; break;
	  ss <<  "Unknown type of element (tag " << type << ")"; break;
	  name = ss.str();
	}
return name;
}


#endif
