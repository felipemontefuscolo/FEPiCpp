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


#ifndef FEPIC_CELL_BASE_HPP
#define FEPIC_CELL_BASE_HPP

#include "Fepic/src/mesh/labelable.hpp"

#include "Fepic/src/mesh/enums.hpp"
#include "Fepic/src/util/assert.hpp"
#include "Fepic/src/util/forward_declarations.hpp"
#include "Fepic/src/util/typedefs.hpp"


class Cell : public _Labelable
{
protected:
  int  _conn_comp_id;  // connected component to which the cell belongs.  
  
public:

  Cell() : _conn_comp_id(-1) {}

  typedef Cell* (*CreatorMemFunPtr)();

  template<class CT>
  static Cell* create()
  {
    return new CT;
  }

  static Cell* create(ECellType fep_tag);

  virtual Cell* clone() const = 0;
  virtual void copy(Cell const&) = 0;
  virtual int dim() const = 0;
  virtual EMshTag getMshTag() const = 0;
//  virtual int  getConnectedComponentId() const = 0;
  virtual int  getCornerId(int corner) const = 0;
  virtual void getCornerNodesId(int f, int *corner_nds) const = 0;
  virtual void getCornerVerticesId(int f, int *vtcs) const = 0;
  virtual int  getEdgeNodeId(int eC, int nCe) const = 0;
  virtual int  getFacetId(int facet) const = 0;
  virtual void getFacetNodesId(int f, int * facet_nds) const = 0;
  virtual void getFacetVerticesId(int f, int * vtcs) const = 0;
  virtual void getFacetCornersId(int f, int * corns) const = 0;
  virtual int  getIncidCell(int facet) const = 0;
  virtual char getIncidCellAnch(int facet) const = 0;
  virtual char getIncidCellPos(int facet) const = 0;
  virtual int  getNodeId(int const ith) const = 0;
  virtual void getNodesId(int * begin) const = 0;
  virtual void getVerticesId(int * begin) const = 0;
  virtual void getFacetsId(int * begin) const = 0;
  virtual void getCornersId(int * begin) const = 0;
  virtual bool inBoundary() const = 0;
  virtual bool isCorner(int const* vtcs, int &lrid) const = 0;
  virtual bool isFacet(int const* vtcs, int &lfid) const = 0;
  virtual bool isParametric() const = 0;
  virtual int  numCorners() const = 0;
  virtual int  numCornersPerFacet() const = 0;
  virtual int  numFacets() const = 0;
  virtual int  numNodes() const = 0;
  virtual int  numNodesPerCorner() const = 0;
  virtual int  numNodesPerFacet() const = 0;
  virtual int  numVertices() const = 0;
  virtual int  numVerticesPerFacet() const = 0;
  virtual int  numVerticesPerCorner() const = 0;
  virtual void reset() = 0;
  virtual void resetFacets() = 0;
  virtual void resetCorners() = 0;
  virtual void resetIncidCells() = 0;
  virtual void resetNodes() = 0;
//  virtual void setConnectedComponentId(int id) = 0;
  virtual void setCornerId(int corner, int cornerid) = 0;
  virtual void setFacetId(int facet, int facetid) = 0;
  virtual void setIncidCell(int facet, int icellid) = 0;
  virtual void setIncidCellAnch(int facet, int anch) = 0;
  virtual void setIncidCellPos(int facet, int pos) = 0;
  virtual void setIncidence(int facet, int icellid, char pos, char anch=0) = 0;
  virtual void setNodeId(int const ith, int const nodeid) = 0;
  virtual void setAllMembers(int const* nodes,
                             int const* corners,
                             int const* facets,
                             int const* ics,
                             int const* ics_pos,
                             int const* ics_ancs,
                             int const* conn_comp_id,
                             int const* tag,
                             int const* flags) = 0;


  virtual ~Cell() = 0;

  /* NON VIRTUAL FUNCTIONS */

  void setConnectedComponentId(int id)
  {
    _conn_comp_id = id;
  }

  int getConnectedComponentId() const
  {
    return _conn_comp_id;
  }

  virtual int table_fC_x_vC(int i, int j) const = 0;
  virtual int table_fC_x_nC(int i, int j) const = 0;
  virtual int table_vC_x_fC(int i, int j) const = 0;
  virtual int table_fC_x_bC(int i, int j) const = 0;
  virtual int table_bC_x_vC(int i, int j) const = 0;
  virtual int table_bC_x_nC(int i, int j) const = 0;
  virtual int table_bC_x_fC(int i, int j) const = 0;
  
};




#endif



