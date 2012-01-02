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

#ifndef FEPIC_SHAPE_FUNCTIONS_HPP
#define FEPIC_SHAPE_FUNCTIONS_HPP

#include "shape_declarations.hpp"
#include "enums_sf.hpp"
#include "../mesh/enums.hpp"
#include "../util/common.hpp"
#include "../util/macros.hpp"
#include "../util/assert.hpp"
#include "../dofhandler/dof_traits.hpp"


class ShapeFunction
{
public:

  enum Properties {
    Part_of_unity = 0x1,
    Discontinuous = 0x2,
    Enriched      = 0x4
  };

  explicit ShapeFunction() : _tags(0) {}

protected:
  explicit ShapeFunction(unsigned tags) : _tags(tags) {}

public:

  static
  EShapeType defaultEShapeType(ECellType ct);

  static
  ShapeFunction* create(ECellType ct, EShapeType sf=UNDEFINED_SHAPE);


  /** @brief evaluates the i-th shape function at coordinate X of unitary element.
   *  @param x the coordinate.
   *  @param ith which function to evaluate.
   *  @return the evaluation result.
   */
  virtual Real operator() (Real const*x, int ith) const = 0;

  /** @brief evaluates the gradient of i-th shape function at coordinate X of unitary element.
   *  @param x the coordinate.
   *  @param ith which function to evaluate.
   *  @param c c-th component of the vector.
   *  @return vector to put the result.
   */
  virtual Real gradL (Real const*x, int ith, int c) const = 0;

  
  inline
  Real eval(Real const*x, int ith) const
  {
    return (*this).operator() (x, ith);
  }


  /** @brief evaluates the gradient of i-th shape function at coordinate X of unitary element.
   *  @param x the coordinate.
   *  @param ith which function to evaluate.
   *  @param[out] vector to put the result.
   */
  void gradL (Real const*x, int ith, Real *grad) const
  {
    for (int i = 0; i < this->getDim(); ++i)
      grad[i] = this->gradL(x, ith, i);
  }

  bool partUnity() const
  {
    return _tags & Part_of_unity;
  }
  
  bool enriched() const
  {
    return _tags & Enriched;
  }
  
  bool discontinuous() const
  {
    return _tags & Discontinuous;
  }

  /** @brief returns the number of degrees of freedom.
   */
  virtual int getNumDof() const = 0;
  virtual int getDim() const = 0;
    
  virtual int numDofsAssociatedToVertice() const = 0;//{return -1;};
  virtual int numDofsAssociatedToCorner()  const = 0;//{return -1;};
  virtual int numDofsAssociatedToFacet()   const = 0;//{return -1;};
  virtual int numDofsAssociatedToCell()    const = 0;//{return -1;};
    
    
  virtual ~ShapeFunction() {};

protected:
  unsigned _tags;

};


// =====================================================
//  ___                  _     ___  _                  
// |  _> ___ ._ _  ___ _| |_  / __>| |_  ___  ___  ___ 
// | <__/ . \| ' |<_-<  | |   \__ \| . |<_> || . \/ ._>
// `___/\___/|_|_|/__/  |_|   <___/|_|_|<___||  _/\___.
//                                           |_| 
// =====================================================
class Shape_CONST : public ShapeFunction, public _DofTraits<Shape_CONST>
{
public:

  Shape_CONST() : ShapeFunction(Discontinuous | Part_of_unity) {}

  typedef Shape_CONST ShapeType;

  //typedef Real (ShapeType::*ShapeMemFunPtr)(Real x) const;

  Real operator() (Real const*, int) const
  {
    return 1.;
  }
  Real gradL (Real const*, int, int) const
  {
    return 0.;
  }

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_CONST() {}
protected:

};


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

class Shape_EDGE_P1 : public ShapeFunction, public _DofTraits<Shape_EDGE_P1>
{
public:
  typedef Shape_EDGE_P1 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x) const;

  Shape_EDGE_P1() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}
  
  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}
  
  ~Shape_EDGE_P1() {}
protected:


};

class Shape_EDGE_P2 : public ShapeFunction, public _DofTraits<Shape_EDGE_P2>
{
public:
  typedef Shape_EDGE_P2 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x) const;

  Shape_EDGE_P2() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_EDGE_P2() {}
protected:



};

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

class Shape_TRI_P1 : public ShapeFunction, public _DofTraits<Shape_TRI_P1>
{
public:
  typedef Shape_TRI_P1 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_P1() : ShapeFunction(Part_of_unity) {}
  
  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_P1() {}
protected:



};


class Shape_TRI_P2 : public ShapeFunction, public _DofTraits<Shape_TRI_P2>
{
public:
  typedef Shape_TRI_P2 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_P2() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_P2() {}
protected:



};


class Shape_TRI_BUBBLE : public ShapeFunction, public _DofTraits<Shape_TRI_BUBBLE>
{
public:
  typedef Shape_TRI_BUBBLE ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_BUBBLE() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;


  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_BUBBLE() {}
protected:

};



class Shape_TRI_P1ph : public ShapeFunction, public _DofTraits<Shape_TRI_P1ph>
{
public:
  typedef Shape_TRI_P1ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_P1ph() : ShapeFunction(Enriched) {}
  
  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_P1ph() {}
protected:

};


class Shape_TRI_P2ph : public ShapeFunction, public _DofTraits<Shape_TRI_P2ph>
{
public:
  typedef Shape_TRI_P2ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_P2ph() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_P2ph() {}
protected:

};


class Shape_TRI_Pm1 : public ShapeFunction, public _DofTraits<Shape_TRI_Pm1>
{
public:
  typedef Shape_TRI_Pm1 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_TRI_Pm1() : ShapeFunction(Discontinuous) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TRI_Pm1() {}
protected:



};





class Shape_QUAD_Q1 : public ShapeFunction, public _DofTraits<Shape_QUAD_Q1>
{
public:
  typedef Shape_QUAD_Q1 ShapeType;
  typedef Real(ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_Q1() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_Q1() {}
protected:

  static const int _nds_pos[n_dof][dim];

};


class Shape_QUAD_Q2 : public ShapeFunction, public _DofTraits<Shape_QUAD_Q2>
{
public:
  typedef Shape_QUAD_Q2 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_Q2() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_Q2() {}
protected:

  static const int _nds_pos[n_dof][dim];

};


class Shape_QUAD_BUBBLE : public ShapeFunction, public _DofTraits<Shape_QUAD_BUBBLE>
{
public:
  typedef Shape_QUAD_BUBBLE ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_BUBBLE() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;


  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_BUBBLE() {}
protected:

};



class Shape_QUAD_Q1ph : public ShapeFunction, public _DofTraits<Shape_QUAD_Q1ph>
{
public:
  typedef Shape_QUAD_Q1ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_Q1ph() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_Q1ph() {}
protected:

};


class Shape_QUAD_Pm1 : public ShapeFunction, public _DofTraits<Shape_QUAD_Pm1>
{
public:
  typedef Shape_QUAD_Pm1 ShapeType;
  typedef Real(ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_Pm1() : ShapeFunction(Discontinuous) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_Pm1() {}
protected:

};



class Shape_QUAD_Q2ser : public ShapeFunction, public _DofTraits<Shape_QUAD_Q2ser>
{
public:
  typedef Shape_QUAD_Q2ser ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y) const;

  Shape_QUAD_Q2ser() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_QUAD_Q2ser() {}
protected:

  static const int _nds_pos[n_dof][dim];

};






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

class Shape_TET_P1 : public ShapeFunction, public _DofTraits<Shape_TET_P1>
{
public:
  typedef Shape_TET_P1 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_P1() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_P1() {}
protected:

  template<int Ith>
  Real Phi_(Real x, Real y, Real z) const;

  template<int Ith, int Coord>
  Real dPhi_(Real x, Real y, Real z) const;

};


class Shape_TET_P2 : public ShapeFunction, public _DofTraits<Shape_TET_P2>
{
public:
  typedef Shape_TET_P2 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_P2() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_P2() {}
protected:

  template<int Ith>
  Real Phi_(Real x, Real y, Real z) const;

  template<int Ith, int Coord>
  Real dPhi_(Real x, Real y, Real z) const;

};


class Shape_TET_BUBBLE : public ShapeFunction, public _DofTraits<Shape_TET_BUBBLE>
{
public:
  typedef Shape_TET_BUBBLE ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_BUBBLE() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;


  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_BUBBLE() {}
protected:

};



class Shape_TET_P1ph : public ShapeFunction, public _DofTraits<Shape_TET_P1ph>
{
public:
  typedef Shape_TET_P1ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_P1ph() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}
  
  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_P1ph() {}
protected:

};


class Shape_TET_P2ph : public ShapeFunction, public _DofTraits<Shape_TET_P2ph>
{
public:
  typedef Shape_TET_P2ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_P2ph() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_P2ph() {}
protected:

};


class Shape_TET_Pm1 : public ShapeFunction, public _DofTraits<Shape_TET_Pm1>
{
public:
  typedef Shape_TET_Pm1 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_TET_Pm1() : ShapeFunction(Discontinuous) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_TET_Pm1() {}
protected:

  template<int Ith>
  Real Phi_(Real x, Real y, Real z) const;

  template<int Ith, int Coord>
  Real dPhi_(Real x, Real y, Real z) const;

};



class Shape_HEX_Q1 : public ShapeFunction, public _DofTraits<Shape_HEX_Q1>
{
public:
  typedef Shape_HEX_Q1 ShapeType;
  typedef Real(ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_Q1() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_Q1() {}
protected:

  static const int _nds_pos[n_dof][dim];

};


// TODO
class Shape_HEX_Q2 : public ShapeFunction, public _DofTraits<Shape_HEX_Q2>
{
public:
  typedef Shape_HEX_Q2 ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_Q2() : ShapeFunction(Part_of_unity), lag() {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_Q2() {}
protected:

  static const int _nds_idx_lag[n_dof][dim];
  
  Shape_EDGE_P2 const lag;
};


class Shape_HEX_BUBBLE : public ShapeFunction, public _DofTraits<Shape_HEX_BUBBLE>
{
public:
  typedef Shape_HEX_BUBBLE ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_BUBBLE() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;


  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_BUBBLE() {}
protected:

};



class Shape_HEX_Q1ph : public ShapeFunction, public _DofTraits<Shape_HEX_Q1ph>
{
public:
  typedef Shape_HEX_Q1ph ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_Q1ph() : ShapeFunction(Enriched) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_Q1ph() {}
protected:

};


class Shape_HEX_Pm1 : public ShapeFunction, public _DofTraits<Shape_HEX_Pm1>
{
public:
  typedef Shape_HEX_Pm1 ShapeType;
  typedef Real(ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_Pm1() : ShapeFunction(Discontinuous) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_Pm1() {}
protected:

};


// TODO
class Shape_HEX_Q2ser : public ShapeFunction, public _DofTraits<Shape_HEX_Q2ser>
{
public:
  typedef Shape_HEX_Q2ser ShapeType;
  typedef Real (ShapeType::*ShapeMemFunPtr)(Real x, Real y, Real z) const;

  Shape_HEX_Q2ser() : ShapeFunction(Part_of_unity) {}

  Real operator() (Real const*x, int ith) const;
  Real gradL (Real const*x, int ith, int c) const;

  int getNumDof() const {return n_dof;}
  int getDim() const {return dim;}

  int numDofsAssociatedToVertice() const {return n_dof_per_vertice;}
  int numDofsAssociatedToCorner()  const {return n_dof_per_corner ;}
  int numDofsAssociatedToFacet()   const {return n_dof_per_facet  ;}
  int numDofsAssociatedToCell()    const {return n_dof_per_cell   ;}

  ~Shape_HEX_Q2ser() {}
protected:

  static const int _nds_pos[n_dof][dim];

};





#endif





