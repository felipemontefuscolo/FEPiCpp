// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <gtest/gtest.h>
#include <Fepic/Mesh>
#include <Fepic/Shape>
#include <limits> // for std::numeric_limits<Real>::epsilon()
#include <cmath>
#include <tr1/functional>
#include <functional>
#include <tr1/array>
#include <tr1/tuple>
#include <iostream>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

using namespace std::tr1::placeholders;

Real const FEP_EPS = std::numeric_limits<Real>::epsilon();
Real const FEP_TOL = 500000*FEP_EPS; // ~ 1.1e-10 para double

typedef std::tr1::tuple<ECellType, EShapeType, const char*> CSN_T;

// calculates the derivative a function f(x)
// f: callable object in form  "Real f(Real const*)"
template<class F>
Real diff_(F const& f, Real const*x, int c, int Dim)
{
	Real h = MAX(abs(x[c]), 1.)*pow(FEP_EPS, 1./3.);
	volatile Real t = h + x[c];
	h = t - x[c];

	Real y[Dim];
	Real z[Dim];
	//Real zz[Dim] = {x[0], x[1]}; // high order
	//Real yy[Dim] = {x[0], x[1]};	

  for (int i = 0; i < Dim; ++i)
  {
    y[i] = x[i];
    z[i] = x[i];
    //zz[Dim] = x[i];
    //yy[Dim] = x[i];
  }
  
	y[c]  += h;
	z[c]  -= h;
	//yy[c] += 2.*h;
	//zz[c] -= 2.*h;

  return (f(y)-f(z))/(2.*h);
	//return (-f(yy) + 8.*(f(y) - f(z)) + f(zz))/(12.*h); // for higher order
}


class ShapeTest : public testing::TestWithParam< CSN_T >
{ // vertex id
public:

  virtual void SetUp() {

    MeshIoMsh msh_reader; 
    
    ECellType   cell_t = std::tr1::get<0>(GetParam());
    EShapeType shape_t = std::tr1::get<1>(GetParam());
    const char* m_name = std::tr1::get<2>(GetParam());
    
    mesh = Mesh::create(cell_t);
    //msh_reader.readFileMsh("meshes/simptri3.msh", mesh);
    
    msh_reader.readFileMsh(m_name, mesh);

    phi = ShapeFunction::create(cell_t, shape_t);

  };

  virtual void TearDown()
  {
    delete mesh;
    delete phi;
  }

  //MeshIoVtk vtk_printer;
  Mesh *mesh;
  ShapeFunction *phi;
};

static const
CSN_T csn_table[] = {CSN_T(EDGE2,        P0    , "meshes/uni_edge.msh"),  //  0
                     CSN_T(EDGE2,        P1    , "meshes/uni_edge.msh"),  //  1
                     CSN_T(EDGE2,        P2    , "meshes/uni_edge.msh"),  //  2
                     CSN_T(TRIANGLE3,    P0    , "meshes/uni_tri.msh"),   //  3
                     CSN_T(TRIANGLE3,    P1    , "meshes/uni_tri.msh"),   //  4
                     CSN_T(TRIANGLE3,    P2    , "meshes/uni_tri.msh"),   //  5
                     CSN_T(TRIANGLE3,    P1ph  , "meshes/uni_tri.msh"),   //  6
                     CSN_T(TRIANGLE3,    P2ph  , "meshes/uni_tri.msh"),   //  7
                     CSN_T(TRIANGLE3,    Pm1   , "meshes/uni_tri.msh"),   //  8
                     CSN_T(TRIANGLE3,    BUBBLE, "meshes/uni_tri.msh"),   //  9
                     CSN_T(QUADRANGLE4,  Q0    , "meshes/uni_quad.msh"),  // 10
                     CSN_T(QUADRANGLE4,  Q1    , "meshes/uni_quad.msh"),  // 11
                     CSN_T(QUADRANGLE4,  Q2    , "meshes/uni_quad.msh"),  // 12
                     CSN_T(QUADRANGLE4,  Q1ph  , "meshes/uni_quad.msh"),  // 13
                     CSN_T(QUADRANGLE4,  Pm1   , "meshes/uni_quad.msh"),  // 14
                     CSN_T(QUADRANGLE4,  Q2ser , "meshes/uni_quad.msh"),  // 15
                     CSN_T(QUADRANGLE4,  BUBBLE, "meshes/uni_quad.msh"),  // 16
                     CSN_T(TETRAHEDRON4, P0    , "meshes/uni_tet.msh"),   // 17
                     CSN_T(TETRAHEDRON4, P1    , "meshes/uni_tet.msh"),   // 18
                     CSN_T(TETRAHEDRON4, P2    , "meshes/uni_tet.msh"),   // 19
                     CSN_T(TETRAHEDRON4, P1ph  , "meshes/uni_tet.msh"),   // 20
                     CSN_T(TETRAHEDRON4, P2ph  , "meshes/uni_tet.msh"),   // 21
                     CSN_T(TETRAHEDRON4, Pm1   , "meshes/uni_tet.msh"),   // 22
                     CSN_T(TETRAHEDRON4, BUBBLE, "meshes/uni_tet.msh"),   // 23
                     CSN_T(HEXAHEDRON8,  Q0    , "meshes/uni_hex.msh"),   // 24
                     CSN_T(HEXAHEDRON8,  Q1    , "meshes/uni_hex.msh"),   // 25
                     CSN_T(HEXAHEDRON8,  Q2    , "meshes/uni_hex.msh"),   // 26
                     CSN_T(HEXAHEDRON8,  Q1ph  , "meshes/uni_hex.msh"),   // 27
                     CSN_T(HEXAHEDRON8,  Pm1   , "meshes/uni_hex.msh"),   // 28
                     CSN_T(HEXAHEDRON8,  Q2ser , "meshes/uni_hex.msh"),   // 29
                     CSN_T(HEXAHEDRON8,  BUBBLE, "meshes/uni_hex.msh")};  // 30


TEST_P(ShapeTest, GradPrecisionTest)
{

  Real const* x;
  int  ith = 1;
  int  component = 0;
  
  int sdim = mesh->spaceDim();
  
  int const n_nodes = mesh->numNodes();
  
  // Shape function placeholder
  std::tr1::function<Real (Real const*)> Func;
  
  for (int i = 0; i < n_nodes; ++i)
  {
    x = mesh->getNodePtr(i)->getCoord();
    
    for (ith = 0; ith < phi->getNumDof(); ++ith)
    {
      Func = std::tr1::bind(std::tr1::mem_fn(&ShapeFunction::eval), phi, _1, ith);
      
      for (component = 0; component < sdim; ++component)
      {
        //ASSERT_DOUBLE_EQ(diff_( Func, x, component, sdim ),  phi->gradL(x, ith, component));
        ASSERT_NEAR(diff_( Func, x, component, sdim ),  phi->gradL(x, ith, component), FEP_TOL)
          << "component = " << component
          << "\nith(shape's local id) = " << ith;
      }
      
    }
    
  }
  
}

//// já está instanciado lá embaixo
//INSTANTIATE_TEST_CASE_P(GradientCheck, ShapeTest, ::testing::ValuesIn(csn_table));

 

TEST_P(ShapeTest, PartitionOfUnityTest)
{

  Real const* x;
  Real acc;
  int  ith = 1;
  int const n_nodes = mesh->numNodes();
  
  if (phi->partUnity())
    for (int i = 0; i < n_nodes; ++i)
    {
      x = mesh->getNodePtr(i)->getCoord();
      
      acc = 0;
      for (ith = 0; ith < phi->getNumDof(); ++ith)
        acc += phi->eval(x, ith);

      ASSERT_NEAR(Real(1),  acc, 3.*FEP_EPS)
        << "at node: " << i;        
      
    }
  
}


TEST(OutroTeste, HaHaHaHa)
{

  MeshIoMsh msh_reader;
  //MeshIoVtk vtk_printer;
  Mesh *mesh;

    
  ECellType   cell_t = TETRAHEDRON4;
  const char* m_name = "meshes/uni_tet.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(m_name, mesh);

  int iVs[128];
  int lalala[128];
  int q = -1;
  int z = -1;
  int n_nodes = mesh->numNodes();
  Point * p;
  for (int k = 0; k < n_nodes; ++k)
  {
    p = mesh->getNodePtr(k);
    int aux = (int)(mesh->connectedVtcs(p, iVs)-iVs);
    int aux2= (int)(mesh->vertexStar(p, iVs, lalala) - iVs);
    if (q < aux)
      q = aux;
    if (z < aux2)
      z = aux2;
  }
  
  //std::cout << "CONNECTED NODES: " << q << std::endl;
  //std::cout << "VERTEX STAR:     " << z << std::endl;

  delete mesh;
  
  EXPECT_TRUE(true);
}


INSTANTIATE_TEST_CASE_P(PartUnityCheck, ShapeTest, ::testing::ValuesIn(csn_table));




