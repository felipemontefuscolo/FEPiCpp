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
#include <Fepic/Shape>
#include <Fepic/Quadrature>
#include <limits> // for std::numeric_limits<Real>::epsilon()
#include <cmath>
#include <tr1/tuple>
#include <iostream>

typedef std::numeric_limits< Real > dbl;

Real const FEP_EPS = std::numeric_limits<Real>::epsilon();
Real const FEP_TOL = 500000*FEP_EPS; // ~ 1.1e-10 para double

//typedef std::tr1::tuple<ECellType, EShapeType, const char*> CSN_t;


TEST(QuadratureTest, WeightTest)
{
  ECellType types[] = {EDGE2, EDGE3, TRIANGLE3, TRIANGLE6, QUADRANGLE4,
                      QUADRANGLE8, QUADRANGLE9, TETRAHEDRON4, TETRAHEDRON10,
                      HEXAHEDRON8, HEXAHEDRON20, HEXAHEDRON27};
  
  Real areas[] = {2,2,0.5,0.5,4,4,4,1./6,1./6,8,8,8};
  
  int n_types = sizeof(types)/sizeof(ECellType);
  Quadrature *quadr;
  Real weight;
  

  for (int i = 0; i < n_types; ++i)
  {
    quadr = Quadrature::create(types[i]);
    
    weight = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      weight += quadr->weight(qp);
    }
    
    ASSERT_NEAR(areas[i], weight, FEP_EPS);
    
    delete quadr;
  }

}

// factorial
int _fac(int n)
{
	return n==0 ? 1 : n*_fac(n-1);
}

TEST(QuadratureTest, EdgeIntegrationTest)
{
  // f(x) = x^p;
  
  Quadrature *quadr = Quadrature::create(EDGE2);
  Real integ, exact, x;
  Real diff(0);
  int  at_p(0);
  
  for (int p = 0; p < quadr->maxOrder(); ++p)
  {
    quadr->setOrder(p);
    
    integ = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      x = *quadr->point(qp);
      integ += quadr->weight(qp) * pow(x,p);
    }
    
    exact = (pow(-1,p) + 1.0)/(p+1);
    
    ASSERT_NEAR(exact, integ, FEP_EPS)
      << "in order = " << FEP_EPS;
    
    
    
    if (diff < fabs(exact-integ))
    {
      diff = fabs(exact-integ);
      at_p = p;
    }
    
  }
  
  delete quadr;
  
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  std::cout << ">> Edge precison test: \n";
  std::cout.precision(dbl::digits10 + 1);
  std::cout << ">> max absolute error: " << diff << "\n";
  std::cout << ">> in function with order: " << at_p << "\n";
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  
};

TEST(QuadratureTest, TriIntegrationTest)
{
  // f(x) = x^p + y^p + (1-x-y)^p;
  
  Quadrature *quadr = Quadrature::create(TRIANGLE3);
  Real integ, exact, x, y, f_qp;
  Real diff(0);
  int  at_p(0);
  
  for (int p = 0; p < quadr->maxOrder(); ++p)
  {
    quadr->setOrder(p);
    
    integ = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      x = quadr->point(qp)[0];
      y = quadr->point(qp)[1];
      
      f_qp = pow(x,p) + pow(y,p) + pow(1.-x-y,p);
      
      integ += quadr->weight(qp) * f_qp;
    }
    
    exact = _fac(p) * 3./_fac(p+2);
    
    ASSERT_NEAR(exact, integ, FEP_EPS)
      << "in order = " << p;
    
    
    
    if (diff < fabs(exact-integ))
    {
      diff = fabs(exact-integ);
      at_p = p;
    }
    
  }
  
  delete quadr;
  
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  std::cout << ">> Triangle precison test: \n";
  std::cout.precision(dbl::digits10 + 1);
  std::cout << ">> max absolute error: " << diff << "\n";
  std::cout << ">> in function with order: " << at_p << "\n";
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  
};

TEST(QuadratureTest, QuadIntegrationTest)
{
  // f(x) = (x*y)^p;
  
  Quadrature *quadr = Quadrature::create(QUADRANGLE4);
  Real integ, exact, x, y, f_qp;
  Real diff(0);
  int  at_p(0);
  
  for (int p = 0; p < quadr->maxOrder(); ++p)
  {
    quadr->setOrder(p);
    
    integ = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      x = quadr->point(qp)[0];
      y = quadr->point(qp)[1];
      
      f_qp = pow(x*y,p);
      
      integ += quadr->weight(qp) * f_qp;
    }
    
    exact = pow( (pow(-1,p) + 1)/(p+1.) , 2);
    
    ASSERT_NEAR(exact, integ, 2.*FEP_EPS)
      << "in order = " << p;
    
    
    
    if (diff < fabs(exact-integ))
    {
      diff = fabs(exact-integ);
      at_p = p;
    }
    
  }
  
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  std::cout << ">> Quadrangle precison test: \n";
  std::cout.precision(dbl::digits10 + 1);
  std::cout << ">> max absolute error: " << diff << "\n";
  std::cout << ">> in function with order: " << at_p << "\n";
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";

  delete quadr;
};

TEST(QuadratureTest, TetIntegrationTest)
{
  // f(x) = x^p + y^p + z^p + (1-x-y-z)^p;
  
  Quadrature *quadr = Quadrature::create(TETRAHEDRON4);
  Real integ, exact, x, y, z, f_qp;
  Real diff(0);
  int  at_p(0);
  
  for (int p = 0; p < quadr->maxOrder(); ++p)
  {
    quadr->setOrder(p);
    
    integ = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      x = quadr->point(qp)[0];
      y = quadr->point(qp)[1];
      z = quadr->point(qp)[2];
      f_qp = pow(x,p) + pow(y,p) + pow(z,p) + pow(1.-x-y-z,p);
      
      integ += quadr->weight(qp) * f_qp;
    }
    
    exact = _fac(p) * 4./_fac(p+3);
    
    ASSERT_NEAR(exact, integ, FEP_EPS)
      << "in order = " << p;
    
    
    
    if (diff < fabs(exact-integ))
    {
      diff = fabs(exact-integ);
      at_p = p;
    }
    
  }
  
  delete quadr;
  
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  std::cout << ">> Tetrahedron precison test: \n";
  std::cout.precision(dbl::digits10 + 1);
  std::cout << ">> max absolute error: " << diff << "\n";
  std::cout << ">> in function with order: " << at_p << "\n";
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  
};

TEST(QuadratureTest, HexIntegrationTest)
{
  // f(x) = (x*y*z)^p;
  
  Quadrature *quadr = Quadrature::create(HEXAHEDRON8);
  Real integ, exact, x, y, z, f_qp;
  Real diff(0);
  int  at_p(0);
  
  for (int p = 0; p < quadr->maxOrder(); ++p)
  {
    quadr->setOrder(p);
    
    integ = 0;
    for (int qp = 0; qp < quadr->numPoints(); ++qp)
    {
      x = quadr->point(qp)[0];
      y = quadr->point(qp)[1];
      z = quadr->point(qp)[2];
      
      f_qp = pow(x*y*z,p);
      
      integ += quadr->weight(qp) * f_qp;
    }
    
    exact = pow( (pow(-1,p) + 1)/(p+1.) , 3);
    
    ASSERT_NEAR(exact, integ, 2.*FEP_EPS)
      << "in order = " << p;
    
    if (diff < fabs(exact-integ))
    {
      diff = fabs(exact-integ);
      at_p = p;
    }
    
  }
  
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";
  std::cout << ">> Hexahedron precison test: \n";
  std::cout.precision(dbl::digits10 + 1);
  std::cout << ">> max absolute error: " << diff << "\n";
  std::cout << ">> in function with order: " << at_p << "\n";
  std::cout << "-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n";

  delete quadr;
};




