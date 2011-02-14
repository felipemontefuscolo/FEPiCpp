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

#include <Fepic/Mesh>
#include <gtest/gtest.h>

typedef DefaultTraits<2, Simplex<2>> Tri2dTraits;
typedef DefaultTraits<3, Simplex<2>> Tri3dTraits;
typedef DefaultTraits<3, Simplex<3>> Tet3dTraits;

typedef Triangle<Tri2dTraits>    Tri2d;
typedef Triangle<Tri3dTraits>    Tri3d;
typedef Tetrahedron<Tet3dTraits> Tet3d;

// esperando ansiosamente a herança de construtor ser implementada no c++0x ...
template<class _Traits> class UserMesh : public iMesh<_Traits> { public:
  template<class ... Args>
  UserMesh(Args ...args) : iMesh<_Traits>(args...) {}
  UserMesh(UserMesh const&) =default;
  UserMesh() {};
};
template<class _Traits> class UserTri : public Triangle<_Traits> { public:
  template<class ...Args>
  UserTri(Args ...args) : Triangle<_Traits>(args...) {}
  UserTri(UserTri const&) = default;
  UserTri() = default;
};
template<class _Traits> class UserTet : public Tetrahedron<_Traits> { public:
  template<class ...Args>
  UserTet(Args ...args) : Tetrahedron<_Traits>(args...) {}
  UserTet(UserTet const&) = default;
  UserTet() = default;
};
template<class _Traits> class UserHalfEdge : public HalfEdge<_Traits> { public:
  template<class ... Args>
  UserHalfEdge(Args ...args) : HalfEdge<_Traits>(args...) {}
  UserHalfEdge(UserHalfEdge const&) =default;
  UserHalfEdge() {};
};
template<class _Traits> class UserHalfFace : public HalfFace<_Traits> { public:
  template<class ... Args>
  UserHalfFace(Args ...args) : HalfFace<_Traits>(args...) {}
  UserHalfFace(UserHalfFace const&) =default;
  UserHalfFace() {};
};
template<class _Traits> class UserPoint : public Point<_Traits> { public:
  template<class ... Args>
  UserPoint(Args ...args) : Point<_Traits>(args...) {}
  UserPoint(UserPoint const&) =default;
  UserPoint() {};
};

// user Traits
class UserTri2dTraits {public:
  typedef UserTri2dTraits Traits;
  typedef UserTri<Traits> CellT;
  typedef UserHalfEdge<Traits> HalfT;
  typedef UserPoint<Traits> PointT;
  typedef UserMesh <Traits> MeshT;
  enum {spacedim = 2};
};
class UserTri3dTraits {public:
  typedef UserTri3dTraits Traits;
  typedef UserTri<Traits> CellT;
  typedef UserHalfEdge<Traits> HalfT;
  typedef UserPoint<Traits> PointT;
  typedef UserMesh <Traits> MeshT;
  enum {spacedim = 3};
};
class UserTetTraits {public:
  typedef UserTetTraits Traits;
  typedef UserTet<Traits> CellT;
  typedef UserHalfEdge<Traits> HalfT;
  typedef UserPoint<Traits> PointT;
  typedef UserMesh <Traits> MeshT;
  enum {spacedim = 3};
};


/* ====================================================================== */

template<class _Traits>
class TriangleTest : public ::testing::Test
{
protected:
  
  virtual void SetUp()
  {
  }

  virtual void TearDown()
  {
  }
  
};

typedef ::testing::Types<Tri2dTraits, Tri3dTraits, UserTri2dTraits, UserTri3dTraits> TriTraitsTypes;
TYPED_TEST_CASE(TriangleTest, TriTraitsTypes);

TYPED_TEST(TriangleTest, Methods) {

  typedef typename TypeParam::CellT CellT;
  
  Eigen::VectorXi v1(3),v2(6),v3(10),vm(3);
  v1 << 1,2,3;
  v2 << 1,2,3,4,5,6;
  v3 << 1,2,3,4,5,6,7,8,9,10;
  vm << -1,-1,-1;
  
  CellT tri1;
  CellT tri2(v1, 1);
  CellT tri3(v2, 2,'a',DISABLED);
  CellT tri4(v3, 3,111, DISABLED | MARKED);

#ifdef FEPIC_DEBUG_ON
  ASSERT_ANY_THROW(CellT(v1, 2));     // incorrect order
  ASSERT_ANY_THROW(CellT(v1, 1,666)); // tag out of range
#endif  
  
  // nodes and vertices
  EXPECT_TRUE( tri1.getVertices() == vm );
  EXPECT_TRUE( tri2.getVertices() == v1 );
  EXPECT_TRUE( tri3.getNodes()    == v2 );
  EXPECT_TRUE( tri3.getVertices() == v1 );
  EXPECT_TRUE( tri4.getVertices() == v1 );
  EXPECT_TRUE( tri4.getNodes()    == v3 );
  
  // labels
  EXPECT_EQ  ( 'a', tri3.getTag() );
  EXPECT_EQ  ( 111, tri4.getTag() );
  EXPECT_TRUE( tri3.disabled() );
  EXPECT_TRUE( tri4.disabled() && tri4.marked() );
  
  // borders nodes and vertices
  EXPECT_TRUE( tri4.getBorderVertices(0) == Eigen::Vector2i(1,2) );
  EXPECT_TRUE( tri4.getBorderVertices(1) == Eigen::Vector2i(2,3) );
  EXPECT_TRUE( tri4.getBorderVertices(2) == Eigen::Vector2i(3,1) );
  EXPECT_TRUE( tri4.getBorderNodes(0) == Eigen::Vector4i(1,2,4,5) );
  EXPECT_TRUE( tri4.getBorderNodes(1) == Eigen::Vector4i(2,3,6,7) );
  EXPECT_TRUE( tri4.getBorderNodes(2) == Eigen::Vector4i(3,1,8,9) );
  
  // neighbor getOppELN
  EXPECT_TRUE( CellT::getOppELN(1) == Eigen::Vector2i(1,0) );
  EXPECT_TRUE( CellT::getOppELN(2) == Eigen::Vector3i(1,2,0) );
  EXPECT_TRUE( CellT::getOppELN(3) == Eigen::Vector4i(1,3,2,0) );
  
}


template<class _Traits>
class TetrahedronTest : public ::testing::Test
{
protected:
  
  virtual void SetUp()
  {
  }

  virtual void TearDown()
  {
  }
  
};

typedef ::testing::Types<Tet3dTraits, UserTetTraits> TetTraitsTypes;
TYPED_TEST_CASE(TetrahedronTest, TetTraitsTypes);

TYPED_TEST(TetrahedronTest, Methods) {

  typedef typename TypeParam::CellT CellT;
  
  Eigen::VectorXi v1(4),v2(10),v3(20),vm(4);
  v1 << 1,2,3,4;
  v2 << 1,2,3,4,5,6,7,8,9,10;
  v3 << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;
  vm << -1,-1,-1,-1;
  
  CellT tet1;
  CellT tet2(v1 ,1);
  CellT tet3(v2, 2,'a',DISABLED);
  CellT tet4(v3, 3,111, DISABLED | MARKED);

#ifdef FEPIC_DEBUG_ON
  ASSERT_ANY_THROW(CellT(v1, 2));   // incorrect order
  ASSERT_ANY_THROW(CellT(v1, 1,666)); // tag out of range
#endif  
  
  // nodes and vertices
  EXPECT_TRUE( tet1.getVertices() == vm );
  EXPECT_TRUE( tet2.getVertices() == v1 );
  EXPECT_TRUE( tet3.getNodes()    == v2 );
  EXPECT_TRUE( tet3.getVertices() == v1 );
  EXPECT_TRUE( tet4.getVertices() == v1 );
  EXPECT_TRUE( tet4.getNodes()    == v3 );
  
  // labels
  EXPECT_EQ  ( 'a', tet3.getTag() );
  EXPECT_EQ  ( 111, tet4.getTag() );
  EXPECT_TRUE( tet3.disabled() );
  EXPECT_TRUE( tet4.disabled() && tet4.marked() );
  
  // borders nodes and vertices
  Eigen::VectorXi v10(10);
  EXPECT_TRUE( tet4.getBorderVertices(0) == Eigen::Vector3i(2,1,3) );
  EXPECT_TRUE( tet4.getBorderVertices(1) == Eigen::Vector3i(1,2,4) );
  EXPECT_TRUE( tet4.getBorderVertices(2) == Eigen::Vector3i(4,3,1) );
  EXPECT_TRUE( tet4.getBorderVertices(3) == Eigen::Vector3i(3,4,2) );
  
  v10 << 2,1,3,6,5,10,9,8,7,17;
  EXPECT_TRUE( tet4.getBorderNodes(0) == v10 );
  v10 << 1,2,4,5,6,16,15,11,12,18;
  EXPECT_TRUE( tet4.getBorderNodes(1) == v10 );
  v10 << 4,3,1,13,14,9,10,12,11,19;
  EXPECT_TRUE( tet4.getBorderNodes(2) == v10 );
  v10 << 3,4,2,14,13,15,16,7,8,20;
  EXPECT_TRUE( tet4.getBorderNodes(3) == v10 );
  
  
  // neighbor getOppELN
  EXPECT_TRUE( CellT::getOppELN(1) == Eigen::Vector2i(1,0) );
  EXPECT_TRUE( CellT::getOppELN(2) == Eigen::Vector3i(1,2,0) );
  EXPECT_TRUE( CellT::getOppELN(3) == Eigen::Vector4i(1,3,2,0) );
  
  // neighbor getOppFLN
  Eigen::Matrix3i order1;
  order1 <<  2,1,0,  // anchor 0
             1,0,2,  // anchor 1
             0,2,1; // anchor 2 
                 
  Eigen::MatrixXi order2(3,6);
  
  order2 << 2,1,0,4,3,5,  // anchor 0
            1,0,2,3,5,4,  // anchor 1
            0,2,1,5,4,3; // anchor 2
  
  Eigen::MatrixXi order3(3,10); 
  order3 << 2,1,0,6,5,4,3,8,7,9,  // anchor 0
            1,0,2,4,3,8,7,6,5,9,  // anchor 1
            0,2,1,8,7,6,5,4,3,9; // anchor 2
  
  
  EXPECT_TRUE( CellT::getOppFLN(1) == order1 );
  EXPECT_TRUE( CellT::getOppFLN(2) == order2 );
  EXPECT_TRUE( CellT::getOppFLN(3) == order3 );
  
  
  
}












