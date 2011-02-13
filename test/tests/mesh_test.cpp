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
  //template<class ... Args>
  //UserMesh(Args ...args) : iMesh<_Traits>(args...) {}
  //UserMesh(UserMesh const&) =delete;
  //UserMesh()=default {};
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
  typedef UserMesh<Traits> MeshT;
  enum {spacedim = 2};
};
class UserTri3dTraits {public:
  typedef UserTri3dTraits Traits;
  typedef UserTri<Traits> CellT;
  typedef UserHalfEdge<Traits> HalfT;
  typedef UserPoint<Traits> PointT;
  typedef UserMesh<Traits> MeshT;
  enum {spacedim = 3};
};
class UserTetTraits {public:
  typedef UserTetTraits Traits;
  typedef UserTet<Traits> CellT;
  typedef UserHalfEdge<Traits> HalfT;
  typedef UserPoint<Traits> PointT;
  typedef UserMesh<Traits> MeshT;
  enum {spacedim = 3};
};


template<class _Traits>
class MeshTest : public ::testing::Test
{
protected:
  
  virtual void SetUp()
  {
    
  }

  virtual void TearDown()
  {
  }
 
  typename _Traits::MeshT mesh1;
};


typedef ::testing::Types<Tri2dTraits, Tri3dTraits, Tet3dTraits, UserTri3dTraits> MeshTraitsTypes;
TYPED_TEST_CASE(MeshTest, MeshTraitsTypes);

TYPED_TEST(MeshTest, ReadFromMsh) {
  
  typedef typename TypeParam::MeshT MeshT;
  typedef typename TypeParam::CellT CellT;
  typedef typename MeshT::CellIterator CellIterator;
  
  const int spacedim = TypeParam::spacedim;
  
  if (spacedim==3)
  {
    this->mesh1.readFileMsh("meshes/simptet.msh");
  }
  else if (spacedim==2)
    this->mesh1.readFileMsh("meshes/simptri.msh");
  else throw;

  if (spacedim==3)
  this->mesh1.setOutputFileName("meshes/out/simptet");
  if (spacedim==2)
  this->mesh1.setOutputFileName("meshes/out/simptri");

  this->mesh1.writeVtk();
  this->mesh1.setFamilyFiles();
  this->mesh1.writeVtk();
  this->mesh1.writeVtk();
  this->mesh1.setFamilyFiles(false);
  
  this->mesh1.writeFsf();
  this->mesh1.setFamilyFiles();
  this->mesh1.writeFsf();
  this->mesh1.writeFsf();
  
}
















