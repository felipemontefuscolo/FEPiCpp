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

typedef DefaultTraits<2, Simplex<2>> Tri2dTraits;
typedef DefaultTraits<3, Simplex<2>> Tri3dTraits;
typedef DefaultTraits<3, Simplex<3>> Tet3dTraits;

template<class _Traits>
class HalfCoreTest : public ::testing::Test
{
protected:
  
  virtual void SetUp()
  {
    cell_id_limit = _Traits::HalfT::cell_id_limit;
    position_limit= _Traits::HalfT::position_limit;
    anchor_limit  = _Traits::HalfT:: anchor_limit;
  }

  virtual void TearDown()
  {
  }
  
  int cell_id_limit;
  int position_limit;
  int anchor_limit;
  
  typename _Traits::HalfT  half;
};


typedef ::testing::Types<Tri2dTraits, Tri3dTraits, Tet3dTraits> TraitsTypes;
TYPED_TEST_CASE(HalfCoreTest, TraitsTypes);

TYPED_TEST(HalfCoreTest, Constructors) {
  
  typedef typename TypeParam::HalfT HalfT;
  
  HalfT H[] = {HalfT(0, 0, 0), HalfT(0, -1, 0),
               HalfT(this->cell_id_limit, this->position_limit, this->anchor_limit),
               HalfT(this->cell_id_limit, -1, this->anchor_limit),
               HalfT(2, 2, 2)};

  EXPECT_EQ(0, H[0].getIncidCell());
  EXPECT_EQ(0, H[0].getPosition());
  EXPECT_EQ(0, H[0].getAnchor());
  
  EXPECT_EQ( 0, H[1].getIncidCell());
  EXPECT_EQ(-1, H[1].getPosition());
  EXPECT_EQ( 0, H[1].getAnchor());
  
  EXPECT_EQ( (int)this->cell_id_limit, H[2].getIncidCell());
  EXPECT_EQ( this->position_limit, H[2].getPosition());
  EXPECT_EQ( (int)this->anchor_limit, H[2].getAnchor());
  
  EXPECT_EQ( (int)this->cell_id_limit, H[3].getIncidCell());
  EXPECT_EQ( -1, H[3].getPosition());
  EXPECT_EQ( (int)this->anchor_limit, H[3].getAnchor());
  
  EXPECT_EQ( 2, H[4].getIncidCell());
  EXPECT_EQ( 2, H[4].getPosition());
  EXPECT_EQ( 2, H[4].getAnchor());
  
}


TYPED_TEST(HalfCoreTest, Sets) {
  
  typedef typename TypeParam::HalfT HalfT;
  
  this->half.setIncidCell(0);
  this->half.setPosition(-1);
  this->half.setAnchor(0);
  
  EXPECT_EQ( 0, this->half.getIncidCell());
  EXPECT_EQ( -1, this->half.getPosition());
  EXPECT_EQ( 0, this->half.getAnchor());
  
  this->half.setIncidCell(this->cell_id_limit);
  this->half.setPosition(this->position_limit);
  this->half.setAnchor(this->anchor_limit);
  
  EXPECT_EQ( (int)this->cell_id_limit, this->half.getIncidCell());
  EXPECT_EQ( this->position_limit,      this->half.getPosition());
  EXPECT_EQ( (int)this->anchor_limit,  this->half.getAnchor());

  this->half.setIncidCell(2);
  this->half.setPosition(2);
  this->half.setAnchor(2);
  
  EXPECT_EQ( 2, this->half.getIncidCell());
  EXPECT_EQ( 2,  this->half.getPosition());
  EXPECT_EQ( 2, this->half.getAnchor());

  /* Repeat with getCompleteId */
  
  this->half.setCompleteId(0,-1,0);
  EXPECT_EQ( 0, this->half.getIncidCell());
  EXPECT_EQ( -1, this->half.getPosition());
  EXPECT_EQ( 0, this->half.getAnchor());

  this->half.setCompleteId(this->cell_id_limit,this->position_limit,this->anchor_limit);
  EXPECT_EQ( (int)this->cell_id_limit, this->half.getIncidCell());
  EXPECT_EQ( this->position_limit,      this->half.getPosition());
  EXPECT_EQ( (int)this->anchor_limit,  this->half.getAnchor());
  
  this->half.setCompleteId(2,2,2);
  EXPECT_EQ( 2, this->half.getIncidCell());
  EXPECT_EQ( 2,  this->half.getPosition());
  EXPECT_EQ( 2, this->half.getAnchor());
  
}



