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
#include <cstdlib>

/* necessary inheritance: can not instantiate _Labelable because it has
 * protected constructor */
class tLabelable : public _Labelable
{
public:
  tLabelable(int tag, int flags=0) : _Labelable(tag, flags)
  {
  }

  tLabelable() : _Labelable(0, 0) {}
};

class LabelableTest : public testing::Test {
protected:

  virtual void SetUp() {
  }

};

#ifdef FEPIC_DEBUG_ON
TEST(LabelableTest, ConstructorErrors) {
  ASSERT_ANY_THROW(tLabelable(256));
  ASSERT_ANY_THROW(tLabelable(-1));
}
#endif


TEST(LabelableTest, Constructors) { // also copy assignment

  srand(123);

  tLabelable L[9];
                    
  int tags[8], idx;
  
  // create
  for (int b0 = 0; b0 < 2; ++b0)
    for (int b1 = 0; b1 < 2; ++b1)
      for (int b2 = 0; b2 < 2; ++b2)
      {
        idx = 4*b0 + 2*b1 + b2;
        tags[idx] = rand()%255;
        L[idx] = tLabelable(tags[idx], b0*DISABLED | b1*MARKED | b2*VISITED);
      }

  //check
  for (int b0 = 0; b0 < 2; ++b0)
    for (int b1 = 0; b1 < 2; ++b1)
      for (int b2 = 0; b2 < 2; ++b2)
      {
        idx = 4*b0 + 2*b1 + b2;
        EXPECT_EQ(tags[idx], L[idx].getTag());
        EXPECT_EQ(b0, L[idx].disabled());
        EXPECT_EQ(b1, L[idx].marked());
        EXPECT_EQ(b2, L[idx].visited());
      }

  // check deafult constructor
  EXPECT_EQ(0, L[8].getTag());
  EXPECT_EQ(0, L[8].marked());
  EXPECT_EQ(0, L[8].disabled());
  EXPECT_EQ(0, L[8].visited());

}


TEST(LabelableTest, SetUps) {

  srand(123);

  tLabelable L[8];
                    
  int tags[8], idx;

  // setups
  for (int b0 = 0; b0 < 2; ++b0)
    for (int b1 = 0; b1 < 2; ++b1)
      for (int b2 = 0; b2 < 2; ++b2)
      {
        idx = 4*b0 + 2*b1 + b2;
        tags[idx] = rand()%255;
        L[idx].setTag(tags[idx]);
        L[idx].disabled(b0);
        L[idx].marked(b1);
        L[idx].visited(b2);
      }


  //check
  for (int b0 = 0; b0 < 2; ++b0)
    for (int b1 = 0; b1 < 2; ++b1)
      for (int b2 = 0; b2 < 2; ++b2)
      {
        idx = 4*b0 + 2*b1 + b2;
        EXPECT_EQ(tags[idx], L[idx].getTag());
        EXPECT_EQ(b0, L[idx].disabled());
        EXPECT_EQ(b1, L[idx].marked());
        EXPECT_EQ(b2, L[idx].visited());
      }



}












