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

/* necessary inheritance: can not instantiate _Labelable because it has
 * protected constructor */
class tLabelable : public _Labelable
{
public:
  template<class... Args>
  tLabelable(Args... args) : _Labelable(args...) {}
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


TEST(LabelableTest, Constructors) {

  tLabelable L[] = {tLabelable(),
                    tLabelable(0),
                    tLabelable(100,DISABLED),
                    tLabelable(123,MARKED),
                    tLabelable(255,DISABLED | MARKED) };

  EXPECT_EQ(0, L[0].marked());
  EXPECT_EQ(0, L[0].disabled());
  EXPECT_EQ(0, L[0].getTag());

  EXPECT_EQ(0, L[1].marked());
  EXPECT_EQ(0, L[1].disabled());
  EXPECT_EQ(0, L[1].getTag());

  EXPECT_EQ(0, L[2].marked());
  EXPECT_EQ(1, L[2].disabled());
  EXPECT_EQ(100, L[2].getTag());

  EXPECT_EQ(1, L[3].marked());
  EXPECT_EQ(0, L[3].disabled());
  EXPECT_EQ(123, L[3].getTag());

  EXPECT_EQ(1, L[4].marked());
  EXPECT_EQ(1, L[4].disabled());
  EXPECT_EQ(255, L[4].getTag());

}


TEST(LabelableTest, SetUps) {

  tLabelable L[8];
  L[0].marked(0); L[0].disabled(0); L[0].setTag(0);
  L[1].marked(0); L[1].disabled(0); L[1].setTag(1);
  L[2].marked(0); L[2].disabled(1); L[2].setTag(0);
  L[3].marked(0); L[3].disabled(1); L[3].setTag(255);
  L[4].marked(1); L[4].disabled(0); L[4].setTag(0);
  L[5].marked(1); L[5].disabled(0); L[5].setTag(255);
  L[6].marked(1); L[6].disabled(1); L[6].setTag(0);
  L[7].marked(1); L[7].disabled(1); L[7].setTag(255);

  EXPECT_EQ(0, L[0].marked());
  EXPECT_EQ(0, L[0].disabled());
  EXPECT_EQ(0, L[0].getTag());

  EXPECT_EQ(0, L[1].marked());
  EXPECT_EQ(0, L[1].disabled());
  EXPECT_EQ(1, L[1].getTag());

  EXPECT_EQ(0,   L[2].marked());
  EXPECT_EQ(1,   L[2].disabled());
  EXPECT_EQ(0,   L[2].getTag());

  EXPECT_EQ(0,   L[3].marked());
  EXPECT_EQ(1,   L[3].disabled());
  EXPECT_EQ(255, L[3].getTag());

  EXPECT_EQ(1, L[4].marked());
  EXPECT_EQ(0, L[4].disabled());
  EXPECT_EQ(0, L[4].getTag());

  EXPECT_EQ(1, L[5].marked());
  EXPECT_EQ(0, L[5].disabled());
  EXPECT_EQ(255, L[5].getTag());

  EXPECT_EQ(1, L[6].marked());
  EXPECT_EQ(1, L[6].disabled());
  EXPECT_EQ(0, L[6].getTag());

  EXPECT_EQ(1,   L[7].marked());
  EXPECT_EQ(1,   L[7].disabled());
  EXPECT_EQ(255, L[7].getTag());

}












