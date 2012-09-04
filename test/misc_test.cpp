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
#include <vector>
#include <deque>
#include <algorithm>
#include <string>
#include <fstream>
#include <Fepic/Mesh>

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::ValuesIn;
using ::testing::Combine;
using ::std::tr1::tuple;
using ::std::tr1::get;



// In this sample, we want to ensure that every test finishes within
// ~5 seconds.  If a test takes longer to run, we consider it a
// failure.
//
// We put the code for timing a test in a test fixture called
// "QuickTest".  QuickTest is intended to be the super fixture that
// other fixtures derive from, therefore there is no test case with
// the name "QuickTest".  This is OK.
//
// Later, we will derive multiple test fixtures from QuickTest.
class QuickTest : public testing::Test {
 protected:
  // Remember that SetUp() is run immediately before a test starts.
  // This is a good place to record the start time.
  virtual void SetUp() {
    start_time_ = time(NULL);
  }

  // TearDown() is invoked immediately after a test finishes.  Here we
  // check if the test was too slow.
  virtual void TearDown() {
    // Gets the time when the test finishes
    const time_t end_time = time(NULL);

    // Asserts that the test took no more than ~5 seconds.  Did you
    // know that you can use assertions in SetUp() and TearDown() as
    // well?
    EXPECT_TRUE(end_time - start_time_ <= 3) << "The test took too long.";
  }

  // The UTC time (in seconds) when the test starts
  time_t start_time_;
};



/*  getExtension Test 
 * 
 *  Nomenclatures:
 * 
 *  WE := With Extension
 *  NE := No Extension
 */

class getExtension1Test : public TestWithParam<tuple<const char*, const char*> > {
protected:
  virtual void SetUp() {
    filename = std::string(get<0>(GetParam())) + std::string(get<1>(GetParam()));
  }
  std::string filename;
};

class getExtension2Test : public TestWithParam<tuple<const char*, const char*> > {
protected:
  virtual void SetUp() {
    filename = std::string(get<0>(GetParam())) + std::string(get<1>(GetParam()));
  }
  std::string filename;
};

const char* directories[] = {"", "./", "bar/", "./bar/", "../", "./../", 
                             "./../bar/", "/bar/", "bar.doc/", "./bar.doc/",
                             "/bar.doc/", "./../bar.doc/"}; // 12

const char* file_yes_ext[] = {"foo.dat", "moe.foo.dat", "moe foo.dat",
                              "moe-foo.dat", ".foo.dat", ".foo.bar.dat"};
                              
const char* file_no_ext[] = {"foo", "moe foo","moe-foo", ".foo", "foo.", "foo..."}; //5
                              

TEST_P(getExtension1Test, FilesWithExtension) {
  /* this folder */
  EXPECT_EQ(".dat", getExtension(filename));
}


TEST_P(getExtension2Test, FilesWithoutExtension) {
  /* this folder */
  EXPECT_EQ("", getExtension(filename));
}


INSTANTIATE_TEST_CASE_P(Misc,
                        getExtension1Test,
                        Combine(ValuesIn(directories), ValuesIn(file_yes_ext)));

INSTANTIATE_TEST_CASE_P(Misc,
                        getExtension2Test,
                        Combine(ValuesIn(directories), ValuesIn(file_no_ext)));

TEST(getRelativePathTest, UsualCases) {
  EXPECT_EQ("./", getRelativePath("foo.dat"));
  EXPECT_EQ("./", getRelativePath("./foo.dat"));
  EXPECT_EQ("../", getRelativePath("../foo.dat"));
  EXPECT_EQ("./../", getRelativePath("./../foo.dat"));
  EXPECT_EQ("../bar/", getRelativePath("../bar/foo.dat"));
  EXPECT_EQ("../bar/jow/moe/", getRelativePath("../bar/jow/moe/foo.dat"));
}

TEST(getBaseNameTest, UsualCases) {
  EXPECT_EQ("foo", getBaseName("foo.dat"));
  EXPECT_EQ("foo", getBaseName("./foo.dat"));
  EXPECT_EQ("foo", getBaseName("../foo.dat"));
  EXPECT_EQ("foo", getBaseName("./../foo.dat"));
  EXPECT_EQ("foo", getBaseName("../bar/foo.dat"));
  EXPECT_EQ("foo", getBaseName("foo"));
  EXPECT_EQ("foo", getBaseName("./foo"));
  EXPECT_EQ("foo", getBaseName("../foo"));
  EXPECT_EQ("foo", getBaseName("./../foo"));
  EXPECT_EQ("foo", getBaseName("../bar/foo"));
  EXPECT_EQ(".foo", getBaseName(".foo"));
  EXPECT_EQ(".foo", getBaseName("./.foo"));
  EXPECT_EQ(".foo", getBaseName("../.foo"));
  EXPECT_EQ(".foo", getBaseName("./../.foo"));
  EXPECT_EQ(".foo", getBaseName("../bar/.foo"));  
}

TEST(stripTrailingSpacesTest, UsualCases) {
  EXPECT_EQ("", stripTrailingSpaces(""));
  EXPECT_EQ("", stripTrailingSpaces(" "));
  EXPECT_EQ("", stripTrailingSpaces("       "));
  EXPECT_EQ("f", stripTrailingSpaces("f"));
  EXPECT_EQ("f", stripTrailingSpaces("f "));
  EXPECT_EQ("f", stripTrailingSpaces("f     "));
}

TEST(itoaTest, UsualCases) {
  EXPECT_EQ("-1", itoa(-1));
  EXPECT_EQ("-123", itoa(-123));
  EXPECT_EQ("123", itoa(123));
  EXPECT_EQ("0", itoa(0));
  EXPECT_EQ("0", itoa(-0));
  EXPECT_EQ("123", itoa(123u));
  EXPECT_EQ("0", itoa(0u));
}

//class itoafill0Test : public QuickTest {};

TEST(itoafill0Test, UsualCases) {
  EXPECT_EQ("1", itoafill0(1,0));
  EXPECT_EQ("1", itoafill0(1,1));
  EXPECT_EQ("01", itoafill0(1,2));
  EXPECT_EQ("001", itoafill0(1,3));
  EXPECT_EQ("0001", itoafill0(1,4));
  EXPECT_EQ("00012", itoafill0(12,5));
  EXPECT_EQ("001001", itoafill0(1001,6));
  EXPECT_EQ("0909090", itoafill0(909090,7));
  EXPECT_EQ("909090", itoafill0(909090,3));
}



// array test

TEST(arrayIsCyclicallyEqualTest, ExpectTrueTest) {
  
  std::vector<int> U, V;
  int size_max = 4;
  U.reserve(size_max); V.reserve(size_max);
  
  for (int size = 1; size <= size_max; ++size)
  {
    U.resize(size); V.resize(size);
    for (int i = 0; i < size; ++i)
    {
      U.at(i) = i+1;
      V.at(i) = U.at(i);
    }
    
    // cyclically
    for (int k = 0; k < size; ++k)
    {
      EXPECT_EQ(true, arrayIsCyclicallyEqual(U.begin(),U.end(),V.begin(),V.end()));
      std::rotate(V.begin(), V.begin()+1, V.end());
    }
    
    std::reverse(V.begin(), V.end());
    
    // anti-cyclically
    for (int k = 0; k < size; ++k)
    {
      EXPECT_EQ(true, arrayIsCyclicallyEqual(U.begin(),U.end(),V.begin(),V.end()));
      std::rotate(V.begin(), V.begin()+1, V.end());
    }
      
  }
  
}

TEST(arrayIsCyclicallyEqualTest, ExpectFalseTest) {
  
  std::vector<int> U, V;
  int size_max = 4;
  U.reserve(size_max); V.reserve(size_max);
  
  for (int size = 1; size <= size_max; ++size)
  {
    U.resize(size); V.resize(size);
    for (int i = 0; i < size; ++i)
    {
      U.at(i) = i+1;
      V.at(i) = U.at(i);
    }
    V.at(0) += 1;
    
    // cyclically
    for (int k = 0; k < size; ++k)
    {
      EXPECT_EQ(0, arrayIsCyclicallyEqual(U.begin(),U.end(),V.begin(),V.end()));
      std::rotate(V.begin(), V.begin()+1, V.end());
    }
    
    std::reverse(V.begin(), V.end());
    
    // anti-cyclically
    for (int k = 0; k < size; ++k)
    {
      EXPECT_EQ(0, arrayIsCyclicallyEqual(U.begin(),U.end(),V.begin(),V.end()));
      std::rotate(V.begin(), V.begin()+1, V.end());
    }
      
  }
  
}


TEST(sameElementsTest, OutputTest) {
  
  int U[] = {1,2,3,4,5}, V[] = {1,2,3,4,5};
  int size = sizeof(U)/sizeof(int);
  
  do {
    EXPECT_TRUE(sameElements(U,U+size,V,V+size));
  } while ( std::next_permutation (V,V+size) );  
  
  V[0] = 99;
  
  do {
    EXPECT_FALSE(sameElements(U,U+size,V,V+size));
  } while ( std::next_permutation (V,V+size) );  
  
}




