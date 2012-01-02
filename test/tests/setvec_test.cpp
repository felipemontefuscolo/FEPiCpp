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


#include <contrib/Loki/set_vector.hpp>
#include <iostream>
#include <sstream>
#include <contrib/Loki/AssocVector.hpp>
#include <utility>

template<class V>
std::string print_v(V const& v) {
  std::stringstream ss;
  
  for (typename V::const_iterator i=v.begin(); i < v.end(); ++i)
    ss << (*i) << " ";
  ss << "\n";
  
  return ss.str();
}

TEST(SetVectorTest, GeneralTest) 
{
  int a[] = {4,1,2,3,5,3,1};
  SetVector<int> v0(a, a+7); // 1,2,3,4,5
  SetVector<int> v1;
  
  v1.reserve(5);
  
  for (int i = 0; i < 5; ++i)
    v1.insert(i+1);
  
  EXPECT_TRUE(v0 == v1) << print_v(v0);
  
  v0.insert(8);
  v0.insert(7);
  v0.insert(6);
  
  int b[] = {1,2,3,4,5,6,7,8};
  EXPECT_TRUE(v0 == SetVector<int>(b, b+8)) << print_v(v0);
  
  v0.pop_back();
  
  EXPECT_TRUE(v0 == SetVector<int>(b, b+7)) << print_v(v0);
  EXPECT_EQ(v0.back() , b[6]);
  
  
};
