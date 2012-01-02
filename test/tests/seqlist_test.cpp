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


#include "Fepic/src/util/list_type.hpp"
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <vector>

using namespace std;

template<class V>
std::string print_color(V const& v) {
  std::stringstream ss;

  for (typename V::const_iterator i = v.begin(); i != v.end(); ++i)
    ss << (*i).getColor() << " ";
  ss << "\n";

  return ss.str();
}

template<class V>
std::string print_link(V const& v)
{
  std::stringstream ss;

  for (typename V::const_iterator cit = v.begin(); cit != v.end(); ++cit) //
  {
    ss << (*cit).getTag() <<" | "<<(*cit).sameColorPrev() << " " << (*cit).sameColorNext() << std::endl;
  }

  return ss.str();
}

template<class V>
int color_size(V const& v, EColor c)
{
  int k=0;
  for (typename V::color_const_iterator i = v.begin(c); i!= v.end(c); ++i)
  {
    ++k;
  }
  return k;
}

template<class V>
unsigned size(V const& v)
{
  int k=0;
  for (typename V::const_iterator i = v.begin(); i!= v.end(); ++i)
    ++k;
  return k;
}

template<class V>
std::string print_vector(V const& v, int size)
{
  std::stringstream ss;
  for (int i = 0; i < size; ++i)
    ss << v[i] << " ";
  ss << "\n";
  return ss.str();
}

class Dummy : public _Colored, public _Labelable
{public:
  Dummy(EColor c=EColor(-1), int tag = 0, bool d = false) : _Colored(c), _Labelable(tag, DISABLED*d), hist(0) {}
  int hist;
};




TEST(SeqListTest, TestStep0)
{
  int a[] = {0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3}; // 6 x 4 = 24
  int a_size = sizeof(a)/sizeof(int);
  SeqList<Dummy> v0;
  SeqList<Dummy>::iterator it;
  SeqList<Dummy>::const_iterator cit;
  SeqList<Dummy>::color_iterator lit;
  SeqList<Dummy>::color_const_iterator clit;


  for (int i=0; i<a_size; ++i)
    v0.push_back(Dummy(EColor(i%4), i)); // v0 = a

  EXPECT_EQ(4, v0.numColors());

  // ================ iterators ===================

  int *p = a;
  for (cit = v0.begin(); cit != v0.end(); ++cit) //
  {
    EXPECT_EQ(*p, (*cit).getColor());
    ++p;
    //std::cout << (*cit).sameColorPrev() << " " << (*cit).sameColorNext() << std::endl;
  }

  // ================ colors iterators ===================

  for (int k = 0; k < v0.numColors(); ++k)
  {
    int t = 0;
    for (lit = v0.begin(EColor(k)); lit != v0.end(EColor(k)); ++lit) //
    {
      EXPECT_EQ(0, (*lit).getTag()%4 - k);
      t += 4;
    }
  }

  // ================== disable ===================

  // disabling colors 3, except the last one
  for (int i = 3; i <a_size-1; i+=4)
    v0.disable(i);

  // v0 = {0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,3};
  EXPECT_EQ(4, v0.numColors());
  EXPECT_EQ(19u, v0.size());

  v0.disable(a_size-1); // v0 = {0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x};

  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(18u, v0.size());


  // ================ colors iterators ===================

  // v0 = {0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x};
  for (int k = 0; k < v0.numColors(); ++k)
  {
    int t = 0;
    for (clit = v0.begin(EColor(k)); clit != v0.end(EColor(k)); ++clit) //
    {
      EXPECT_EQ(0, (*clit).getTag()%4 - k);
      t += 4;
    }
  }

  // ================== disable ===================


  // nothing
  v0.disable(a_size-1); // v0 = {0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x};

  v0.disable(1); // v0 = {0,x,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x};

  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(5, v0.colorSize(EColor(1)));
  EXPECT_EQ(5, color_size(v0,EColor(1)));
  EXPECT_EQ(17u, v0.size());

  v0.disable(21); // v0 = {0,x,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,1,2,x,  0,x,2,x};

  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(4, v0.colorSize(EColor(1)));
  EXPECT_EQ(4, color_size(v0,EColor(1)));

  v0.disable(9); // v0 = {0,x,2,x,  0,1,2,x,  0,x,2,x,  0,1,2,x,  0,1,2,x,  0,x,2,x};
  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(3, v0.colorSize(EColor(1)));
  EXPECT_EQ(3, color_size(v0,EColor(1)));


  // =============== check link =========================

  // v0 = {0,x,2,x,  0,1,2,x,  0,x,2,x,  0,1,2,x,  0,1,2,x,  0,x,2,x};
  std::vector<int> ids(3), iids(3);
  ids[0] = 5; ids[1]=13; ids[2]=17;
  clit = v0.begin(EColor(1));
  iids[0] = (*clit++).getTag();
  iids[1] = (*clit++).getTag();
  iids[2] = (*clit++).getTag();
  EXPECT_TRUE( sameElements(ids.data(), ids.data()+3, iids.data(), iids.data()+3) );

  ids.resize(6); ids[0]=0;ids[1]=4;ids[2]=8;ids[3]=12;ids[4]=16;ids[5]=20;
  iids.resize(6);
  clit = v0.begin(EColor(0));
  iids[0] = (*clit++).getTag();
  iids[1] = (*clit++).getTag();
  iids[2] = (*clit++).getTag();
  iids[3] = (*clit++).getTag();
  iids[4] = (*clit++).getTag();
  iids[5] = (*clit++).getTag();
  EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, iids.data(), iids.data()+6) );

  ids[0]=2;ids[1]=6;ids[2]=10;ids[3]=14;ids[4]=18;ids[5]=22;
  iids.resize(6);
  clit = v0.begin(EColor(2));
  iids[0] = (*clit++).getTag();
  iids[1] = (*clit++).getTag();
  iids[2] = (*clit++).getTag();
  iids[3] = (*clit++).getTag();
  iids[4] = (*clit++).getTag();
  iids[5] = (*clit++).getTag();
  EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, iids.data(), iids.data()+6) );

  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(6, v0.colorSize(EColor(0)));
  EXPECT_EQ(6, color_size(v0,EColor(0)));
  EXPECT_EQ(3, v0.colorSize(EColor(1)));
  EXPECT_EQ(3, color_size(v0,EColor(1)));
  EXPECT_EQ(6, v0.colorSize(EColor(2)));
  EXPECT_EQ(6, color_size(v0,EColor(2)));
  EXPECT_EQ(0, v0.colorSize(EColor(3)));
  EXPECT_EQ(0, color_size(v0,EColor(3)));

  // =====================================================

  // v0 = {0,x,2,x,  0,1,2,x,  0,x,2,x,  0,1,2,x,  0,1,2,x,  0,x,2,x};
  v0.disable(5);
  v0.disable(13);
  v0.disable(17);

  // v0 = {0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x};
  EXPECT_EQ(2, v0.numColors());
  EXPECT_EQ(6, v0.colorSize(EColor(0)));
  EXPECT_EQ(6, color_size(v0,EColor(0)));
  EXPECT_EQ(0, v0.colorSize(EColor(1)));
  EXPECT_EQ(0, color_size(v0,EColor(1)));
  EXPECT_EQ(6, v0.colorSize(EColor(2)));
  EXPECT_EQ(6, color_size(v0,EColor(2)));
  EXPECT_EQ(0, v0.colorSize(EColor(3)));
  EXPECT_EQ(0, color_size(v0,EColor(3)));

  // ================== link again =======================

  ids.resize(6); ids[0]=0;ids[1]=4;ids[2]=8;ids[3]=12;ids[4]=16;ids[5]=20;
  iids.resize(6);
  clit = v0.begin(EColor(0));
  iids[0] = (*clit++).getTag();
  iids[1] = (*clit++).getTag();
  iids[2] = (*clit++).getTag();
  iids[3] = (*clit++).getTag();
  iids[4] = (*clit++).getTag();
  iids[5] = (*clit++).getTag();
  EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, iids.data(), iids.data()+6) );

  ids[0]=2;ids[1]=6;ids[2]=10;ids[3]=14;ids[4]=18;ids[5]=22;
  iids.resize(6);
  clit = v0.begin(EColor(2));
  iids[0] = (*clit++).getTag();
  iids[1] = (*clit++).getTag();
  iids[2] = (*clit++).getTag();
  iids[3] = (*clit++).getTag();
  iids[4] = (*clit++).getTag();
  iids[5] = (*clit++).getTag();
  EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, iids.data(), iids.data()+6) );


  // ================== insert =======================

  // v0 = {0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x};
  int id;
  id = v0.insert(Dummy(EColor(3))); // v0 = {0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,x,  0,x,2,3};

  EXPECT_EQ(23, id);
  EXPECT_EQ(13u, v0.size());
  EXPECT_EQ(24u, v0.total_size());
  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(6, v0.colorSize(EColor(0)));
  EXPECT_EQ(6, color_size(v0,EColor(0)));
  EXPECT_EQ(0, v0.colorSize(EColor(1)));
  EXPECT_EQ(0, color_size(v0,EColor(1)));
  EXPECT_EQ(6, v0.colorSize(EColor(2)));
  EXPECT_EQ(6, color_size(v0,EColor(2)));
  EXPECT_EQ(1, v0.colorSize(EColor(3)));
  EXPECT_EQ(1, color_size(v0,EColor(3)));

  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 21, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 19, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 17, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 15, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 13, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ( 11, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ(  9, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ(  7, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ(  5, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ(  3, id);
  id = v0.insert(Dummy(EColor(3))); EXPECT_EQ(  1, id);

  // v0 = {0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3};

  EXPECT_EQ(3, v0.numColors());
  EXPECT_EQ(24u, v0.size());
  EXPECT_EQ(24u, v0.total_size());
  EXPECT_EQ(6, v0.colorSize(EColor(0)));
  EXPECT_EQ(6, color_size(v0,EColor(0)));
  EXPECT_EQ(0, v0.colorSize(EColor(1)));
  EXPECT_EQ(0, color_size(v0,EColor(1)));
  EXPECT_EQ(6, v0.colorSize(EColor(2)));
  EXPECT_EQ(6, color_size(v0,EColor(2)));
  EXPECT_EQ(12, v0.colorSize(EColor(3)));
  EXPECT_EQ(12, color_size(v0,EColor(3)));

  id = v0.insert(Dummy(EColor(5))); EXPECT_EQ( 24, id);
  id = v0.insert(Dummy(EColor(5))); EXPECT_EQ( 25, id);
  id = v0.insert(Dummy(EColor(5))); EXPECT_EQ( 26, id);

  // v0 = {0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3,  0,3,2,3, 5,5,5};
  EXPECT_EQ( 4, v0.numColors());
  EXPECT_EQ(27u, v0.size());
  EXPECT_EQ(27u, v0.total_size());
  EXPECT_EQ( 6, v0.colorSize(EColor(0)));
  EXPECT_EQ( 6, color_size(v0,EColor(0)));
  EXPECT_EQ( 0, v0.colorSize(EColor(1)));
  EXPECT_EQ( 0, color_size(v0,EColor(1)));
  EXPECT_EQ( 6, v0.colorSize(EColor(2)));
  EXPECT_EQ( 6, color_size(v0,EColor(2)));
  EXPECT_EQ(12, v0.colorSize(EColor(3)));
  EXPECT_EQ(12, color_size(v0,EColor(3)));
  EXPECT_EQ( 3, v0.colorSize(EColor(5)));
  EXPECT_EQ( 3, color_size(v0,EColor(5)));

};


TEST(SeqListTest, TestStep1)
{
  int a[] = {0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3}; // 6 x 4 = 24
  int a_size = sizeof(a)/sizeof(int);
  SeqList<Dummy> v;
  SeqList<Dummy>::iterator it;
  SeqList<Dummy>::const_iterator cit;
  SeqList<Dummy>::color_iterator lit;
  SeqList<Dummy>::color_const_iterator clit;

  v.resize(24);

  for (int i=0; i<a_size; ++i)
    v[i] = Dummy(EColor(i%4), i); // v = a

  v.linkColors();

  EXPECT_EQ(4, v.numColors());
  EXPECT_EQ(24u, v.size());
  EXPECT_EQ(24u, v.total_size());
  EXPECT_EQ( 6, v.colorSize(EColor(0)));
  EXPECT_EQ( 6, color_size(v,EColor(0)));
  EXPECT_EQ( 6, v.colorSize(EColor(1)));
  EXPECT_EQ( 6, color_size(v,EColor(1)));
  EXPECT_EQ( 6, v.colorSize(EColor(2)));
  EXPECT_EQ( 6, color_size(v,EColor(2)));
  EXPECT_EQ( 6, v.colorSize(EColor(3)));
  EXPECT_EQ( 6, color_size(v,EColor(3)));

  // check link
  std::vector<int> ids(6), eds(6);
  for (int k = 0; k < v.numColors(); ++k)
  { 
    int i;
    for (i = 0; i < 6; ++i)
      ids[i] = 4*i + k;
    lit = v.begin(EColor(k));
    i=0;
    for (; lit != v.end(EColor(k)); lit++)
      eds[i++] = (*lit).getTag();
    EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, eds.data(), eds.data()+6) )
      << "ids: " << print_vector(ids, 6) << "eds: " << print_vector(eds,6);    
  }
    

};


TEST(SeqListTest, TestStep2)
{
  int a[] = {1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7}; // 6 x 4 = 24
  //int a[] = {0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3}; // 6 x 4 = 24
  int a_size = sizeof(a)/sizeof(int);
  SeqList<Dummy> v;
  SeqList<Dummy>::iterator it;
  SeqList<Dummy>::const_iterator cit;
  SeqList<Dummy>::color_iterator lit;
  //SeqList<Dummy>::color_const_iterator clit;

  v.resize(24);

  for (int i=0; i<a_size; ++i)
    v[i] = Dummy(EColor(a[i]), i); // v = a

  v.linkColors();

  EXPECT_EQ(4, v.numColors());
  EXPECT_EQ(24u, v.size());
  EXPECT_EQ(24u, v.total_size());
  EXPECT_EQ( 0, v.colorSize(EColor(0)));
  EXPECT_EQ( 0, color_size(v,EColor(0))); // erro
  EXPECT_EQ( 6, v.colorSize(EColor(1)));
  EXPECT_EQ( 6, color_size(v,EColor(1)));
  EXPECT_EQ( 0, v.colorSize(EColor(2)));
  EXPECT_EQ( 0, color_size(v,EColor(2)));
  EXPECT_EQ( 6, v.colorSize(EColor(3)));
  EXPECT_EQ( 6, color_size(v,EColor(3)));
  EXPECT_EQ( 0, v.colorSize(EColor(4)));
  EXPECT_EQ( 0, color_size(v,EColor(4)));
  EXPECT_EQ( 6, v.colorSize(EColor(5)));
  EXPECT_EQ( 6, color_size(v,EColor(5)));
  EXPECT_EQ( 0, v.colorSize(EColor(6)));
  EXPECT_EQ( 0, color_size(v,EColor(6)));
  EXPECT_EQ( 6, v.colorSize(EColor(7)));
  EXPECT_EQ( 6, color_size(v,EColor(7)));

  // check link
  std::vector<int> ids(6), eds(6);
  for (int k = 1; k < v.maxColor()+1; k += 2)
  { 
    int i;
    for (i = 0; i < 6; ++i)
      ids[i] = 4*i + k/2;
    lit = v.begin(EColor(k));
    i=0;
    for (; lit != v.end(EColor(k)); lit++)
      eds[i++] = (*lit).getTag();
    EXPECT_TRUE( sameElements(ids.data(), ids.data()+6, eds.data(), eds.data()+6) )
      << "ids: " << print_vector(ids, 6) << "eds: " << print_vector(eds,6);    
  }
    

  v.disable(0);
  // {x,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  EXPECT_EQ(1, v.begin()->getTag());

  v.disable(0);
  // {x,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  EXPECT_EQ(1, v.begin()->getTag());

  v.disable(1);
  // {x,x,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  EXPECT_EQ(2, v.begin()->getTag());

  v.disable(2);
  // {x,x,x,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  EXPECT_EQ(3, v.begin()->getTag());

  v.disable(3);
  // {x,x,x,x,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  EXPECT_EQ(4, v.begin()->getTag());

  EXPECT_EQ(20u, v.size());
  EXPECT_EQ(24u, v.total_size());

  int id;
  id = v.insert(Dummy(EColor(7), 3)); EXPECT_EQ(3, id);
  id = v.insert(Dummy(EColor(5), 2)); EXPECT_EQ(2, id);
  id = v.insert(Dummy(EColor(3), 1)); EXPECT_EQ(1, id);
  id = v.insert(Dummy(EColor(1), 0)); EXPECT_EQ(0, id);
  // {1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7};
  
  // ================ iterators ===================

  int *p = a;
  for (cit = v.begin(); cit != v.end(); ++cit) //
  {
    EXPECT_EQ(*p, (*cit).getColor());
    ++p;
  }
  
  v.disable(0); v.disable(20);
  v.disable(1); v.disable(21);
  v.disable(2); v.disable(22);
  v.disable(3); v.disable(23);
  // {x,x,x,x,  1,3,5,7,  1,3,5,7,  1,3,5,7,  1,3,5,7,  x,x,x,x};
  
  p = a;
  for (cit = v.begin(); cit != v.end(); ++cit) //
  {
    EXPECT_EQ(*p, (*cit).getColor());
    ++p;
  }  
  p = a;
  for (it = v.begin(); it != v.end(); ++it) //
  {
    EXPECT_EQ(*p, (*it).getColor());
    ++p;
  }
  
  EXPECT_EQ(0, v.contiguousId(4));

  v.disable(5);
  v.disable(11);
  v.disable(18);
  v.disable(15);
  // {x,x,x,x,  1,x,5,7,  1,3,5,x,  1,3,5,x,  1,3,x,7,  x,x,x,x};
  
  EXPECT_EQ(12u, v.size());
  EXPECT_EQ(12u, size(v));
  EXPECT_EQ(24u, v.total_size());  
  EXPECT_EQ(0, v.contiguousId(4));  EXPECT_EQ(6, v.contiguousId(12));
  EXPECT_EQ(1, v.contiguousId(6));  EXPECT_EQ(7, v.contiguousId(13));
  EXPECT_EQ(2, v.contiguousId(7));  EXPECT_EQ(8, v.contiguousId(14));
  EXPECT_EQ(3, v.contiguousId(8));  EXPECT_EQ(9, v.contiguousId(16));
  EXPECT_EQ(4, v.contiguousId(9));  EXPECT_EQ(10, v.contiguousId(17));
  EXPECT_EQ(5, v.contiguousId(10)); EXPECT_EQ(11, v.contiguousId(19));

  int u[] = {4,6,7,8,9,10,12,13,14,16,17,19};
  int w[] = {0,1,2,3,4,5,6,7,8,9,10,11};

  v.contiguousIds(u, u+12, u);
  
  EXPECT_TRUE(sameElements(u, u+12, w, w+12));
};


TEST(SeqListTest, ParallelItersTest0)
{

  int a[] = {0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,
             0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,
             0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3,  0,1,2,3}; // 18 x 4 = 72
  int a_size = sizeof(a)/sizeof(int);
  SeqList<Dummy> v;
  SeqList<Dummy>::iterator it;
  SeqList<Dummy>::const_iterator cit;
  SeqList<Dummy>::color_iterator lit;
  SeqList<Dummy>::color_const_iterator clit;

  v.resize(72);

  for (int i=0; i<a_size; ++i)
    v[i] = Dummy(EColor(a[i]), i); // v = a

  v.linkColors();
  
  // ============================ iterators ==============================
  
  int nthreads = 3;
  int tid;
  //    simulating 3 threads
  for (tid = 0; tid<nthreads; ++tid)
    for (it = v.begin(tid, nthreads); it != v.end(tid, nthreads); ++it)
      (*it).hist += 1;
      
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0;
  }
  
  // ============================ colors iterators =========================
  
  for (int k = 0; k < v.numColors(); ++k)
  {
    //    simulating 3 threads
    for (tid = 0; tid<nthreads; ++tid)
      for (lit = v.begin(EColor(k), tid, nthreads); lit != v.end(EColor(k), tid, nthreads); ++lit)
        (*lit).hist += 1;
        
  }
  
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0;
  }
  
  // =======================================================================
  
  v.disable(0); v.disable(68);
  v.disable(1); v.disable(69);
  v.disable(2); v.disable(70);
  v.disable(3); v.disable(71);
  
  EXPECT_EQ(64u, size(v));
  EXPECT_EQ(64u, v.size());
  
  //    simulating 3 threads
  for (tid = 0; tid<nthreads; ++tid)
    for (it = v.begin(tid, nthreads); it != v.end(tid, nthreads); ++it)
      (*it).hist += 1;
      
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0; // reseting
  }
  
  // ============================ colors iterators =========================
  
  for (int k = 0; k < v.numColors(); ++k)
  {
    //    simulating 3 threads
    for (tid = 0; tid<nthreads; ++tid)
      for (lit = v.begin(EColor(k), tid, nthreads); lit != v.end(EColor(k), tid, nthreads); ++lit)
        (*lit).hist += 1;
        
  }
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0;
  }
  
};





TEST(SeqListTest, ParallelItersTest1)
{

  int a[] = {0,1,2,3,  0,1,2,3}; // 2 x 4 = 8
  int a_size = sizeof(a)/sizeof(int);
  SeqList<Dummy> v;
  SeqList<Dummy>::iterator it;
  SeqList<Dummy>::const_iterator cit;
  SeqList<Dummy>::color_iterator lit;
  SeqList<Dummy>::color_const_iterator clit;

  v.resize(8);

  for (int i=0; i<a_size; ++i)
    v[i] = Dummy(EColor(a[i]), i); // v = a

  v.linkColors();
  
  // ============================ iterators ==============================
  
  int nthreads = 10;
  int tid;
  //    simulating 10 threads
  for (tid = 0; tid<nthreads; ++tid)
    for (it = v.begin(tid, nthreads); it != v.end(tid, nthreads); ++it)
      (*it).hist += 1;
      
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0;
  }
  
  // ============================ colors iterators =========================
  
  for (int k = 0; k < v.numColors(); ++k)
  {
    //    simulating 10 threads
    for (tid = 0; tid<nthreads; ++tid)
      for (lit = v.begin(EColor(k), tid, nthreads); lit != v.end(EColor(k), tid, nthreads); ++lit)
        (*lit).hist += 1;
        
  }
  
  // checking
  for (it = v.begin(); it != v.end(); ++it)
  {
    EXPECT_EQ(1, (*it).hist);
    (*it).hist = 0;
  }
  
  // =======================================================================
  

  
};
  
