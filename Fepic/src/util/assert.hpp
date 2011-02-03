// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.

#ifndef FEPIC_ASSERT_HPP
#define FEPIC_ASSERT_HPP

void assertion_failed(char const* expr, char const* function, char const* file, long line, const char* msg);
bool verify_failed(char const* expr, char const* function, char const* file, long line, const char* msg);




#ifdef FEPIC_DEBUG_ON
  #define FEPIC_ASSERT(ok, msg) if(!(ok)) \
                                  assertion_failed(#ok, __PRETTY_FUNCTION__, __FILE__, __LINE__, msg)
#else
  #define FEPIC_ASSERT(ok,msg) ((void)0)
#endif




#define FEPIC_VERIFY(ok, msg) ((ok)? false : verify_failed(#ok, __PRETTY_FUNCTION__, __FILE__, __LINE__, msg))


void assertion_failed(char const* expr, char const* function, char const* file, long line, const char* msg)
{
  std::cout << "ERROR: " <<file<<":"<<line<<": "<<msg<<"\n"
            << "ERROR: " <<file<<":"<<line<<": assertion '"<<expr<<"' failed\n"
            << "ERROR: " <<file<<": in "<<"'"<<function<<"'"<<std::endl;
  throw;
}

bool verify_failed(char const* expr, char const* function, char const* file, long line, const char* msg)
{
  std::cout << "WARNING: " <<file<<":"<<line<<": "<<msg<<"\n"
            << "WARNING: " <<file<<":"<<line<<": assertion '"<<expr<<"' failed\n"
            << "WARNING: " <<file<<": in "<<"'"<<function<<"'"<<std::endl;
  return true;
}

#endif
