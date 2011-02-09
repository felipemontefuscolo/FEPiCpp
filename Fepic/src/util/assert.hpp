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

#ifndef FEPIC_CHECK_HPP
#define FEPIC_CHECK_HPP

template<class Exception>
void assertion_failed(std::string const& expr, std::string const& msg);


// see boost
#ifdef FEPIC_DEBUG_ON
  #define FEPIC_CHECK(ok, msg, except) if(!(ok)) assertion_failed<except>(#ok, msg)
#else
  #define FEPIC_CHECK(ok,msg, except) ((void)0)
#endif

#define FEPIC_ASSERT(ok, msg, except) if(!(ok)) assertion_failed<except>(#ok, msg)

template<class Exception>
void assertion_failed(std::string const& expr, std::string const& msg)
{
  std::string what_arg ="\nERROR: "__FILE__":"+std::string(itoa(__LINE__))+": "+msg+"\n"
                          "ERROR: "__FILE__":"+std::string(itoa(__LINE__))+": assertion '"+expr+"' failed\n"
                          "ERROR: "__FILE__": in '"+std::string(__PRETTY_FUNCTION__)+"' \n";
  throw Exception(what_arg);
}


#endif
