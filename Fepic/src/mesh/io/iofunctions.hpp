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

#ifndef FEPIC_IOFUNCTIONS_HPP
#define FEPIC_IOFUNCTIONS_HPP

template<class T>
static void printData(T const& data, int n_data, std::string const& filename)
{
  std::ofstream Fout(filename.data());
  for (int i = 0; i < n_data; i++)
  {
    Fout << data[i] << std::endl;
  }
  Fout.close();
}

template<class T>
static void readData(T & data, int n_data, std::string const& filename)
{
  std::ifstream Fin(filename.data());
  FEPIC_ASSERT(!Fin.fail(), "can not find mesh file", std::invalid_argument);
  
  for (int i = 0; i < n_data; i++)
  {
    Fin >> data[i];
  }
  Fin.close();
}


#endif

