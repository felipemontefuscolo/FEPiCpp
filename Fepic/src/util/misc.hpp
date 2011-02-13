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

#ifndef FEPIC_MISC_HPP
#define FEPIC_MISC_HPP

/* requirements
 * <vector>
 * <deque>
 * <string>
 * <fstream>
 * 
 */ 

//namespace Fepic {
  
  typedef std::vector<std::vector<int> > matrixi;

  typedef std::vector<int> vectori;
  typedef std::deque<int>  dequei;

  /* vetores para representar os espaços R² e R³ são os do Eigen*/

namespace Fepic {
  template <class T>
  const T& max ( const T& a) {
    return a;
  }

  template <class T, class ... O>
  const T& max ( const T& a, const O& ... b ) {
    auto& c(max(b...));
    return (c<a)?a:c;
  }

}


/** função parecida com a copy do STL. \n
 * copia o conteúdo de uma sequência `a' para outra `c' no intervalo [c_begin, c_end).
 * @param a_begin iterador na posição inicial da sequência `a'
 * @param c_begin iterador na posição inicial da sequência `c' (c de cópia)
 * @param c_end iterador na posição final da seguência (excluindo o último elemento)
 */
template<class In_it, class Out_it>
inline
void copy_from(In_it a_begin, Out_it c_begin, Out_it c_end)
{
  while(c_begin != c_end) *c_begin++ = *a_begin++;
}


/** Verifica se dois vetores são (anti)ciclicamente iguais.
 *  TESTAR, POIS FIZ MUDANÇAS
 */
template<class T>
bool arrayIsCyclicallyEqual(std::vector<T> const& v1, std::vector<T> const& v2)
{
  std::vector<T> vaux(v2);
  int vsize = v2.size();
  auto vaux_begin = vaux.begin(), vaux_end = vaux.end();

  for (int i = 0; i < vsize; ++i)
  {
    if (v1 == vaux)
    {
      return true;
    }
    std::rotate(vaux_begin, vaux_begin+1, vaux_end);
  }

  std::reverse(vaux_begin, vaux_end);
  if(vaux == v1) return true;

  for (int i = 1; i < vsize; ++i)
  {
    std::rotate(vaux_begin, vaux_begin+1, vaux_end);
    if (v1 == vaux)
    {
      return true;
    }
  }

  return false;
}


template<class T>
inline
void printArray(T const& v, std::ostream &o, int N)
{
  o << v[0];
  for (int i=1; i<N; ++i)
  {
    o << " " << v[i];
  }
}

template<class T>
inline
void printMatrix(T const& mat, std::ostream &o, int N, int M)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < M; j++)
    {
      o << mat[i][j] << " ";
    }
    o << std::endl;
  }
}


/** Retorna o n-ésimo bit da constante c passada.
 */
template<class T>
inline bool get_bit(T const& c, unsigned bitoffset)
{
    T bitmask = ( 1 << bitoffset) ;
    return ((c & bitmask)!=0) ? 1 : 0;
}

/** Seta o n-ésimo bit da constante c passada.
 */
template<class T>
inline void set_bit(T & c, unsigned bitoffset)
{
    T bitmask = 1 << bitoffset;
    c = (c | bitmask);
}

/** Zera o n-ésimo bit da constante c passada.
 */
template<class T>
inline void unset_bit(T & c, unsigned bitoffset)
{
    T bitmask = 1 << bitoffset;
    c = (c & (~bitmask));
}


// files

/* procura uma palavra e deixa/retorna o indicador de leitura nesta posição */
inline
size_t search_word(std::ifstream &File, const char* word) {
  std::string s;
  std::string word_s = std::string(word);

    do
  {
    s="";
    File >> s;
    if (File.peek()==EOF)
    {
      return static_cast<size_t>(-1);
    }
  } while(s!=word_s);

  return size_t(File.tellg());

}

/**
 * C++ version 0.4 std::string style "itoa":
 * Contributions from Stuart Lowe, Ray-Yuan Sheu,
 * Rodrigo de Salvo Braz, Luc Gallant, John Maloney
 * and Brian Hunt
 * 
 * modified: Felipe Montefuscolo
 */
std::string itoa(int value)
{

	std::string buf;

	enum { kMaxDigits = 35, base=10 };
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.

	int quotient = value;

	// Translating number to string with base:
	do {
		buf += "0123456789"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );

	// Append the negative sign
	if ( value < 0) buf += '-';

	std::reverse( buf.begin(), buf.end() );
	return buf;
}

std::string itoafill0(int value, int fill)
{
	std::string buf;

	enum { kMaxDigits = 35, base=10 };
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.

  fill = fill < kMaxDigits ? fill : kMaxDigits;

	int quotient = value;

	// Translating number to string with base:
	do {
		buf += "0123456789"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
  
  int i(fill - buf.size());
  
  while(i>0)
  {
    buf += '0';
    i=fill - buf.size();
  }
  
  std::reverse( buf.begin(), buf.end() );
  return buf;
}

// remove espaços finais de strings
inline
std::string stripTrailingSpaces(std::string name)
{
  if (name.empty())
    return name;
  while( *(name.end()-1) == ' ')
    name.erase(name.end()-1);
  return name;
}

/** @example <tt>user/file.dat</tt> returns <tt>user/</tt>
 *  @note aassumes that is a valid file name.
 */ 
inline
std::string getRelativePath(std::string const& name)
{
  size_t bar = name.rfind("/");
  
  if (bar==std::string::npos)
    return "./";
  else
    return name.substr(0,bar+1);
}

/** @example <tt>user/file.dat</tt> returns <tt>.dat</tt>
 *  @note assumes that is a valid file name.
 */ 
inline
std::string getExtension(std::string const& name)
{
  int size = name.size();
  if (size==0)
    return name;
  
  size_t dot = name.rfind(".");
  size_t bar = name.rfind("/");
  if (dot==std::string::npos)
    return "";
  else if ((bar==dot-1) || (name[dot-1]=='.') || (dot==0) || 
           ((bar!=std::string::npos) && (bar>dot)) || (dot==name.size()-1) )
    return "";
  else 
    return name.substr(dot);
}

/** @example <tt>user/file.dat</tt> returns <tt>file</tt>
 *  @note assumes that is a valid file name.
 */ 
inline
std::string getBaseName(std::string name)
{
  size_t bar = name.rfind("/");
  int npath;
  
  if (bar==std::string::npos)
    npath=0;
  else
    npath = name.substr(0,bar+1).size();

  int nexte = getExtension(name).size();
  
  for (int i = 0; i < nexte; ++i)
    name.erase(name.end()-1);
  for (int i = 0; i < npath; ++i)
    name.erase(name.begin());
  
  return name;
}


#endif
