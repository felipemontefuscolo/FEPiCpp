#ifndef MISC_HPP
#define MISC_HPP

#include <algorithm>
#include <vector>
#include <array>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>

//typedef unsigned int uint;

namespace Fepic
{
	typedef std::vector<std::vector<int> > 			matrixi;
	typedef std::vector<std::vector<uint> > 		matrixUi;
	typedef std::vector<std::vector<double> > 		matrixd;
    typedef std::vector<std::vector<bool> >         matrixb;

    typedef std::vector<double> vectord;
	typedef std::vector<int> 	vectori;
    typedef std::vector<bool>   vectorb;
	typedef std::vector<uint> 	vectorui;
	typedef std::deque<uint>  	dequeui;
	
	/* vetores para representar os espaços R² e R³ são os do Eigen*/

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


/* prefixo M_ := Meta, de Metaprogramação */

/* verifica se a classe C1 é igual a C2 */
template<class C1, class C2>
class M_Compare {
public:
	static const bool isEqual = false;
};

template<class C1>
class M_Compare<C1, C1> {
public:
	static const bool isEqual = true;
};






/** Verifica se dois vetores são (anti)ciclicamente iguais.
 *  TESTAR, POIS FIZ MUDANÇAS
 */ 
template<class T>
bool arrayIsCyclicallyEqual(std::vector<T> const& v1, std::vector<T> const& v2)
{
	std::vector<T> vaux(v2);
    uint vsize = v2.size();
    auto vaux_begin = vaux.begin(), vaux_end = vaux.end();
	
	for (uint i = 0; i < vsize; ++i)
	{
		if (v1 == vaux)
		{
			return true;
		}
		rotate(vaux_begin, vaux_begin+1, vaux_end);
	}
	
	reverse(vaux_begin, vaux_end);
	if(vaux == v1) return true;
	
	for (uint i = 1; i < vsize; ++i)
	{
		rotate(vaux_begin, vaux_begin+1, vaux_end);
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
			std::cout << "Incorrect file format" << std::endl;
			return -1;
		}
	} while(s!=word_s);
	
	return size_t(File.tellg());
	
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





/* Meta: controle de fluxo if-then-else */

template<bool C, class Ta, class Tb>
class M_IfThenElse;

template<class Ta, class Tb>
class M_IfThenElse<true, Ta, Tb> {
public:
	typedef Ta resultT;
};

template<class Ta, class Tb>
class M_IfThenElse<false, Ta, Tb> {
public:
	typedef Tb resultT;
};

#endif
