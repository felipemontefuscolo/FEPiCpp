#include <iostream>
#include <Fepic/misc/misc.hpp>

/* 
 * 
 * Testando para as funções do arquivo Fepic/misc/misc.hpp
 * 
 * */

using namespace std;

template<class T, class F, class ... Args>
bool TEST(T valor_esperado, F func, Args ... args);

/* test functions */
int TESTarrayIsCyclicallyEqual();
int TESTmax();
int TESTbit();

int main()
{
	int cc=0;
	
	cc += TESTarrayIsCyclicallyEqual();
	cc += TESTmax();
	cc += TESTbit();
	
	cout << endl;
	cout << "TOTAL NUM OF FAILS: " << cc << endl;
	return 0;
}



/* testing tool */
template<class T>
bool TEST(T valor_esperado, T valor_obtido)
{
	if (valor_obtido == valor_esperado)
	{
		cout << "OK!" << endl;
		return 0;
	}
	else
	{
		cout << "FAIL!" << endl;
		return 1;
	}
}


int TESTarrayIsCyclicallyEqual()
{
	/* 
	 * Testando a função array arrayIsCyclicallyEqual() 
	 * com lista de inicialização;
	 * 
	 * */
	
	int cc  = 0;
	
	cout << "testando arrayIsCyclicallyEqual(): " << endl;
	
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1},         vectori{1}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2},       vectori{1,2}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2},       vectori{2,1}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{1,2,3}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{3,1,2}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{2,3,1}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{3,2,1}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{1,3,2}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3},     vectori{2,1,3}));
	cc += TEST(true, arrayIsCyclicallyEqual<int> (vectori{1,2,3,4},   vectori{1,2,3,4}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1,2,4,3}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1,4,2,3}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1,2}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1,2,3}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4},   vectori{1,2,3,4,5}));
	cc += TEST(false, arrayIsCyclicallyEqual<int>(vectori{1,2,3,4,5}, vectori{1,2,3,4}));
	cout << "number of fails: " << cc << endl;
	cout << endl;

	return cc;
	
}


int TESTmax()
{
	
	int cc = 0;
	
	cout << "testando Fepic::max(): " << endl;
	
	cc += TEST(1,  Fepic::max(1));
	cc += TEST(2,  Fepic::max(1,2));
	cc += TEST(2,  Fepic::max(2,1));
	cc += TEST(3,  Fepic::max(1,2,3));
	cc += TEST(3,  Fepic::max(1,3,2));
	cc += TEST(3,  Fepic::max(3,2,1));
	cc += TEST(1,  Fepic::max(1,1));
	cc += TEST(1,  Fepic::max(1,1,1));
	cc += TEST(1,  Fepic::max(1,1,1,1,1,1));
	cc += TEST(5,  Fepic::max(1,2,3,4,5));
	cc += TEST(10, Fepic::max(1,2,3,4,5,10,6,7,8,9,0));
	cc += TEST(-1, Fepic::max(-5,-1,-3,-2));
	cout << "number of fails: " << cc << endl;
	cout << endl;

	return cc;	
	
}


int TESTbit()
{

	cout << "testando as funções get_bit, set_bit e unset_bit: " << endl;
	
	int cc = 0;
	
	char C = 0;
	
	cc += TEST(false, get_bit(C,0) );
	set_bit(C,0);
	cc += TEST(true, get_bit(C,0) );
	unset_bit(C,0);
	cc += TEST(false, get_bit(C,0) );
	
	set_bit(C,1);
	cc += TEST(true, get_bit(C,1) );
	unset_bit(C,1);
	cc += TEST(false, get_bit(C,1) );
	
	set_bit(C,2);
	cc += TEST(true, get_bit(C,2) );
	unset_bit(C,2);
	cc += TEST(false, get_bit(C,2) );
	
	set_bit(C,3);
	cc += TEST(true, get_bit(C,3) );
	unset_bit(C,3);
	cc += TEST(false, get_bit(C,3) );
	
	set_bit(C,4);
	cc += TEST(true, get_bit(C,4) );
	unset_bit(C,4);
	cc += TEST(false, get_bit(C,4) );
	
	set_bit(C,7);
	cc += TEST(true, get_bit(C,7) );
	unset_bit(C,7);
	cc += TEST(false, get_bit(C,7) );
	


	int A=0;
	cc += TEST(false, get_bit(A,0) );
	set_bit(A,0);
	cc += TEST(true, get_bit(A,0) );
	unset_bit(A,0);
	cc += TEST(false, get_bit(A,0) );
	
	set_bit(A,1);
	cc += TEST(true, get_bit(A,1) );
	unset_bit(A,1);
	cc += TEST(false, get_bit(A,1) );
	
	set_bit(A,2);
	cc += TEST(true, get_bit(A,2) );
	unset_bit(A,2);
	cc += TEST(false, get_bit(A,2) );
	
	set_bit(A,3);
	cc += TEST(true, get_bit(A,3) );
	unset_bit(A,3);
	cc += TEST(false, get_bit(A,3) );
	
	set_bit(A,4);
	cc += TEST(true, get_bit(A,4) );
	unset_bit(A,4);
	cc += TEST(false, get_bit(A,4) );
	
	set_bit(A,7);
	cc += TEST(true, get_bit(A,7) );
	unset_bit(A,7);
	cc += TEST(false, get_bit(A,7) );	
	
	cout << "number of fails: " << cc << endl;
	cout << endl;

	return cc;	
	
}
































































