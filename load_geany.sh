# Carrega os arquivos do projeto no geany
# inicie o geany antes de executar esse arquivo
#                           remove arquivos ocultos e coloca o path absoluto
A=`find . -iname *.hpp    | sed '/\/\./d' | sed 's=\.\/='"$PWD"'\/=g'` 
B=`find . -iname *.cpp    | sed '/\/\./d' | sed 's=\.\/='"$PWD"'\/=g'` 
C=`find . -iname *.rst    | sed '/\/\./d' | sed 's=\.\/='"$PWD"'\/=g'` 
D=`find . -iname *.txt    | sed '/\/\./d' | sed 's=\.\/='"$PWD"'\/=g'` 
E=`find . -iname makefile | sed '/\/\./d' | sed 's=\.\/='"$PWD"'\/=g'` 
ALLFILES="$A $B $C $D $E"
geany $ALLFILES

