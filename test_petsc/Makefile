EXE			+= teste
SOURCES		+= teste.cpp
OBJECTS		+= teste.o

#onde é procurado as bibliotecas
#SLINKER		+= -Wl,-rpath,/home/felipe/Slibs/smesh -L/home/felipe/Slibs/smesh

#linka a biblioteca
#SLIBS		+= -lsmesh

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PETSC_INCLUDE += -I${EIGEN_DIR}

teste: teste.o chkopts
	-${CLINKER} -o $(EXE) $(OBJECTS) ${PETSC_KSP_LIB}
	${RM} *.o

include ${PETSC_DIR}/conf/test
