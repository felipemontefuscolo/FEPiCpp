# Um exemplo de makefile:
# Neste exemplo, são simplesmente compilados alguns sources
# ----------------------------------------------#
# lembre-se que a opção "-pg" deve ser colocado na hora de linkar e compilar

include $(FEPIC_DIR)/conf/variables

CXX	    = g++
#CFLAGS	= -Wall -Winline -g -DEIGEN_NO_DEBUG -O3 -msse2
CFLAGS	= -Wall -Wextra -g
LD	    = $(CXX)
LDFLAGS	= -L${FEP_LIBS_DIR} -lfepic -fopenmp
#LDFLAGS	=
RM	    = rm
INCPATH = -I$(EIGEN_DIR) -I$(FEPIC_DIR)

EXE	  = test.exe
SRCS	= test.cpp

# pega SRCS e substitui .c por .o
OBJS	= ${SRCS:.cpp=.o}


# default target by convention is ``all''
all : $(OBJS) $(EXE) clean

# define a suffix rule for .c -> .o
$(OBJS) :
	$(CXX) $(INCPATH) $(CFLAGS) -c $(SRCS)

# colocar -pg aqui também
$(EXE) : $(OBJS)
	$(LD) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean :
	-$(RM) -f *.o
	
run:
	valgrind --leak-check=yes ./$(EXE)

# No caso geral, $@ se refere ao alvo





# ----------------------------------------------#


# exemplo de compilar uma biblioteca estática

