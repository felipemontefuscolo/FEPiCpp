# Um exemplo de makefile:
# Neste exemplo, são simplesmente compilados alguns sources
# ----------------------------------------------#
# lembre-se que a opção "-pg" deve ser colocado na hora de linkar e compilar

CXX	    = g++
#CFLAGS	= -O2 -march=native -msse2
CFLAGS	= -Wall -Wextra -g
LD	    = $(CXX)
LDFLAGS	= -L$(FEPIC_DIR)/Fepic/libs -lfepic
LDFLAGS+= -fopenmp
RM	    = rm
INCPATH = -I$(EIGEN_DIR) -I$(FEPIC_DIR)

EXE	  = reading.exe
SRCS	= reading.cpp

# pega SRCS e substitui .c por .o
OBJS	= ${SRCS:.cpp=.o}


# default target by convention is ``all''
all : $(OBJS) $(EXE) clean

# define a suffix rule for .c -> .o
$(OBJS) :
	$(CXX) $(INCPATH) $(CFLAGS) -c $(SRCS) $(LDFLAGS)

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

