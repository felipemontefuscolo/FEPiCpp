# Um exemplo de makefile:
# Neste exemplo, são simplesmente compilados alguns sources
# ----------------------------------------------#
# lembre-se que a opção "-pg" deve ser colocado na hora de linkar e compilar

CXX	    = g++-4.5
#CFLAGS	= -Wall -Winline -g -DEIGEN_NO_DEBUG
CFLAGS	= -std=c++0x -Wall -g
LD	    = $(CXX)
LDFLAGS	=
RM	    = rm
INCPATH = -I $(EIGEN_DIR) -I$(FEPIC_DIR)
#LIBS	= -Wl,-rpath,$(EFEM_DIR)/Fepic/mesh -L$(EFEM_DIR)/Fepic/mesh -lmesh

EXE	    = teste.exe
SRCS	= teste.cpp

# pega SRCS e substitui .c por .o
OBJS	= ${SRCS:.cpp=.o}


# default target by convention is ``all''
all : $(OBJS) $(EXE) clean

# define a suffix rule for .c -> .o
$(OBJS) :
	$(CXX) $(INCPATH) $(CFLAGS) -c $(SRCS)

# colocar -pg aqui também
$(EXE) : $(OBJS)
	$(LD) $(CFLAGS) -o $@ $(OBJS)

clean :
	-$(RM) -f *.o
	
run:
	valgrind --leak-check=yes ./$(EXE)

# No caso geral, $@ se refere ao alvo





# ----------------------------------------------#


# exemplo de compilar uma biblioteca estática

