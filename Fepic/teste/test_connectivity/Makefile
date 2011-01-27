# Um exemplo de makefile:
# Neste exemplo, são simplesmente compilados alguns sources
# ----------------------------------------------#

CXX	    = g++
#CFLAGS	= -Wall -Winline -g -DEIGEN_NO_DEBUG
CFLAGS	= -std=c++0x -Wall -Winline -g -pg
LD	    = $(CXX)
LDFLAGS	=
RM	    = rm
INCPATH = -I$(EIGEN_DIR) -I$(EFEM_DIR)
#LIBS	= -Wl,-rpath,$(EFEM_DIR)/Fepic/mesh -L$(EFEM_DIR)/Fepic/mesh -lmesh

EXE	    = teste
SRCS	= teste.cpp

# pega SRCS e substitui .c por .o
OBJS	= ${SRCS:.cpp=.o}


# define a suffix rule for .c -> .o
$(OBJS) :
	$(CXX) $(INCPATH) $(CFLAGS) -c $(SRCS)

# default target by convention is ``all''
all : $(OBJS) $(EXE) clean

$(EXE) : $(OBJS)
	$(LD) -o $@ $(OBJS)

clean :
	-$(RM) -f *.o
	
run:
	valgrind --leak-check=yes ./$(EXE)

# No caso geral, $@ se refere ao alvo





# ----------------------------------------------#


# exemplo de compilar uma biblioteca estática

