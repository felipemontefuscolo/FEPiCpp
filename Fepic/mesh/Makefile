CXX	    = g++
CFLAGS	= -std=c++0x -Wall -Winline -fPIC -shared
LD	    = $(CXX)
LDFLAGS	=
RM	    = rm
INCPATH = -I$(EIGEN_DIR) -I$(EFEM_DIR)
SRCS	= mesh.cpp
LIBNAME = mesh
OBJS    = ${SRCS:.cpp=.o}

all: $(OBJS) lib

lib : $(OBJS)
	$(CXX) $(INCPATH) $(CFLAGS) $(OBJS) -o lib$(LIBNAME).so

$(OBJS) :
	$(CXX) $(INCPATH) $(CFLAGS) -c $(SRCS)

clean:
	rm -f *.o
	rm -f *.so
