include $(FEPIC_DIR)/conf/variables

RM	  = rm

CPPSOURCES = $(wildcard *.cpp)
OBJS	= ${SRCS:.cpp=.o}

all: $(CPPSOURCES:.cpp=.o)

%.o: %.cpp
	$(FEP_CXX) -c $< $(FEP_CXXFLAGS) -o $@
	
clean:
	$(FEP_RM) *.o *~
	
remade:
	$(MAKE) clean
	$(MAKE)

-include $(CPPSOURCES:.cpp=.d)

%.d: %.cpp
	$(FEP_CXX) $< -MM -MT '$*.o $*.d ' -MD $(FEP_CXXFLAGS)
