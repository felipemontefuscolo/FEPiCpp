include $(FEPIC_DIR)/conf/variables

SUBDIRS = custom_eigen mesh mypetsc mesh_tools quadrature shapefunctions util dofhandler

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

all: subdirs build ok_msg

$(SUBDIRS):
	$(MAKE) all -C $@
	
clean:
	$(FEP_RM) `find . -iname *.o`
	
# pay attention to the wildcard depth.
build: subdirs
	ar -cr libfepic.a $(wildcard */*.o) $(wildcard */*/*.o) $(wildcard */*/*/*.o)
	mv -f libfepic.a $(FEP_LIBS_DIR)
	ranlib $(FEP_LIBS_DIR)/libfepic.a

ok_msg: build
	@echo
	@echo compiled successfully
	@echo
	
	
	
