include $(FEPIC_DIR)/conf/variables

.PHONY: src

all: src ok_msg

src:
	$(MAKE) all -C Fepic/src
	
clean:
	$(MAKE) clean -C Fepic/src
	
ok_msg: src
	@echo
	@echo successfully compiled
	@echo
	
