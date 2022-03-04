LIB_DIR = lib
HOSTCOMPUTER := $(shell hostname | cut -d . -f 1)

default: milos

clean-lib:
	rm -f $(LIB_DIR)/*.o 
	rm -f $(LIB_DIR)/*.a *.so 
	rm -f $(LIB_DIR)/*.so 
	rm -f $(LIB_DIR)/milos 
	rm -f milos.$(HOSTCOMPUTER).x 

clean-cython:
	rm -f *.so
	rm -Rf build
	rm -f pymilos.c

clean: clean-lib clean-cython

.PHONY: milos
milos: clean-lib
	make -C $(LIB_DIR) default

.PHONY: build
build: clean-lib
	make -C $(LIB_DIR) default
	python3 setup.py build_ext --inplace 

.PHONY: cython-build
cython-build: clean
	make -C $(LIB_DIR) default
	python3 setup.py build_ext --inplace  --use-cython 

