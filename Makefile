.PHONY: gvar

OS_NAME = $(shell uname)
PLATFORM = $(shell uname -p)
ARCH = $(OS_NAME)_$(PLATFORM)

########################################################################
all: gvar

gvar:
	cd makefiles/$(ARCH); make -f Makefile$(GCC_VER) ;
	cp makefiles/$(ARCH)/gvar . 

clean:
	cd makefiles/$(ARCH); make -f Makefile$(GCC_VER) clean
	rm -f gvar

