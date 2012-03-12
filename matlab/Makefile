

include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk


ifdef OCTAVE
  MEX ?= mkoctfile --mex
  MATLAB ?= octave
  MEXEXT ?= mex
else
  MEX ?= mex
  MATLAB ?= matlab -nodisplay

  UNAME := $(shell uname)

  # determine the MEX extension

  ifeq ($(UNAME), Linux)
    ifeq ($(shell getconf LONG_BIT),64)
	MEXEXT ?= mexa64
    else
	MEXEXT ?= mexglx
    endif
  else
    ifeq ($(UNAME), Darwin)
	ifeq ($(shell getconf LONG_BIT),64)
	  MEXEXT ?= mexmaci
	else
	  MEXEXT ?= mexmaci64
	endif
    else
      # assume windows
	MEXEXT ?= mexw32
    endif
  endif
endif

# add Fortran run-time library
LIBS += $(FRTLIB) 

ifdef RPATH
  LIBS += -Wl,-R$(RPATH)
endif

%.$(MEXEXT): %.c
	$(MEX)  -g $<  -L.. -loak $(LIBS)

.SUFFIXES: $(SUFFIXES) $(MEXEXT)

all: mex

mex: oak_init.$(MEXEXT) oak_done.$(MEXEXT) oak_analysis.$(MEXEXT)

check:
	(cd ../matlab; cat ../Mex/oak_check_mex.m | $(MATLAB))

clean: 
	rm -f *.$(MEXEXT)