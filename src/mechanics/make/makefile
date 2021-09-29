
FC = gfortran

# Set optimization level
ARCH = $(shell uname -m)

ifeq (${ARCH},x86_64)
	OPT_LEVEL = 2
else
	OPT_LEVEL = 0
endif

FCbasic = -Wno-uninitialized -fno-range-check -fmax-errors=7 $(FCbasic2)
FCbasic2 = -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan
FCstatic =
FCopenmp = -fopenmp
FCstandard = -std=f2008 -Wno-error=tabs
FCfree = -ffree-form -ffree-line-length-none -x f95-cpp-input
FCwarn = -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination
FCimpno = -fimplicit-none
FCopt = -O${OPT_LEVEL} -ftree-vectorize
FCchecks = -fbounds-check -Wuninitialized -Warray-bounds
FCdebug = -ggdb


COMPILE_BASIC_FLAGS = $(FCbasic) $(FCopenmp) $(FCstatic) $(FCstandard)
COMPILE_BASIC = $(FC) $(COMPILE_BASIC_FLAGS)
COMPILE_TO_DEPLOY = $(COMPILE_BASIC) $(FCwarn) $(FCimpno) $(FCopt) -c
COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree)
COMPILE_CMD = $(COMPILE)
LOADER = $(FC)
WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FCfree) -c

ifndef ORBIT
	ORBIT = orbit
endif

SRCS = \
		const_def.f90 \
		mechanics_def.f90 \
		mechanics_lib.f90

OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
OBJS_FULL = $(OBJS) run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src

$(ORBIT) : $(OBJS_FULL)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(ORBIT) $(OBJS_FULL)


run.o : $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<

const_def.o : $(WORK_SRC_DIR)/const_def.f90
	$(WORK_COMPILE) $<

mechanics_def.o : $(WORK_SRC_DIR)/mechanics_def.f90
	$(WORK_COMPILE) $<

mechanics_lib.o : $(WORK_SRC_DIR)/mechanics_lib.f90
	$(WORK_COMPILE) $<

clean :
	-@rm -f *.o *.mod $(WORK_DIR)/$(ORBIT)
