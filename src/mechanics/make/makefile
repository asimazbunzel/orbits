
# MESA_DIR needed for use its numerical library
ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif

include $(MESA_DIR)/utils/makefile_header

LOAD_LIBS = -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)

ORBIT_WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FCfree) \
							-I$(MESA_INCLUDE_DIR) -c

ifndef ORBIT
	ORBIT = orbit
endif

SRCS = \
		const_def.f90 \
		mechanics_def.f90 \
		ctrls_io.f90 \
		mechanics_lib.f90

OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
OBJS_FULL = $(OBJS) run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src

$(ORBIT) : $(OBJS_FULL)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(ORBIT) $(OBJS_FULL) $(LOAD_LIBS)

run.o : $(WORK_SRC_DIR)/run.f90
	$(ORBIT_WORK_COMPILE) $<

const_def.o : $(WORK_SRC_DIR)/const_def.f90
	$(ORBIT_WORK_COMPILE) $<

mechanics_def.o : $(WORK_SRC_DIR)/mechanics_def.f90
	$(ORBIT_WORK_COMPILE) $<

mechanics_lib.o : $(WORK_SRC_DIR)/mechanics_lib.f90
	$(ORBIT_WORK_COMPILE) $<

ctrls_io.o : $(WORK_SRC_DIR)/ctrls_io.f90
	$(ORBIT_WORK_COMPILE) $<

clean :
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(ORBIT)
