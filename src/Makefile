# Include configuration file 
#CN=LAPTOP
#C=ifort
#LIB=mkl
#MPI=impi+omp
-include ../config/${C}+${LIB}+${MPI}_$(CN).cnf

# Commands
GIT_VERSION := $(shell (tail -1 ../../.git/logs/HEAD | awk '{print $$2}'))
DATE := $(shell date +%d-%b-%G)

MKDIR_P = mkdir -p
OUT_DIR = ../bin

QMECHA_EXE = $(OUT_DIR)/QMeCha

# List of directories
-include directories

# List of executables to be built within the package
-include fileslist

# "make" builds all
all: builddirs qmecha 

builddirs: 
	$(MKDIR_P) $(OUT_DIR)

qmecha: builddirs version QMeCha.o $(QMECHA_FLS) 
	$(FC) $(FCFLAGS) QMeCha.o $(QMECHA_FLS) $(LDFLAGS) -o $(QMECHA_EXE) 

version: 
	sed "s/LLLDATELLL/$(DATE)/g" $(UTILS_DIR)/qmecha_version.tmp > $(UTILS_DIR)/qmecha_version.f90

-include dependencies

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -rf *.o *.mod */*.mod */*.o */*/*.mod */*/*.o

veryclean: clean
	rm -rf $(OUT_DIR)
