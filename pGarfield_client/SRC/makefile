OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include/Garfield
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

BOOST_ROOT=/lustre/sw/xc40/cle7/boost/1.76.0

INTEL_PATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64

# Compiler flags
CFLAGS += -Wall -Wextra -Wno-long-long \
  `root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) \
	-I$(INTEL_PATH)/include64 -DMPICH_IGNORE_CXX_SEEK

#CFLAGS += -I/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/include
CFLAGS += -I$(INTEL_PATH)/include
CFLAGS += -I$(BOOST_ROOT)/include
CFLAGS += -DPARALLEL -DMPICH_IGNORE_CXX_SEEK

#CFLAGS += -I/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/include

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
LDFLAGS += -L$(BOOST_ROOT)/lib -lboost_program_options
LDFLAGS += -Wl,--no-relax

#LDFLAGS += -lpthread
#LDFLAGS += -L/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/lib -lmpich 
#LDFLAGS += -L$(INTEL_PATH)/lib -lmpicxx -lmpifort


CXX = CC

.PHONY:	gem


gem: gem.C
	$(CXX) $(CFLAGS) gem.C
	$(CXX) -o gem gem.o $(LDFLAGS)
	rm gem.o
 
gem2B: gem2B.C
	$(CXX) $(CFLAGS) gem2B.C
	$(CXX) -o gem2B gem2B.o $(LDFLAGS)
	rm gem2B.o

clean:
	rm *.o gem
