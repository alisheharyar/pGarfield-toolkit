OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

BOOST_ROOT=/lustre/sw/xc40ac/boost/1.58.0/cnl5.2_gnu5.2.0

# Compiler flags
CFLAGS += -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) \
	-I/opt/intel/ics/2013/028/impi/4.1.0.024/intel64/include/ -DMPICH_IGNORE_CXX_SEEK

CFLAGS += -I/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/include
CFLAGS += -I$(BOOST_ROOT)/include
CFLAGS += -DPARALLEL -DMPICH_IGNORE_CXX_SEEK
CFLAGS += -I/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/include

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
LDFLAGS += -L$(BOOST_ROOT)/lib -lboost_program_options
#LDFLAGS += -lpthread
LDFLAGS += -L/opt/cray/mpt/7.3.1/gni/mpich-gnu/4.9/lib -lmpich 


#CXX = mpicc

.PHONY:	gem


gem: gem.C
	$(CXX) $(CFLAGS) gem.C
	$(CXX) -o gem gem.o $(LDFLAGS)
	rm gem.o

clean:
	rm *.o gem
