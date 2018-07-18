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
	-I$(INCDIR) -I$(HEEDDIR)

CFLAGS += -I$(BOOST_ROOT)/include
CFLAGS += -DPARALLEL


LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
LDFLAGS += -L$(BOOST_ROOT)/lib -lboost_program_options


CXX = mpicxx

.PHONY:	gem


gem: gem.C
	$(CXX) $(CFLAGS) gem.C
	$(CXX) -o gem gem.o $(LDFLAGS)
	rm gem.o

clean:
	rm *.o gem
