#	Makefile
#
#	Assumes Root is installed in $HOME/root.
#
#	TARGETS:
#		all	build the usual program (default)
#		sptest	build a test program for SimplePlot
#
#	NOTE: This has only been tested on Mac OS X 10.9.5 (Mavericks).

# list source files here 
SRC = *.c

ROOT = $(HOME)/root

CFLAGS = $(shell $(ROOT)/bin/root-config --cflags)
LIBS =   $(shell $(ROOT)/bin/root-config --libs)
RPATH = -Wl,-rpath,$(ROOT)/lib

all:	program sptest

program: $(SRC) SimplePlot.cc SimplePlot.hh
	g++ $(CFLAGS) -o program $(SRC) SimplePlot.cc $(RPATH) $(LIBS)

sptest:	SimplePlot.hh SimplePlot.cc
	g++ $(CFLAGS) -DTEST -o sptest SimplePlot.cc $(RPATH) $(LIBS)


clean:
	rm -f *.o program sptest
