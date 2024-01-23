ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX  = g++
CXX += -I./	
CXX += -I./obj/

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += -Wreorder

OUTLIB = ./obj/

#----------------------------------------------------#

all: run_dbscan

run_dbscan: run_dbscan.cpp obj/dbscan.o
	$(CXX) -o run_dbscan run_dbscan.cpp $(OUTLIB)*.o $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS)

obj/dbscan.o: src/dbscan.cpp src/dbscan.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)/dbscan.o $<

clean:
	rm -f run_dbscan
	rm -f *~
	rm -f src/*~
	rm -f $(OUTLIB)*.o
