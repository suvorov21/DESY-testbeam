ROOTCFLAGS 	= $(shell root-config --cflags)
ROOTLIBS   	= $(shell root-config --libs)

all: dEdx_test SR_test

dEdx_test: dEdx_test.cxx
	$(CXX) $^ -o $@.exe $(ROOTCFLAGS) $(ROOTLIBS)
SR_test: SR_test.cxx
	$(CXX) $^ -o $@.exe $(ROOTCFLAGS) $(ROOTLIBS)