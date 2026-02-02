CXX=g++
CXXFLAGS=-I.
LDFLAGS=-lstdc++

.PHONY: all
all: testPercPerBC percPerBC.o

percPerBC.o: percPerBC.h percPerBC.cpp
	$(CXX) $(CXXFLAGS) -c percPerBC.cpp -o percPerBC.o

testPercPerBC: testPercPerBC.cpp percPerBC.o
	$(CXX) $(CXXFLAGS) -o testPercPerBC testPercPerBC.cpp percPerBC.o $(LDFLAGS)

.PHONY: clean
clean:
	rm -f *.o testPercPerBC
