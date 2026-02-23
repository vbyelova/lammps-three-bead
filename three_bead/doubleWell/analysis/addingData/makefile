CXX=g++
CXXFLAGS=-I. -fPIC
LDFLAGS=-lstdc++

.PHONY: all
all: addingData addingData.o

addingData.o: addingData.h addingData.cpp percPerBC.h
	 $(CXX) $(CXXFLAGS) -c addingData.cpp -o addingData.o $(LDFLAGS)

.PHONY: clean
clean:
	rm -f *.o addingData
