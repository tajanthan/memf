# Makefile

CPP = g++
LDFLAGS= -lm -lstdc++
#CPPFLAGS = -g -Wall
CPPFLAGS= -DIL_STD -msse4 -O3 -s -ffast-math -funroll-loops -ansi -Wall -DNDEBUG -std=c++11

all: memf
memf: *.cpp *.h *.inc util/*.h util/*.cpp
	$(CPP) $(CPPFLAGS) *.cpp util/*.cpp -o memf $(LDFLAGS)

clean:
	rm -f memf
