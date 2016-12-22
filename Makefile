# Default prefix. Can change this to system directories if needed
PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

# Intel compiler - uncomment if you have icpc and mkl.
#CXX=icpc
#CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
#D_LDLIBS=-L$(PREFIX)/lib -larmadillo -lboost_program_options -mkl
# gcc
#CXXFLAGS=-Wall -O3 -std=c++11
# gcc test
CXXFLAGS=-Wall -g -O0 -std=c++11
D_LDLIBS=-L$(PREFIX)/lib -larmadillo -lboost_program_options -llapack -lblas

CPPFLAGS=-I$(PREFIX)/include

PROGRAMS=tajima

OBJECTS=cmdLine.o tajima.o

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

tajima: $(OBJECTS)
	$(LINK.cpp) $^ $(D_LDLIBS) -o $@

.PHONY: all clean install

