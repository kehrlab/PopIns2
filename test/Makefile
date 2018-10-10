
CXX=g++ -std=c++14
CC=$(CXX)

LDLIBS = -lbifrost -pthread -lz -lrt -rdynamic

CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result
CXXFLAGS+= -march=native
CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

all: test_popins2

test_popins2:test_popins2.o ../build/ColoredDeBruijnGraph.o ../build/UnitigExtension.o
test_popins2.o: test_popins2.cpp $(HEADERS)

clean:
	rm -f *.o test_popins2

purge:
	rm -f *.o test_popins2 *.gfa *_colors
