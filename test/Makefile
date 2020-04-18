
CXX=g++ -std=c++14
CC=$(CXX)

LDLIBS= -lbifrost -pthread -lz -lrt -rdynamic
LDLIBS+= -DMAX_KMER_SIZE=64

CXXFLAGS+= -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result
CXXFLAGS+= -march=native
CXXFLAGS+= -DMAX_KMER_SIZE=64
CXXFLAGS+= -g -O0 -DDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

all: test_popins2

test_popins2:test_popins2.o ../build/ColoredDeBruijnGraph.o ../build/UnitigExtension.o ../build/Traceback.o ../build/LECC_Finder.o ../build/Setcover.o
test_popins2.o: test_popins2.cpp $(HEADERS)

clean:
	rm -f *.o test_popins2

purge:
	rm -f *.o test_popins2 *.gfa *.bfg_colors *.csv *.log
