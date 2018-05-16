CXX=g++ -std=c++14
CC=$(CXX)

LDLIBS = -lbifrost -pthread -lz -lrt -rdynamic

CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result
CXXFLAGS+= -march=native
CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

all: test_color_print

test_color_print:test_color_print.o ../build/ColoredCDBG_Data_extension.o ../build/ColoredCDBG_Graph_extension.o
test_color_print.o: test_color_print.cpp $(HEADERS)

test_CDBG:test_CDBG.o ../build/CompactedDBG_Data_extension.o ../build/CompactedDBG_Graph_extension.o
test_CDBG.o: test_CDBG.cpp $(HEADERS)

test_popins2:test_popins2.o ../build/CDBG_Data_extension.o ../build/ExtendedCDBG.o ../build/argument_parsing.o
test_popins2.o: test_popins2.cpp $(HEADERS)

clean:
	rm -f *.o test_popins2 test_CDBG test_color_print
