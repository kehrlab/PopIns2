TARGET = popins2
BUILD_DIR = ./build
SRC_DIR = ./src

SRCS := $(shell find $(SRC_DIR) -type f -name *.cpp)
OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.cpp=.o))

# Compiler
CXX = g++ -std=c++14
CC = $(CXX)

-include popins2.config
TOOLS=-DSAMTOOLS=\"$(SAMTOOLS)\" -DBWA=\"$(BWA)\" -DSICKLE=\"$(SICKLE)\" -DVELVETH=\"$(VELVETH)\" -DVELVETG=\"$(VELVETG)\" -DMINIA=\"$(PWD)/$(MINIA)\"

# Date and version number from git
DATE := on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
VERSION := 0.9.2-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)
CXXFLAGS += -DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Compiler flags
CXXFLAGS += -DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK
CXXFLAGS += -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result
CXXFLAGS += -march=native -DMAX_KMER_SIZE=128

# Linker flags
LDLIBS = -lbifrost -pthread -lz -lrt -rdynamic -DMAX_KMER_SIZE=128

# DEBUG build
#CXXFLAGS += -g -pg -O0 -DDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS += -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) $(TOOLS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

purge:
	rm -f $(OBJS) $(TARGET) *.gfa *.bfg_colors *.fasta *.csv popins2*log

metaclean:
	rm -f $(TARGET) *.gfa *.bfg_colors *.fasta *.csv popins2*log
