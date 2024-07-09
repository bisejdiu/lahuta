# Define the compiler to use
# CXX = g++ -mavx2 -mfma 
CXX = g++ -g

# Define any compile-time flags
# CXXFLAGS = -I$(CONDA_PREFIX)/include/rdkit -I$(CONDA_PREFIX)/include/gemmi -I/home/bisejdiu/raw_dev/lahuta_cpp/gemmi/lahuta/SimSIMD/include -Wall -O0
#
#

ROOT = /home/bisejdiu/raw_dev/lahuta_cpp/gemmi

# CXXFLAGS = -DDYNAMIC_CRC_TABLE -I$(ROOT)/include -I$(ROOT)/third_party -I$(ROOT)/third_party/zlib -I$(CONDA_PREFIX)/include/rdkit -I/home/bisejdiu/raw_dev/lahuta_cpp/gemmi/lahuta/SimSIMD/include 

CXXFLAGS = -DDYNAMIC_CRC_TABLE -I$(ROOT)/lahuta/bond_order -I$(ROOT)/include -I$(ROOT)/third_party -I$(ROOT)/third_party/zlib -I$(CONDA_PREFIX)/include/rdkit -I/home/bisejdiu/tmp/dpnblist/only_cpu/ -O3


# GEMMI_BASE = /home/bisejdiu/raw_dev/lahuta_cpp/gemmi
# CXXFLAGS = -I$(CONDA_PREFIX)/include/rdkit -I$(GEMMI_BASE)/include/gemmi -I$(GEMMI_BASE)/include/gemmi/third_party -I$(GEMMI_BASE)/third_party -I$(GEMMI_BASE)/include/gemmi/third_party/tao -I$(GEMMI_BASE)/include/gemmi/third_party/tao/pegtl -I$(GEMMI_BASE)/third_party/zlib -I/home/bisejdiu/raw_dev/lahuta_cpp/gemmi/lahuta/SimSIMD/include -Wall -O0

# Define any directories containing libraries
LDFLAGS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib

# Define any libraries to link into executable
# LDLIBS = -lrdkit -lRDKit -lRDGeneral -lRDGeometryLib -lRDBoost -lFileParsers -lGraphMol -lDataStructs -lSmilesParse -lSubstructMatch -lDescriptors -lChemTransforms -lChemReactions -lChemTransforms -lChemReactions -lChemInformatics
# LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitDetermineBonds -lgemmi_cpp
LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitDetermineBonds -lRDKitRDGeometryLib -lRDKitRDBoost -lRDKitSubgraphs -lRDKitSubstructMatch -lRDKitSmilesParse


# LIBC = /home/bisejdiu/raw_dev/lahuta_cpp/gemmi/lahuta/SimSIMD/c/lib.c

SRCS = main.cpp bonds.cpp conv.cpp nsgrid.cpp elements.cpp bond_utils.cpp kekulize.cpp bitvec.cpp bond_order.cpp $(wildcard $(ROOT)/src/*.cpp) $(ROOT)/prog/options.cpp $(wildcard /home/bisejdiu/tmp/dpnblist/only_cpu/*.cpp)
ZLIB_SRCS = $(wildcard $(ROOT)/third_party/zlib/*.c)
OBJS = $(SRCS:.cpp=.o) $(ZLIB_SRCS:.c=.o)


# Define the executable file
TARGET = rk

# Define the source files
# SOURCES = init.cpp bonds.cpp $(GEMMI_BASE)/src/* $(GEMMI_BASE)/third_party/zlib/*.c

# Define the object files
# OBJECTS = $(SOURCES:.cpp=.o)

#
# Default target
# all: $(TARGET)
#
# $(TARGET): $(SOURCES)
# 	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS) $(LDLIBS)
#
# # To remove generated files
# clean:
# 	rm -f $(TARGET)


# # Define rules
all: $(TARGET)

$(TARGET): main.o $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) main.o $(TARGET)

.PHONY: all clean
