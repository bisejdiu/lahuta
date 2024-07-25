# Define the compiler to use
# CXX = g++ -mavx2 -mfma 
CXX = g++ -g

ROOT = /home/bisejdiu/raw_dev/lahuta_cpp/gemmi

CXXFLAGS = -DDYNAMIC_CRC_TABLE -I$(ROOT)/lahuta/lahuta_cpp/bond_table -I$(ROOT)/lahuta/bond_order -I$(ROOT)/include -I$(ROOT)/third_party -I$(ROOT)/third_party/zlib -I$(CONDA_PREFIX)/include/rdkit -I/home/bisejdiu/tmp/dpnblist/only_cpu/ -O0

# Define any directories containing libraries
LDFLAGS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib

# Define any libraries to link into executable
# LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitDetermineBonds -lRDKitRDGeometryLib -lRDKitRDBoost -lRDKitSubgraphs -lRDKitSubstructMatch -lRDKitSmilesParse
LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib -lRDKitSubstructMatch -lRDKitSmilesParse

SRCS = main.cpp bonds.cpp conv.cpp nsgrid.cpp elements.cpp bond_utils.cpp kekulize.cpp bitvec.cpp bond_order.cpp bond_table/bonds.cpp bond_table/tokens_table.cpp $(wildcard $(ROOT)/src/*.cpp) $(ROOT)/prog/options.cpp $(wildcard /home/bisejdiu/tmp/dpnblist/only_cpu/*.cpp)
ZLIB_SRCS = $(wildcard $(ROOT)/third_party/zlib/*.c)
OBJS = $(SRCS:.cpp=.o) $(ZLIB_SRCS:.c=.o)

# Define the executable file
TARGET = rk

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
