CXX = g++
CC = gcc

ROOT = /home/bisejdiu/raw_dev/lahuta_cpp/gemmi
INCLUDE_DIR = $(ROOT)/include

CXXFLAGS = -DDYNAMIC_CRC_TABLE -I$(INCLUDE_DIR) -Ibond_table -I$(ROOT)/third_party -I$(ROOT)/third_party/zlib -I$(CONDA_PREFIX)/include/rdkit -O3

LDFLAGS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib -lRDKitSubstructMatch -lRDKitSmilesParse

SRCS = $(wildcard ./*cpp) $(wildcard bond_table/*.cpp) $(wildcard $(ROOT)/src/*.cpp) 
ZLIB_SRCS = $(wildcard $(ROOT)/third_party/zlib/*.c)
OBJS = $(SRCS:.cpp=.o) $(ZLIB_SRCS:.c=.o)

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
