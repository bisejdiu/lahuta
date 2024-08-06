CXX = g++ -std=c++17
CC = gcc

GEMMI = external/gemmi
OB = external/ob

CXXFLAGS = -DDYNAMIC_CRC_TABLE -Iinclude -I$(GEMMI)/include -I$(OB) -I$(GEMMI)/third_party -I$(GEMMI)/third_party/zlib -I$(CONDA_PREFIX)/include/rdkit -O3

LDFLAGS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
LDLIBS =  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib -lRDKitSubstructMatch -lRDKitSmilesParse

SRCS = $(wildcard src/*cpp src/bonds/*.cpp src/ob/*.cpp) $(wildcard $(GEMMI)/src/*.cpp) $(wildcard $(OB)/*.cpp) 
ZLIB_SRCS = $(wildcard $(GEMMI)/third_party/zlib/*.c)
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
