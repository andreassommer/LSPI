CXX := g++
CXXFLAGS := -std=c++11 -Wall -O2 -I.

TARGET_1D = test_1D_local_interpolation
TARGET_2D = test_2D_local_interpolation

.PHONY: all clean dim1 dim2

all: $(TARGET_1D) $(TARGET_2D)

dim1: $(TARGET_1D)
dim2: $(TARGET_2D)

$(TARGET_1D): test_1D_local_interpolation.cpp LSPI_tools.hpp LSPI_1D.hpp tictoc.hpp
	@echo "Building $@"
	$(CXX) $(CXXFLAGS) -o $@ $<

$(TARGET_2D): test_2D_local_interpolation.cpp LSPI_tools.hpp LSPI_2D.hpp tictoc.hpp
	@echo "Building $@"
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f *.o $(TARGET_1D) $(TARGET_2D)

.PHONY: run