TARGET = fluid_sim
CXX = g++-14
CXXFLAGS = -O3 -fopenmp -Wall -std=c++17
LDFLAGS =  -lomp -L/opt/homebrew/opt/libomp/lib

all:
	$(CXX) src/main.cpp src/simulation.cpp src/geometry.cpp src/optimal_transport.cpp src/utils.cpp lbfgs.c \
	    -o $(TARGET) $(CXXFLAGS) $(LDFLAGS) -I. -Isrc

clean:
	rm -f $(TARGET) fluid_sim_*.png

run: all
	./$(TARGET)
