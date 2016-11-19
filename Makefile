all: main.cpp
	$(CXX) -std=c++11 -lGLU -lglut -lGL -fopenmp -O3 -march=native -Wall -o tutorial $<

clean:
	rm tutorial
