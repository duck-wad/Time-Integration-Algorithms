#include <iostream>
#include <vector>

#include "Utils.h"
#include "Algorithm.h"

int main() {

	std::string fileName = "INPUT10s.txt";
	std::vector<double> force;
	std::vector<double> time;
	double delta_t;
	int numsteps;
	std::vector<double> displacement;
	std::vector<double> velocity;
	std::vector<double> acceleration;

	ReadFile(fileName, force, time, delta_t, numsteps);

	std::cout << "Hello world" << std::endl;

	return 0;
}