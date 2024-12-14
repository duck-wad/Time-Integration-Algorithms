#include <iostream>
#include <vector>

#include "Utils.h"
#include "Algorithm.h"

int main() {

	//input data from input file
	std::string fileName = "INPUT.txt";
	std::vector<double> force;
	std::vector<double> time;
	double delta_t;
	int numsteps;

	double initial_d = 0.0;
	double initial_v = 0.0;

	double mass = 1.0;
	double damping = 1.0;
	double stiffness = 1.0;

	ReadFile(fileName, force, time, delta_t, numsteps);

	std::vector<std::vector<double>> results;
	results.push_back(time);
	results.push_back(force);

	std::vector<double> displacement;
	std::vector<double> velocity;
	std::vector<double> acceleration;

	ConstantAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness);

	results.push_back(displacement);
	results.push_back(velocity);
	results.push_back(acceleration);

	AverageAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness);

	results.push_back(displacement);
	results.push_back(velocity);
	results.push_back(acceleration);

	LinearAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness);

	results.push_back(displacement);
	results.push_back(velocity);
	results.push_back(acceleration);

	writeMatrixToCSV(results, "OUTPUT.csv");

	return 0;
}