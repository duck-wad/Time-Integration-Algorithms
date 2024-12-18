#include <iostream>
#include <vector>

#include "Utils.h"
#include "Algorithm.h"

void TestSDOF() {
	//input data from input file
	std::string fileName = "INPUT_SDOF.txt";
	std::vector<double> force;
	std::vector<double> time;
	double delta_t;
	int numsteps;

	ReadFile(fileName, force, time, delta_t, numsteps);

	double initial_d = 0.0;
	double initial_v = 0.0;

	double mass = 1.0;
	double damping = 1.0;
	double stiffness = 1.0;

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

	writeMatrixToCSV(transpose(results), "OUTPUT_SDOF.csv");
}

//3 node system. force is applied to one node (last in vector)
void TestMDOF() {
	//input data from input file, for one forced node only
	std::string fileName = "INPUT_MDOF.txt";
	std::vector<double> force_node;
	std::vector<double> time;
	double delta_t;
	int numsteps;

	ReadFile(fileName, force_node, time, delta_t, numsteps);

	//initialize the force vector. force applied to the third node
	std::vector<std::vector<double>> force(numsteps, std::vector<double>(3, 0.0));
	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		force[i][2] = force_node[i];
	}
	//initial conditions
	std::vector<double> initial_d = { 0.0, 0.0, 0.0 };
	std::vector<double> initial_v = { 0.0, 0.0, 0.0 };
	//mass, damping, stiffness matrices. row major order
	std::vector<std::vector<double>> mass = {
		{100.0, 0.0, 0.0},
		{0.0, 100.0, 0.0},
		{0.0, 0.0, 50.0}
	};
	std::vector<std::vector<double>> stiffness = {
	{2.0e7, -1.0e7, 0.0},
	{-1.0e7, 2.5e7, -0.5e7},
	{0.0, -0.5e7, 0.5e7}
	};
	std::vector<std::vector<double>> damping = {
	{5000.0, 0.0, 0.0},
	{0.0, 2500.0, -1000.0},
	{0.0, -1000.0, 1000.0}, 
	};

	std::vector<std::vector<double>> displacement(numsteps, std::vector<double>(3, 0.0));
	std::vector<std::vector<double>> velocity(numsteps, std::vector<double>(3, 0.0));
	std::vector<std::vector<double>> acceleration(numsteps, std::vector<double>(3, 0.0));
	int nodes = 3;

	std::vector<std::vector<double>> results;
	results.push_back(time);
	results.push_back(force_node);

	ConstantAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness, nodes);

	results.push_back(transpose(displacement)[0]);
	results.push_back(transpose(displacement)[1]);
	results.push_back(transpose(displacement)[2]);

	AverageAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness, nodes);

	results.push_back(transpose(displacement)[0]);
	results.push_back(transpose(displacement)[1]);
	results.push_back(transpose(displacement)[2]);

	LinearAccelerationMethod(displacement, velocity, acceleration, force, initial_d, initial_v, numsteps, delta_t, mass, damping, stiffness, nodes);

	results.push_back(transpose(displacement)[0]);
	results.push_back(transpose(displacement)[1]);
	results.push_back(transpose(displacement)[2]);

	writeMatrixToCSV(transpose(results), "OUTPUT_MDOF.csv");
}

int main() {

	TestSDOF();
	TestMDOF();

	return 0;
}