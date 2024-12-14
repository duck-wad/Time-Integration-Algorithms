#include <fstream>

#include "Utils.h"

void ReadFile(std::string fileName, std::vector<double>& f, std::vector<double>& t, double& delta_t, int& numsteps) {

	std::string junk;

	std::ifstream infile(fileName);
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
	}

	infile >> junk >> numsteps >> junk >> delta_t;

	t.resize(numsteps, 0.0);
	f.resize(numsteps, 0.0);

	for (size_t i = 0; i < numsteps; i++) {
		infile >> junk >> t[i] >> junk >> f[i];
	}
}