#pragma once

#include <iostream>
#include <vector>

void ReadFile(std::string fileName, std::vector<double>& f, std::vector<double>& t, double& delta_t, int& numsteps);

template<typename T>
void writeVectorToCSV(const std::vector<T>& vector, const std::string& filename) {

	std::string folderPath = "data";
	// Open the file stream

	std::ofstream file(folderPath + "\\" + filename);

	if (!file.is_open()) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	for (size_t i = 0; i < vector.size(); i++) {
		file << vector[i] << "\n";
	}

	file.close();

	std::cout << "Vector written to " << filename << " succesfully." << std::endl;
}