#pragma once

#include <iostream>
#include <vector>
#include <fstream>

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

template<typename T>
void writeMatrixToCSV(const std::vector<std::vector<T>>& matrix, const std::string& filename) {

	// Open the file stream

	std::ofstream file(filename);

	if (!file.is_open()) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	size_t rows = matrix[0].size();
	size_t cols = matrix.size();

	// Write the matrix to the file in row-major order
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			file << matrix[j][i];
			if (j < cols - 1) { // Add a comma unless it's the last column
				file << ",";
			}
		}
		file << "\n"; // Newline after each row
	}

	file.close();

	std::cout << "Matrix written to " << filename << " successfully." << std::endl;
}