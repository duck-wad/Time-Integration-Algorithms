#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>

void ReadFile(std::string fileName, std::vector<double>& f, std::vector<double>& t, double& delta_t, int& numsteps);

template<typename T>
void writeVectorToCSV(const std::vector<T>& vector, const std::string& filename) {

	std::ofstream file(filename);

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
	for (size_t i = 0; i < cols; ++i) {
		for (size_t j = 0; j < rows; ++j) {
			file << matrix[i][j];
			if (j < rows - 1) { // Add a comma unless it's the last column
				file << ",";
			}
		}
		file << "\n"; // Newline after each row
	}

	file.close();

	std::cout << "Matrix written to " << filename << " successfully." << std::endl;
}

//tranpose an nxm matrix
template<typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& matrix) {

	size_t rows = matrix.size();
	assert(rows != 0 && "Matrix must be populated to be transposed");
	size_t cols = matrix[0].size();

	std::vector<std::vector<T>> output(cols, std::vector<T>(rows));

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			output[j][i] = matrix[i][j];
		}
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}

	size_t mat1row = mat1.size();
	size_t mat1col = mat1[0].size();
	size_t mat2row = mat2.size();
	size_t mat2col = mat2[0].size();

	if (mat1col != mat2row || mat1row != mat2col) {
		throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
	}

	std::vector<std::vector<T>> output(mat2col, std::vector<T>(mat1row, T()));

	for (size_t i = 0; i < mat2col; i++) {
		for (size_t j = 0; j < mat1row; j++) {
			for (size_t k = 0; k < mat1col; k++) {
				output[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
	return output;
}

template<typename T>
std::vector<T> operator* (const std::vector<std::vector<T>>& mat, const std::vector<T>& vec) {
	if (mat.empty() || vec.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat.size() != vec.size() || mat[0].size() != vec.size()) {
		throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
	}
	std::vector<T> output(vec.size(), T());

	for (size_t i = 0; i < vec.size(); i++) {
		for (size_t j = 0; j < vec.size(); j++) {
			output[i] += mat[i][j] * vec[j];
		}
	}
	return output;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T c) {
	if (vec.empty()) {
		throw std::invalid_argument("Vector cannot be empty");
	}
	std::vector<T> output(vec.size());
	for (size_t i = 0; i < vec.size(); i++) {
		output[i] += vec[i] * c;
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	std::vector<std::vector<T>> output(mat1.size(), std::vector<T>(mat1[0].size(), T()));

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			output[i][j] = mat1[i][j] + mat2[i][j];
		}
	}
	return output;
}

template<typename T>
std::vector<std::vector<T>>& operator+=(std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {

	if (mat1.empty() || mat2.empty()) {
		throw std::invalid_argument("Matrices must not be empty");
	}
	if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
		throw std::invalid_argument("Matrices must be the same size for addition");
	}

	for (size_t i = 0; i < mat1.size(); i++) {
		for (size_t j = 0; j < mat1[0].size(); j++) {
			mat1[i][j] += mat2[i][j];
		}
	}
	return mat1;
}

template<typename T>
std::vector<T>& operator+=(std::vector<T>& v1, const std::vector<T>& v2) {

	if (v1.empty() || v2.empty()) {
		throw std::invalid_argument("Vectors must not be empty");
	}
	if (v1.size() != v2.size()) {
		throw std::invalid_argument("Vectors must be the same size for addition");
	}

	for (size_t i = 0; i < v1.size(); i++) {
		v1[i] += v2[i];
	}
	return v1;
}

template<typename T>
std::vector<std::vector<T>>& operator*= (std::vector<std::vector<T>>& matrix, T scalar) {
	for (auto& row : matrix) {
		for (auto& value : row) {
			value *= scalar;
		}
	}
	return matrix;
}