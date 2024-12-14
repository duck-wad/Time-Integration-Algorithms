#pragma once

#include <vector>

void ConstantAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);

void AverageAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);

void LinearAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);
