#pragma once

#include <vector>

//SDOF algorithms
void ConstantAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k);
void AverageAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k);
void LinearAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k);

//MDOF algorithms
void ConstantAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);
void AverageAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);
void LinearAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);