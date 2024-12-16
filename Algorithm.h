#pragma once

#include <vector>

//SDOF algorithms
void ConstantAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);
void AverageAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);
void LinearAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k);

//MDOF algorithms
void ConstantAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& F, std::vector<double>& D_i, std::vector<double>& V_i, int numsteps, double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);
void AverageAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, std::vector<double>& D_i, std::vector<double>& V_i, int numsteps, double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);
void LinearAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, std::vector<double>& D_i, std::vector<double>& V_i, int numsteps, double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);