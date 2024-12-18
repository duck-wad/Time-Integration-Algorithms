#include <iostream>
#include <cmath>

#include "Algorithm.h"
#include "Utils.h"

/* IMPLEMENTATION OF CONSTANT, AVERAGE, AND LINEAR ACCELERATION METHODS FOR SINGLE DEGREE OF FREEDOM SYSTEM */

//constant acceleration for a single degree of freedom system
void ConstantAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.assign(numsteps, 0.0);
	v.assign(numsteps, 0.0);
	a.assign(numsteps, 0.0);

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
		}
		else {
			d[i] = a[i - 1] / 2.0 * delta_t * delta_t + v[i - 1] * delta_t + d[i - 1];
			v[i] = a[i - 1] * delta_t + v[i - 1];
		}
		a[i] = 1 / m * (f[i] - k * d[i] - c * v[i]);
	}
}

//average acceleration method for SDOF
void AverageAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.assign(numsteps, 0.0);
	v.assign(numsteps, 0.0);
	a.assign(numsteps, 0.0);

	double k_ = (4.0 / (delta_t * delta_t) * m) + (2.0 / delta_t * c) + k;
	double p_ = 0.0;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
			a[i] = 1.0 / m * (f[i] - k * d[i] - c * v[i]);
		}
		else {
			p_ = f[i] + m * (4.0 / (delta_t * delta_t) * d[i - 1] + 4.0 / delta_t * v[i - 1] + a[i - 1])
				+ c * (2.0 / delta_t * d[i - 1] + v[i - 1]);
			d[i] = p_ / k_;
			v[i] = 2.0 / delta_t * (d[i] - d[i - 1]) - v[i - 1];
			a[i] = 4.0 / (delta_t * delta_t) * (d[i] - d[i - 1] - delta_t * v[i - 1]) - a[i - 1];
		}
	}
}

//linear acceleration method for SDOF
void LinearAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, const double d_i, const double v_i, const int numsteps, const double delta_t, const double m, const double c, const double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.assign(numsteps, 0.0);
	v.assign(numsteps, 0.0);
	a.assign(numsteps, 0.0);

	double k_ = (6.0 / (delta_t * delta_t) * m) + (3.0 / delta_t * c) + k;
	double p_ = 0.0;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
			a[i] = 1.0/ m * (f[i] - k * d[i] - c * v[i]);
		}
		else {
			p_ = f[i] + m * (6.0 / (delta_t * delta_t) * d[i - 1] + 6.0 / delta_t * v[i - 1] + 2.0 * a[i - 1])
				+ c * (3.0 / delta_t * d[i - 1] + 2.0 * v[i - 1] + delta_t / 2.0 * a[i-1]);
			d[i] = p_ / k_;
			v[i] = 3.0 / delta_t * (d[i] - d[i - 1]) - 2.0 * v[i - 1] - delta_t / 2.0 * a[i-1];
			a[i] = 6.0 / (delta_t * delta_t) * (d[i] - d[i - 1] - delta_t * v[i - 1]) - 2.0 * a[i-1];
		}
	}
}

/* IMPLEMENTATION OF MULTIPLE DEGREE OF FREEDOM ALGORITHMS */

void ConstantAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes) {
	if (F.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	D.assign(numsteps, std::vector<double>(nodes, 0.0));
	V.assign(numsteps, std::vector<double>(nodes, 0.0));
	A.assign(numsteps, std::vector<double>(nodes, 0.0));

	std::vector<std::vector<double>> inverseM = invertMatrix(M);

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			D[i] = D_i;
			V[i] = V_i;
		}
		else {
			D[i] = (A[i - 1] / 2.0 * delta_t * delta_t) + (V[i - 1] * delta_t) + D[i - 1];
			V[i] = (A[i - 1] * delta_t) + V[i - 1];
		}

		A[i] = inverseM * (F[i] - (K * D[i]) - (C * V[i]));
	}
}

void AverageAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes) {

	if (F.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	D.assign(numsteps, std::vector<double>(nodes, 0.0));
	V.assign(numsteps, std::vector<double>(nodes, 0.0));
	A.assign(numsteps, std::vector<double>(nodes, 0.0));

	std::vector<std::vector<double>> K_ = M * (4.0 / (delta_t * delta_t)) + C * (2.0 / delta_t) + K;
	std::vector<double> P_;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			D[i] = D_i;
			V[i] = V_i;
			A[i] = invertMatrix(M) * (F[i] - K * D[i] - C * V[i]);
		}
		else {
			P_ = F[i] + M * (D[i - 1] * (4.0 / (delta_t * delta_t)) + V[i - 1] * (4.0 / delta_t) + A[i - 1])
				+ C * (D[i-1] * (2.0 / delta_t) + V[i - 1]);
			D[i] = invertMatrix(K_) * P_;
			V[i] = (D[i] - D[i - 1]) * (2.0 / delta_t) - V[i - 1];
			A[i] = (D[i] - D[i - 1] - V[i-1] * delta_t) * (4.0 / (delta_t * delta_t)) - A[i - 1];
		}
	}
}

void LinearAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes) {

	if (F.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	D.assign(numsteps, std::vector<double>(nodes, 0.0));
	V.assign(numsteps, std::vector<double>(nodes, 0.0));
	A.assign(numsteps, std::vector<double>(nodes, 0.0));

	std::vector<std::vector<double>> K_ = M * (6.0 / (delta_t * delta_t)) + C * (3.0 / delta_t) + K;
	std::vector<double> P_;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			D[i] = D_i;
			V[i] = V_i;
			A[i] = invertMatrix(M) * (F[i] - K * D[i] - C * V[i]);
		}
		else {
			P_ = F[i] + M * (D[i - 1] * (6.0 / (delta_t * delta_t)) + V[i - 1] * (6.0 / delta_t) + A[i - 1] * 2.0) + C * (D[i - 1] * (3.0 / delta_t) + V[i - 1] * 2.0 + A[i-1] * (delta_t / 2.0));
			D[i] = invertMatrix(K_) * P_;
			V[i] = (D[i] - D[i - 1]) * (3.0 / delta_t) - V[i - 1] * 2.0 - A[i-1] * (delta_t / 2.0);
			A[i] = (D[i] - D[i - 1] - V[i - 1] * delta_t) * (6.0 / (delta_t * delta_t)) - A[i-1] * 2.0;
		}
	}
}