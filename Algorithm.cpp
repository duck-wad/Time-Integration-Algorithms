#include <iostream>
#include <cmath>

#include "Algorithm.h"

//implementation of constant, average, and linear acceleration methods of the general Newmark-Beta framework

//constant acceleration for a single degree of freedom system
void ConstantAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.resize(numsteps, 0.0);
	v.resize(numsteps, 0.0);
	a.resize(numsteps, 0.0);

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
		}
		else {
			d[i] = a[i - 1] / 2 * delta_t * delta_t + v[i - 1] * delta_t + d[i - 1];
			v[i] = a[i - 1] * delta_t + v[i - 1];
		}
		a[i] = 1 / m * (f[i] - k * d[i] - c * v[i]);
	}
}

//average acceleration method for SDOF
void AverageAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.resize(numsteps, 0.0);
	v.resize(numsteps, 0.0);
	a.resize(numsteps, 0.0);

	double k_ = (4 / (delta_t * delta_t) * m) + (2 / delta_t * c) + k;
	double p_ = 0.0;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
			a[i] = 1 / m * (f[i] - k * d[i] - c * v[i]);
		}
		else {
			p_ = f[i] + m * (4 / (delta_t * delta_t) * d[i - 1] + 4 / delta_t * v[i - 1] + a[i - 1])
				+ c * (2 / delta_t * d[i - 1] + v[i - 1]);
			d[i] = p_ / k_;
			v[i] = 2 / delta_t * (d[i] - d[i - 1]) - v[i - 1];
			a[i] = 4 / (delta_t * delta_t) * (d[i] - d[i - 1] - delta_t * v[i - 1]) - a[i - 1];
		}
	}
}

//linear acceleration method for SDOF
void LinearAccelerationMethod(std::vector<double>& d, std::vector<double>& v, std::vector<double>& a,
	const std::vector<double>& f, double d_i, double v_i, int numsteps, double delta_t, double m, double c, double k) {

	if (f.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	d.resize(numsteps, 0.0);
	v.resize(numsteps, 0.0);
	a.resize(numsteps, 0.0);

	double k_ = (6 / (delta_t * delta_t) * m) + (3 / delta_t * c) + k;
	double p_ = 0.0;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			d[i] = d_i;
			v[i] = v_i;
			a[i] = 1 / m * (f[i] - k * d[i] - c * v[i]);
		}
		else {
			p_ = f[i] + m * (6 / (delta_t * delta_t) * d[i - 1] + 6 / delta_t * v[i - 1] + 2 * a[i - 1])
				+ c * (3 / delta_t * d[i - 1] + 2 * v[i - 1] + delta_t / 2 * a[i-1]);
			d[i] = p_ / k_;
			v[i] = 3 / delta_t * (d[i] - d[i - 1]) - 2 * v[i - 1] - delta_t / 2 * a[i-1];
			a[i] = 6 / (delta_t * delta_t) * (d[i] - d[i - 1] - delta_t * v[i - 1]) - 2 * a[i-1];
		}
	}
}