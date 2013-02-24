//=====================================================================
//
//Author: Jacob Rawling
//Date:	05/02/2013
//
//---------------------------------------------------------------------
//
//Aim:
//	1. Simulate tragectory of a ball under gravity. 
//	2. Use 4th order Runge-Katta to do this. 
//	3. Output tragectory into a file as a table of time and position. 
//	   (t, x, y)
//
//---------------------------------------------------------------------
//
//This Document:
//	1. Define functions for: derivative of a function & one RK4 step
//
//=====================================================================


//header guards
#ifndef NUMERICALMETHODS_HPP
#define NUMERICALMETHODS_HPP

//useful headers
#include <cmath>
#include <vector>
#include <complex>
using namespace std;

vector<complex<double>> RK4_step(vector<complex<double>>(*compute_f)(vector<complex<double>>& input,vector<complex<double>>& output), 
	const double h,	vector<complex<double>> &y_vector ) {
	//determine number of dimensions
	int n = y_vector.size();
	vector<complex<double>> y(n), f(n), k1(n), k2(n), k3(n), k4(n);
	//get the initial function thing 
	compute_f(y_vector, k1);
	
	//determine the coefficients
	for(int i = 0; i < n; i++)
		k1[i] = h * k1[i];
	
	//second coefficient
	for (int i=0; i < n; i++)
		y[i] = y_vector[i] + 0.5 * k1[i];
	compute_f(y, f);	
	for (int i = 0; i < n; i++)
		k2[i] = h * f[i];

	//third coefficient
	for (int i=0; i < n; i++)
		y[i] = y_vector[i] + 0.5 * k2[i];
	compute_f(y, f);
	for (int i = 0; i < n; i++)
		k3[i] = h * f[i];
	
	//last coefficient
	for (int i=0; i < n; i++)
		y[i] = y_vector[i] + k3[i];
	compute_f(y, f);
	for (int i = 0; i < n; i++)
		k4[i] = h * f[i];
	
	//compute new value for y
	for (int i = 0; i < n; i++)
		y[i] = y[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

	return y;
}

#endif