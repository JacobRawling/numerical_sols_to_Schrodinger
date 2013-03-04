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

//REDO
vector<complex<double>> RK4_step(vector<complex<double>>compute_f(const double& t, vector<complex<double>>& input), 
	const double& h, const double& t, vector<complex<double>> &input_y ) {

	//determine number of dimensions
	int n = input_y.size();
	vector<complex<double>> y(n),						//the new y to be evaluated 
							f(n),						//store variable for compute f
							k1(n), k2(n), k3(n), k4(n);	//coefficients used in RK4
	
	//get the initial coefficient
	k1 = compute_f(t, input_y);
	
	//second coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.5 * h *k1[i];
	k2 = compute_f(t + 05*h, y);	

	//third coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.5 * h * k2[i];
	k3 = compute_f(t + 0.5*h,y);
	
	
	//last coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + h*k3[i];
	k4 = compute_f(t + h,y);
	
	//compute new value for y
	for (int i = 0; i < n; i++)
		input_y[i] = input_y[i] + h*(k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

	return input_y;
}
#endif