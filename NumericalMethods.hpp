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
#include <algorithm>
using namespace std;

//REDO - perhaps use Bogacki–Shampine method and store the error in there too 
vector<complex<double>> RK4_step(vector<complex<double>>compute_f(const double& t, vector<complex<double>>& input), 
	double& h, const double& t, vector<complex<double>> &input_y ) {

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
	k2 = compute_f(t + 0.5*h, y);	

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



//Cash-Karp - 5th order. We are adjusting the step size in this pass of the algorithm to accomadate a maximum error 
//http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf
vector<complex<double>> BSRK4_step_adaptive(vector<complex<double>>compute_f(const double& t, vector<complex<double>>& input), 
	double& h, const double& t, vector<complex<double>> &input_y,double &max_error) {

	//determine number of dimensions
	int n = input_y.size();
	double h_new = h;
	vector<complex<double>> y(n),y2(n),				//the new ys to be evaluated at 4th and 5th order respectively 
							f(n),						//store variable for compute f
							k1(n), k2(n), k3(n), k4(n),k5(n),k6(n);	//coefficients used in algorithim
	
	//get the initial coefficient
	k1 = compute_f(t, input_y);
	
	//second coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.5 *h* k1[i];
	k2 = compute_f(t + 0.5*h, y);	

	//third coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.75*h*k2[i];
	k3 = compute_f(t + 0.75*h,y);
	
	//4th coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + (2.0/9.0)*h*k1[i]+(1.0/3.0)*h*k2[i]+(4.0/9.0)*h*k3[i];
	k4 = compute_f(t + 0.6*h,y);
		
	
	//compute new value for the input_y(5th order) * y2(4th order)
	double tolerance = 0.001;
	for (int i = 0; i < n; i++){
		y2[i]	   = input_y[i] + (2.0/9.0)*h*k1[i] + (1.0/4.0)*h*k2[i] + (4.0/9.0)*h*k3[i];
		input_y[i] = input_y[i] + (7.0/24.0)*h*k1[i] + (1.0/4.0)*h*k2[i] + (1.0/3.0)*h*k3[i] + (1.0/8.0)*h*k4[i];
		//compute the error
		max_error = abs(y2[i] - input_y[i]);
		h_new = 0.9*h*(tolerance/max_error);
	}

	//adjust size as required
	h = min(h_new,h);

	return input_y;
}

//Cash-Karp - 5th order. We are adjusting the step size in this pass of the algorithm to accomadate a maximum error 
//http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf
vector<complex<double>> RK4_step_adaptive(vector<complex<double>>compute_f(const double& t, vector<complex<double>>& input), 
	double& h, const double& t, vector<complex<double>> &input_y,double &error) {

	//determine number of dimensions
	int n = input_y.size();
	double h_new = h;
	vector<complex<double>> y(n),y2(n),				//the new ys to be evaluated at 4th and 5th order respectively 
							f(n),						//store variable for compute f
							k1(n), k2(n), k3(n), k4(n),k5(n),k6(n);	//coefficients used in algorithim
	
	//get the initial coefficient
	k1 = compute_f(t, input_y);
	
	//second coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.2 *h* k1[i];
	k2 = compute_f(t + 0.2*h, y);	

	//third coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + (3.0/40.0)*h*k1[i]+(9.0/40.0) *h*k2[i];
	k3 = compute_f(t + 0.3*h,y);
	
	//4th coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + 0.3*h*k1[i]-0.9*h*k2[i]+1.2*h*k3[i];
	k4 = compute_f(t + 0.6*h,y);
		
	//5th coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] -(11.0/54.0)*h*k1[i]+2.5*h*k2[i]-(70.0/27.0)*h*k3[i]+(35.0/27.0)*h*k4[i];
	k5 = compute_f(t + h,y);
		
	//6th coefficient
	for (int i=0; i < n; i++)
		y[i] = input_y[i] + (1631.0/55296.0)*h*k1[i]+(175.0/512.0)*h*k2[i]-(575.0/13824.0)*h*k3[i]+(44275.0/110592.0)*h*k4[i]+(253.0/409.0)*h*k5[i];
	k6 = compute_f(t + (7.0/8.0)*h,y);
	
	//
	double tolerance = 5e-123;
	for (int i = 0; i < n; i++){
		y2[i]	   = input_y[i] + (2825.0/27648.0)*h*k1[i] + (18575.0/48384.0)*h*k3[i] +(13525.0/55296.0)*h*k4[i] + (277.0/14336.0)*h*k5[i] + 0.25*h*k6[i];
		input_y[i] = input_y[i] + (37.0/378.0)*h*k1[i] + (256.0/621.0)*k3[i] + (125.0/594.0)*h*k4[i] + (512.0/1771.0)*h*k6[i];
		
		//compute the error
		error = abs(y2[i] - input_y[i]);
		h_new = 0.9*h*(tolerance/error);
	}

	//adjust size as required
	h = h_new;

	return input_y;
}

#endif