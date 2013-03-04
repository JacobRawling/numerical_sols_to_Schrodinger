//=====================================================================
//
//Author: Jacob Rawling
//Date:	05/02/2013
//
//---------------------------------------------------------------------
//
//Aims:
//	1. Use 4th order Runge-Katta to investigate TDSE for arbitrary 
//	   potentials
//	2. Use dimension less units for the wavefunction and time variables. 
//		x->x~ where x~ = x/a. a being the width of the grid spacings.
//		t->t~ where t~ = t/T. T = (ma^2)/(hbar^2)
//	   From here on x and t are refering to x~ and t~, simply for convience.
//
//---------------------------------------------------------------------
//
//This Document:
//	1. Main entry point of programme
//	2. Handling wavefunction as n dimnsional vector, where n is number
//	   of segments of space. 
//
//=====================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>
#include "NumericalMethods.hpp"
using namespace std;

//define the space and it's size
const double a = 1e-10;
const double full_width = 6*a;
const int	 no_dimensions = 401;
const double grid_width = (2*full_width)/no_dimensions;

//physics consts 
const double pi = 3.14159265359; //M_PI  
const double hbar = 1;
const double Me = 1;
const complex<double> imag_unit(0.,1.);

//functions	
vector<complex<double>> psi_initial();	//the initial wavefunctions form	
vector<complex<double>> ModifiedHamiltonian(const double& t,vector<complex<double>>& Psi_In);	//how the wavefunction 
																												//evolves with time
double Potential(double x);				//potential field that the electron is in


int main(){
	//define output file
	ofstream output_file;
	output_file.precision(5); //output to 5 sig fig

	//define time step and range  
	double t = 0;										//current place in time
	double T = 1;										//conversion to dimesionless time 
	double dt = 0.05;									//time step
	double no_steps = 2;

	//wavefunction
	vector<complex<double>> psi;
	psi = psi_initial();	
	
	//perform runge-kutta iterations and output results
	for(int i = 0; i < no_steps; i++){
		//open the file
		stringstream name;
		name.precision(3);

		//output previous iteration
		if(i%1 == 0){ //every n steps
			name << "C:\\Users\\Jacob\\Documents\\Theory Computing\\Schroedinger WE\\test\\time" << i << ".txt";
			output_file.open(name.str());	
			for(int j = 0; j < no_dimensions;j++) //output the probability density 
				output_file  << j << "\t" << abs(psi[j])*abs(psi[j]) << "\t" << Potential(j*grid_width) << "\n";

			output_file.close();
		}//close if i%

		//iterate over a step
		psi = RK4_step(ModifiedHamiltonian,dt,t,psi);
		
		//increase the time
		t += dt;
	}//end for

	return 0;
}//close main

//For use in RK4 - what the time derivative of wave the function is equal to.
vector<complex<double>> ModifiedHamiltonian(const double& t,vector<complex<double>>& Psi_In){
	//returned variable
	vector<complex<double>> Psi_Out = Psi_In;

	//coefficients 
	double r = grid_width/a; 
	complex<double> C = imag_unit/(2*r*r);

	//cycle through the wavefunction, and work out the second space derivative and potential 
	for(int i = 1; i < no_dimensions-1;i++){	//excludes boundaries.
		Psi_Out[i] = C*(Psi_In[i+1] - 2.0*Psi_In[i] + Psi_In[i-1]) - (imag_unit)*Potential(i*grid_width)*Psi_In[i];
	}

	return Psi_Out;
}//close Modified Hamiltonian

//the intial distribution of the wavefunction
vector<complex<double>> psi_initial(){
	//the returned wavefunction
	vector<complex<double>> psi(no_dimensions);

	//constants
	int n = 1;												//principle quantum number
	complex<double> total_width = (no_dimensions*grid_width)/a;	//total width of the system
	complex<double> k = (n*pi)/total_width;					//wave vector 
	complex<double> x_0 = total_width/2.0;					//inital displacement from origin
	complex<double> Amplitude = 1;							//normalization perhaps?..

	//create and fill the vector with approriate values
	for(int i = 0; i < no_dimensions;i++){
		complex<double> x  = complex<double>((i*grid_width)/a);
		psi[i] = Amplitude*exp(-k*(x-x_0)*(x-x_0));//*exp(imag_unit*k*x);
	}//close for

	return psi;
}//close psi_initial

//defines the real potential field that the system is in.
double Potential(double x){
	//return ((x-(grid_width*no_dimensions/2))*(x-(grid_width*no_dimensions/2)));
	return 0;	
}//close Potential	