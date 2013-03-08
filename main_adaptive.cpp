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
#define OPEN			false				//defines whether space is open(true) or closed(false), that is if the boundaries loop
											//round to each other or held at zero respectively.
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
double Potential(double t, double x);				//potential field that the electron is in


int main(){
	//define output file
	ofstream output_file, error_file;
	output_file.precision(5); //output to 5 sig fig
	error_file.precision(5);
	error_file.open("C:\\Users\\Jacob\\Documents\\Theory Computing\\Schroedinger WE\\test2\\error.txt");

	//define time step and range  
	double t = 0;										//current place in time
	double T = 1;										//conversion to dimesionless time 
	double dt = 5e-4;									//time step
	double no_steps = 200000;								//60000gives a nice complete picture of what happens
						
	//wavefunction
	vector<complex<double>> psi;
	psi = psi_initial();	
	int output_mod = 10;
	double error = 0.0;

	//perform runge-kutta iterations and output results
	for(int i = 0; i < no_steps; i++){
		//open the file
		stringstream name;
		name.precision(3);

		//output previous iteration
		if(i%output_mod == 0){ //every n steps
			name << "C:\\Users\\Jacob\\Documents\\Theory Computing\\Schroedinger WE\\test2\\time" << i << ".txt";
			output_file.open(name.str());	
			for(int j = 0; j < no_dimensions;j++) //output the probability density 
				output_file  << j << "\t" << abs(psi[j])*abs(psi[j]) << "\t" << (Potential(t, j*grid_width)/100.0) << "\n";

			output_file.close();
		}//close if i%

		//iterate over a step
		//psi = BSRK4_step_adaptive(ModifiedHamiltonian,dt,t,psi,error);
		psi = RK4_step(ModifiedHamiltonian,dt,t,psi);
		//error_file << error << "\t" << dt << "\t" << t << "\t" << endl;
		//increase the time
		t += dt;
	}//end for
	//close the error file
	error_file.close();
	return 0;
}//close main

//For use in RK4 - what the time derivative of wave the function is equal to.
vector<complex<double>> ModifiedHamiltonian(const double& t,vector<complex<double>>& Psi_In){
	//returned variable
	vector<complex<double>> Psi_Out(Psi_In.size());

	//coefficients 
	double r = grid_width/a; 
	complex<double> C = imag_unit/(2*r*r);

	//cycle through the wavefunction, and work out the second space derivative and potential 
	for(int i = 0; i < no_dimensions;i++){	//excludes boundaries - i.e 
		if(OPEN){
			Psi_Out[i] = C*( (i == no_dimensions-1 ? Psi_In[0] : Psi_In[i+1]) - 2.0*Psi_In[i] +(i == 0 ? Psi_In[no_dimensions-1] : Psi_In[i-1]))
				-  (imag_unit)*Potential(t,i*grid_width)*Psi_In[i];	//looped space
		}else{
			Psi_Out[i] =  C*( (i == no_dimensions-1 ? 0.0 : Psi_In[i+1]) - 2.0*Psi_In[i] +(i == 0 ? 0.0 : Psi_In[i-1]))
				-  (imag_unit)*Potential(t,i*grid_width)*Psi_In[i]; //non looped space
		}
	}

	return Psi_Out;
}//close Modified Hamiltonian

//the intial distribution of the wavefunction
vector<complex<double>> psi_initial(){
	//the returned wavefunction
	vector<complex<double>> psi(no_dimensions);

	//constants
	double wdth = 0.8;											//width of the gaussian
	complex<double> total_width = (no_dimensions*grid_width)/a;	//total width of the system
	complex<double> x_0 = total_width/3.0;						//inital displacement from origin
	complex<double> Amplitude = 0.8;							//normalization perhaps?..
	complex<double> p_0		  = 15.0;							//momentum

	//create and fill the vector with approriate values	
	for(int i = 0; i < no_dimensions;i++){
		complex<double> x  = complex<double>((i*grid_width)/a);
		psi[i] = Amplitude*exp(-0.5*(x-x_0)*(x-x_0)/(wdth*wdth))*exp(imag_unit*p_0*(x-x_0));
	}//close for

	if(!OPEN){
		psi[0]				 = 0.0;
		psi[no_dimensions-1] = 0.0;
	}
	return psi;
}//close psi_initial

//defines the real potential field that the system is in.
double Potential(double t, double x){
	if(x >= 7*a
		&& x <= 9*a)
		return 90;
	//return 0.99;
	//return 5*(((x-(grid_width*no_dimensions/2.0))*(x-(grid_width*no_dimensions/2.0)))/((grid_width*no_dimensions/(2.0))*(grid_width*no_dimensions/2.0)));
	return 0;	
}//close Potential	