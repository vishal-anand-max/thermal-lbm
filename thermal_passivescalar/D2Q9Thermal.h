/* Header for Basic (Hydrodynamic) D2Q9 on a Uniform rectangular grid */
#ifndef D2Q9THERMAL_H
#define D2Q9THERMAL_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

class D2Q9Thermal
{
	static const double w[9];
	static const int e[9][2];

public:
	//Values directly from user input
	int Nx; // no. of grid points in x-direction
	int Ny; // no. of points in y-direction
	double viscosity,Re,Pr; //Viscosity, Reynolds Number, Prandtl Number
	double alpha;

	//Values calculated 
	double inletVelocity; //inlet velocity
	double omegaM,omegaT; //Relaxattion paramter (in Non Dimensional case inverse of relaxation time)

	double dx,dy,dt; // temporal and spatial increments (Presently unused)

	typedef struct {
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double f5;
		double f6;
		double f7;
		double f8;
	}distfunc; // distribution function

	double **rho,**ux,**uy, **phi;
	distfunc **df,**feq, **dg, **geq;
	

public:
	D2Q9Thermal(int nx, int ny, double nu, double re, double pr);

	void init();

	void collide(distfunc** f, distfunc** eq, double omega);
	void stream(distfunc**);

	void equilibriumComputation();
	void boundaryCondition();
	void calculateRhoUVPhi();

	void foutput(int time_step,string filetype);
};
#endif