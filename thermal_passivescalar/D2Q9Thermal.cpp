#include "D2Q9Thermal.h"

using namespace std;

const double D2Q9Thermal::w[9]={4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};
const int D2Q9Thermal::e[9][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

// constructor
D2Q9Thermal::D2Q9Thermal(int nx, int ny, double nu, double re, double pr)
{
	Nx = nx;
	Ny = ny;
	viscosity = nu;
	Re = re;
	Pr = pr;
	alpha = viscosity/Pr;

	inletVelocity = (Re*viscosity)/Ny;
	omegaM = 1.0/(3.0*viscosity+0.5);
	omegaT = 1.0/(3.0*alpha+0.5);

	// Dynamic memory allocation for Density, x-, y-velocity, normalized temperature
	rho = new double *[Nx];
	for(int i=0; i<Nx; i++)
		rho[i] = new double [Ny];

	ux = new double *[Nx];
	for(int i=0; i<Nx; i++)
		ux[i] = new double [Ny];

	uy = new double *[Nx];
	for(int i=0; i<Nx; i++)
		uy[i]=new double [Ny];

	phi = new double *[Nx];
	for(int i=0; i<Nx; i++)
		phi[i]=new double [Ny];

	//Distribution functions for fluid flow
	df = new distfunc *[Nx];
	for(int i=0; i<Nx; i++)
		df[i] = new distfunc [Ny];

	feq = new distfunc *[Nx];
	for(int i=0; i<Nx; i++)
		feq[i]=new distfunc [Ny];

	//Distribution functions for temperature
	dg = new distfunc *[Nx];
	for(int i=0; i<Nx; i++)
		dg[i]=new distfunc [Ny];

	geq = new distfunc *[Nx];
	for(int i=0; i<Nx; i++)
		geq[i]=new distfunc [Ny];
}

// initializing distribution functions for the LBM
void D2Q9Thermal::init()
{
	/*all values initialized to feq computed with
	 * rho = 1.0
	 * ux = 0.0
	 * uy = 0.0
	 * phi = 0.0
	 */

	double a0,a1,a2,a3,a4,a5,a6,a7,a8,u_squared;
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			rho[i][j] = 1.0;
			ux[i][j]  = 0.0;
			uy[i][j]  = 0.0;
			phi[i][j] = 0.0;
			
			a0 = 3.0*(e[0][0]*ux[i][j]+e[0][1]*uy[i][j]);
			a1 = 3.0*(e[1][0]*ux[i][j]+e[1][1]*uy[i][j]);
			a2 = 3.0*(e[2][0]*ux[i][j]+e[2][1]*uy[i][j]);
			a3 = 3.0*(e[3][0]*ux[i][j]+e[3][1]*uy[i][j]);
			a4 = 3.0*(e[4][0]*ux[i][j]+e[4][1]*uy[i][j]);
			a5 = 3.0*(e[5][0]*ux[i][j]+e[5][1]*uy[i][j]);
			a6 = 3.0*(e[6][0]*ux[i][j]+e[6][1]*uy[i][j]);
			a7 = 3.0*(e[7][0]*ux[i][j]+e[7][1]*uy[i][j]);
			a8 = 3.0*(e[8][0]*ux[i][j]+e[8][1]*uy[i][j]);
			u_squared=ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j];
			
			//Density distribution
			df[i][j].f0 = feq[i][j].f0 = w[0]*rho[i][j]*(1.+a0+0.5*a0*a0-1.5*u_squared);
			df[i][j].f1 = feq[i][j].f1 = w[1]*rho[i][j]*(1.+a1+0.5*a1*a1-1.5*u_squared);
			df[i][j].f2 = feq[i][j].f2 = w[2]*rho[i][j]*(1.+a2+0.5*a2*a2-1.5*u_squared);
			df[i][j].f3 = feq[i][j].f3 = w[3]*rho[i][j]*(1.+a3+0.5*a3*a3-1.5*u_squared);
			df[i][j].f4 = feq[i][j].f4 = w[4]*rho[i][j]*(1.+a4+0.5*a4*a4-1.5*u_squared);
			df[i][j].f5 = feq[i][j].f5 = w[5]*rho[i][j]*(1.+a5+0.5*a5*a5-1.5*u_squared);
			df[i][j].f6 = feq[i][j].f6 = w[6]*rho[i][j]*(1.+a6+0.5*a6*a6-1.5*u_squared);
			df[i][j].f7 = feq[i][j].f7 = w[7]*rho[i][j]*(1.+a7+0.5*a7*a7-1.5*u_squared);
			df[i][j].f8 = feq[i][j].f8 = w[8]*rho[i][j]*(1.+a8+0.5*a8*a8-1.5*u_squared);
			
			//Temperature dsitribution
			dg[i][j].f0 = geq[i][j].f0 = w[0]*phi[i][j]*(1.+a0);
			dg[i][j].f1 = geq[i][j].f1 = w[1]*phi[i][j]*(1.+a1);
			dg[i][j].f2 = geq[i][j].f2 = w[2]*phi[i][j]*(1.+a2);
			dg[i][j].f3 = geq[i][j].f3 = w[3]*phi[i][j]*(1.+a3);
			dg[i][j].f4 = geq[i][j].f4 = w[4]*phi[i][j]*(1.+a4);
			dg[i][j].f5 = geq[i][j].f5 = w[5]*phi[i][j]*(1.+a5);
			dg[i][j].f6 = geq[i][j].f6 = w[6]*phi[i][j]*(1.+a6);
			dg[i][j].f7 = geq[i][j].f7 = w[7]*phi[i][j]*(1.+a7);
			dg[i][j].f8 = geq[i][j].f8 = w[8]*phi[i][j]*(1.+a8);
		}
	}
}

// function to calculate the f-equilibrium
void D2Q9Thermal::equilibriumComputation()
{
	double a0, a1, a2, a3, a4, a5, a6, a7, a8, u_squared;

	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			a0 = 3.0*(e[0][0]*ux[i][j]+e[0][1]*uy[i][j]);
			a1 = 3.0*(e[1][0]*ux[i][j]+e[1][1]*uy[i][j]);
			a2 = 3.0*(e[2][0]*ux[i][j]+e[2][1]*uy[i][j]);
			a3 = 3.0*(e[3][0]*ux[i][j]+e[3][1]*uy[i][j]);
			a4 = 3.0*(e[4][0]*ux[i][j]+e[4][1]*uy[i][j]);
			a5 = 3.0*(e[5][0]*ux[i][j]+e[5][1]*uy[i][j]);
			a6 = 3.0*(e[6][0]*ux[i][j]+e[6][1]*uy[i][j]);
			a7 = 3.0*(e[7][0]*ux[i][j]+e[7][1]*uy[i][j]);
			a8 = 3.0*(e[8][0]*ux[i][j]+e[8][1]*uy[i][j]);
			u_squared = ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j];

			//Density Equilibrium
			feq[i][j].f0 = w[0]*rho[i][j]*(1.+a0+0.5*a0*a0-1.5*u_squared);
			feq[i][j].f1 = w[1]*rho[i][j]*(1.+a1+0.5*a1*a1-1.5*u_squared);
			feq[i][j].f2 = w[2]*rho[i][j]*(1.+a2+0.5*a2*a2-1.5*u_squared);
			feq[i][j].f3 = w[3]*rho[i][j]*(1.+a3+0.5*a3*a3-1.5*u_squared);
			feq[i][j].f4 = w[4]*rho[i][j]*(1.+a4+0.5*a4*a4-1.5*u_squared);
			feq[i][j].f5 = w[5]*rho[i][j]*(1.+a5+0.5*a5*a5-1.5*u_squared);
			feq[i][j].f6 = w[6]*rho[i][j]*(1.+a6+0.5*a6*a6-1.5*u_squared);
			feq[i][j].f7 = w[7]*rho[i][j]*(1.+a7+0.5*a7*a7-1.5*u_squared);
			feq[i][j].f8 = w[8]*rho[i][j]*(1.+a8+0.5*a8*a8-1.5*u_squared);

			//Temperature Equilibrium
			geq[i][j].f0 = w[0]*phi[i][j]*(1.+a0);
			geq[i][j].f1 = w[1]*phi[i][j]*(1.+a1);
			geq[i][j].f2 = w[2]*phi[i][j]*(1.+a2);
			geq[i][j].f3 = w[3]*phi[i][j]*(1.+a3);
			geq[i][j].f4 = w[4]*phi[i][j]*(1.+a4);
			geq[i][j].f5 = w[5]*phi[i][j]*(1.+a5);
			geq[i][j].f6 = w[6]*phi[i][j]*(1.+a6);
			geq[i][j].f7 = w[7]*phi[i][j]*(1.+a7);
			geq[i][j].f8 = w[8]*phi[i][j]*(1.+a8);
		}
	}
}

// function to calculate the collision terms
void D2Q9Thermal::collide(distfunc** f, distfunc** eq, double omega)
{
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			f[i][j].f0 -= omega*(f[i][j].f0-eq[i][j].f0);
			f[i][j].f1 -= omega*(f[i][j].f1-eq[i][j].f1);
			f[i][j].f2 -= omega*(f[i][j].f2-eq[i][j].f2);
			f[i][j].f3 -= omega*(f[i][j].f3-eq[i][j].f3);
			f[i][j].f4 -= omega*(f[i][j].f4-eq[i][j].f4);
			f[i][j].f5 -= omega*(f[i][j].f5-eq[i][j].f5);
			f[i][j].f6 -= omega*(f[i][j].f6-eq[i][j].f6);
			f[i][j].f7 -= omega*(f[i][j].f7-eq[i][j].f7);
			f[i][j].f8 -= omega*(f[i][j].f8-eq[i][j].f8);
		}										
	}
}

// function to calculate the streaming terms
void D2Q9Thermal::stream(distfunc** f)
{
	for(int j=0; j<Ny; j++)
	{
		// f1  from right to left
		for(int i=Nx-1; i>0; i--)
			f[i][j].f1 = f[i-e[1][0]][j-e[1][1]].f1;
		
		// f3 from left to right
		for(int i=0; i<Nx-1; i++)
			f[i][j].f3 = f[i-e[3][0]][j-e[3][1]].f3;
	}

	// f2, f5, f6 from top to bottom
	for(int j=Ny-1; j>0; j--)
	{
		// f2
		for(int i=0; i<Nx; i++)
			f[i][j].f2 = f[i-e[2][0]][j-e[2][1]].f2;
		// f5
		for(int i=Nx-1; i>0; i--)
			f[i][j].f5 = f[i-e[5][0]][j-e[5][1]].f5;
		// f6
		for(int i=0; i<Nx-1; i++)
			f[i][j].f6 = f[i-e[6][0]][j-e[6][1]].f6;
	}

	// f4, f7, f8 from bottom to top
	for(int j=0;j<Ny-1; j++)
	{
		// f4
		for(int i=0; i<Nx; i++)
			f[i][j].f4 = f[i-e[4][0]][j-e[4][1]].f4;
		// f7
		for(int i=0; i<Nx-1; i++)
			f[i][j].f7 = f[i-e[7][0]][j-e[7][1]].f7;
		// f8
		for(int i=Nx-1; i>0; i--)
			f[i][j].f8 = f[i-e[8][0]][j-e[8][1]].f8;
	}
}

// function to apply the boundary conditions
void D2Q9Thermal::boundaryCondition()
{
	/*
	 * Fluid boundary conditions 
	 */
	
	/*
	// West boundary (Dirichlett Velcoity)
	double rhoW=1.0, uW=inletVelocity, vW=0.0;
	for(int j=0; j<Ny; j++)
	{
		rhoW = 1./(1-uW)*(df[0][j].f0+df[0][j].f2+df[0][j].f4+
						  2.*(df[0][j].f3+df[0][j].f6+df[0][j].f7));

		df[0][j].f1 = df[0][j].f3+(2./3.)*rhoW*uW;
		df[0][j].f5 = df[0][j].f7-0.5*(df[0][j].f2-df[0][j].f4)+rhoW*uW/6.+0.5*rhoW*vW;
		df[0][j].f8 = df[0][j].f6+0.5*(df[0][j].f2-df[0][j].f4)+rhoW*uW/6.-0.5*rhoW*vW;
	}

	//East boundary (Dirichlett Pressure)
	double rhoE=0.9, uE=0.0, vE=0.0;
	for(int j=0; j<Ny; j++)
	{
		uE = -1.+(df[Nx-1][j].f0+df[Nx-1][j].f2+df[Nx-1][j].f4+
				  2.*(df[Nx-1][j].f1+df[Nx-1][j].f5+df[Nx-1][j].f8))/rhoE;
		
		df[Nx-1][j].f3 = df[Nx-1][j].f1-(2/3)*rhoE*uE;
		df[Nx-1][j].f6 = df[Nx-1][j].f5+0.5*(df[Nx-1][j].f2-df[Nx-1][j].f4)-(1/6)*rhoE*uE;
		df[Nx-1][j].f7 = df[Nx-1][j].f8-0.5*(df[Nx-1][j].f2-df[Nx-1][j].f4)-(1/6)*rhoE*uE;
	}

	//South surface (Wall - Bounce Back)
	for(int i=0; i<Nx; i++)
	{
		df[i][0].f2 = df[i][0].f4;
		df[i][0].f5 = df[i][0].f7;
		df[i][0].f6 = df[i][0].f8;
	}

	//North surface (Wall - Bounce Back)
	double rhoN = 1.0, uN = inletVelocity;
	for(int i=0; i<Nx; i++)
	{
		df[i][Ny-1].f4 = df[i][Ny-1].f2;
		df[i][Ny-1].f7 = df[i][Ny-1].f5;
		df[i][Ny-1].f8 = df[i][Ny-1].f6;
	}
	*/
	
	// West boundary (Wall - Bounce back)
	for(int j=0; j<Ny; j++)
	{
		df[0][j].f1 = df[0][j].f3;
		df[0][j].f5 = df[0][j].f7;
		df[0][j].f8 = df[0][j].f6;
	}

	//East boundary (Wall - Bounce back)
	for(int j=0; j<Ny; j++)
	{
		df[Nx-1][j].f3 = df[Nx-1][j].f1;
		df[Nx-1][j].f7 = df[Nx-1][j].f5;
		df[Nx-1][j].f6 = df[Nx-1][j].f8;
	}

	//South surface (Wall - Bounce back)
	for(int i=0; i<Nx; i++)
	{
		df[i][0].f2 = df[i][0].f4;
		df[i][0].f5 = df[i][0].f7;
		df[i][0].f6 = df[i][0].f8;
	}

	//North surface (Moving Lid)
	double rhoN = 1.0, uN = inletVelocity;
	for(int i=0; i<Nx; i++)
	{
		rhoN = df[i][Ny-1].f0+df[i][Ny-1].f1+df[i][Ny-1].f3+
			   2*(df[i][Ny-1].f2+df[i][Ny-1].f6+df[i][Ny-1].f8);

		df[i][Ny-1].f4 = df[i][Ny-1].f2;
		df[i][Ny-1].f7 = df[i][Ny-1].f5+0.5*(df[i][Ny-1].f1-df[i][Ny-1].f3)-0.5*rhoN*uN;
		df[i][Ny-1].f8 = df[i][Ny-1].f6-0.5*(df[i][Ny-1].f1-df[i][Ny-1].f3)+0.5*rhoN*uN;
	}
	
	// Corners
	double d=1.0; // corner density
	df[0][0].f6=df[0][0].f8=0.5*(d-(df[0][0].f0+df[0][0].f1+df[0][0].f2+df[0][0].f3
								+df[0][0].f4+df[0][0].f5+df[0][0].f7));
	df[Nx-1][0].f5=df[Nx-1][0].f7=0.5*(d-(df[Nx-1][0].f0+df[Nx-1][0].f1+df[Nx-1][0].f2+
									df[Nx-1][0].f3+df[Nx-1][0].f4+df[Nx-1][0].f6+df[Nx-1][0].f8));
	df[Nx-1][Ny-1].f6=df[Nx-1][Ny-1].f8=0.5*(d-(df[Nx-1][Ny-1].f0+df[Nx-1][Ny-1].f1+
									df[Nx-1][Ny-1].f2+df[Nx-1][Ny-1].f3+df[Nx-1][Ny-1].f4+
									df[Nx-1][Ny-1].f5+df[Nx-1][Ny-1].f7));
	df[0][Ny-1].f5=df[0][Ny-1].f7=0.5*(d-(df[0][Ny-1].f0+df[0][Ny-1].f1+df[0][Ny-1].f2+
									df[0][Ny-1].f3+df[0][Ny-1].f4+df[0][Ny-1].f6+df[0][Ny-1].f8));

	/*
	 * Thermal boundary conditions 
	 */
	// West boundary (phiW = 0) 
	double phiW = 0.0;
	for(int j=0; j<Ny; j++)
	{
		dg[0][j].f1 = -dg[0][j].f3;
		dg[0][j].f5 = -dg[0][j].f7;
		dg[0][j].f8 = -dg[0][j].f6;
	}
		
	//East boundary (phiW = 0.0)
	for(int j=0; j<Ny; j++)
	{
		dg[Nx-1][j].f3 = -dg[Nx-1][j].f1;
		dg[Nx-1][j].f7 = -dg[Nx-1][j].f5;
		dg[Nx-1][j].f6 = -dg[Nx-1][j].f8;
	}

	//Bottom surface (Adiabatic)
	for(int i=0; i<Nx; i++)
	{
		dg[i][0].f0 = dg[i][0].f0;
		dg[i][0].f1 = dg[i][0].f1;
		dg[i][0].f2 = dg[i][0].f2;
		dg[i][0].f3 = dg[i][0].f3;
		dg[i][0].f4 = dg[i][0].f4;
		dg[i][0].f5 = dg[i][0].f5;
		dg[i][0].f6 = dg[i][0].f6;
		dg[i][0].f7 = dg[i][0].f7;
		dg[i][0].f8 = dg[i][0].f8;
	}

	//Top surface (Dirichlett Temperature phiT = 1.0)
	double phiT = 1.0;
	for(int i=0; i<Nx; i++)
	{
		dg[i][Ny-1].f4 = phiT*(w[2]+w[4])-dg[i][Ny-1].f2;
		dg[i][Ny-1].f7 = phiT*(w[5]+w[7])-dg[i][Ny-1].f5;
		dg[i][Ny-1].f8 = phiT*(w[6]+w[8])-dg[i][Ny-1].f6;
		
		dg[i][Ny-1].f1 = phiT*(w[1]+w[3])-dg[i][Ny-1].f3;
	}

	// Corners
	double phiC = 0.0; // corner density
	phiC = 0.0;
	dg[0][0].f6=dg[0][0].f8=0.5*(phiC-(dg[0][0].f0+dg[0][0].f1+dg[0][0].f2+dg[0][0].f3+
									   dg[0][0].f4+dg[0][0].f5+dg[0][0].f7));
	dg[Nx-1][0].f5=dg[Nx-1][0].f7=0.5*(phiC-(dg[Nx-1][0].f0+dg[Nx-1][0].f1+dg[Nx-1][0].f2+dg[Nx-1][0].f3+
											 dg[Nx-1][0].f4+dg[Nx-1][0].f6+dg[Nx-1][0].f8));

	phiC = 1.0;
	dg[Nx-1][Ny-1].f6=dg[Nx-1][Ny-1].f8=0.5*(phiC-(dg[Nx-1][Ny-1].f0+dg[Nx-1][Ny-1].f1+dg[Nx-1][Ny-1].f2+
												   dg[Nx-1][Ny-1].f3+dg[Nx-1][Ny-1].f4+dg[Nx-1][Ny-1].f5+
												   dg[Nx-1][Ny-1].f7));
	dg[0][Ny-1].f5=dg[0][Ny-1].f7=0.5*(phiC-(dg[0][Ny-1].f0+dg[0][Ny-1].f1+dg[0][Ny-1].f2+
									         dg[0][Ny-1].f3+dg[0][Ny-1].f4+dg[0][Ny-1].f6+dg[0][Ny-1].f8));
}

//calculate the new field variables from the updated distribution function
void D2Q9Thermal::calculateRhoUVPhi()
{
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			rho[i][j]=df[i][j].f0+df[i][j].f1+df[i][j].f2+df[i][j].f3+df[i][j].f4+
					  df[i][j].f5+df[i][j].f6+df[i][j].f7+df[i][j].f8;
			ux[i][j]=(1./rho[i][j])*(df[i][j].f1-df[i][j].f3+df[i][j].f5-
									 df[i][j].f6-df[i][j].f7+df[i][j].f8);
			uy[i][j]=(1./rho[i][j])*(df[i][j].f2-df[i][j].f4+df[i][j].f5+
									 df[i][j].f6-df[i][j].f7-df[i][j].f8);

			phi[i][j]=dg[i][j].f0+dg[i][j].f1+dg[i][j].f2+dg[i][j].f3+dg[i][j].f4+
					  dg[i][j].f5+dg[i][j].f6+dg[i][j].f7+dg[i][j].f8;
		}
	}
}

// write to file function
void D2Q9Thermal::foutput(int time,string ftype)
{
	ofstream outputfile;
	stringstream ss;
	ss<<"C:/Users/Yashwanth Yadavalli/Dropbox/LBM Solver Thermal/Data/"<<ftype<<"-D2Q9-"<<time<<".csv";
	string filename=ss.str();
	cout<<endl<<filename<<endl;
	outputfile.open(filename,ios::out);
	outputfile<<"TITLE = \"Lid Driven Cavity\""<<endl;
	outputfile<<"variables = x , y, Density, U, V, Velocity, Phi"<<endl;
	outputfile<<"ZONE"<<" i = "<<Nx<<" j = "<<Ny<<" F=POINT "<<endl;
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			double u = sqrt(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]);
			outputfile<<i<<","<<j<<","<<rho[i][j]<<","<<ux[i][j]<<","<<uy[i][j]<<","<<u<<","<<phi[i][j]<<endl;
		}
	}
	outputfile.close();
}