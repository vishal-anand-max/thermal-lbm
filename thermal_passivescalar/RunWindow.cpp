#include "D2Q9Thermal.h"

int main (void)
{
	/*Code for file input
	 *
	int Dimension=2;
	double *values =new double [Dimension+3];
	
	string fileIn="C:/Users/Raja Banerjee/Documents/Research/CFD/LBM/LBM_AAMd/Input.dat";
	
	void fReadFile(double *,int);

	fReadFile(values,Dimension+7);

	D2Q9Thermal twod_flow((int) values[0],(int) values[1], values[2], values[3], values[4]);
	*/
	
	D2Q9Thermal twod_flow(500,500,0.1,1000.0,0.71);

	int no_of_timesteps=50000, autosave=5000, counter=0;
	string filetype="LidDrivenCavity";
	
	//intialization
	twod_flow.init();

	// calling lattice boltzmann method iterations
	for(int j=0; j<=no_of_timesteps; j++)
	{
		twod_flow.equilibriumComputation();
		twod_flow.collide(twod_flow.df,twod_flow.feq,twod_flow.omegaM);
		twod_flow.collide(twod_flow.dg,twod_flow.geq,twod_flow.omegaT);
		twod_flow.stream(twod_flow.df);
		twod_flow.stream(twod_flow.dg);
		twod_flow.boundaryCondition();
		twod_flow.calculateRhoUVPhi();
		if(j==counter*autosave)
		{
			twod_flow.foutput(j,filetype);
			counter++;
		}
		cout<<j<<endl;
	}
	twod_flow.foutput(no_of_timesteps,filetype);
}
/*
void fReadFile(string fileIn,double *values, int n)
{
	ifstream fileInput(fileIn);
	string ReadData;
	
	for(int i=0;i<n;i++)
	{
		do{
			fileInput>>ReadData;
		}while(ReadData.compare(":")!=0);
		fileInput>>values[i];
	}
}
*/