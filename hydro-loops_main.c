#include "hydro-loops.h"

int main(int argc, char *argv[])
{
	//Define global variables
	RSOL = 6.9550e+10; 	//Radius of the Sun in cm 
	GSOL = 27395.;		//Gravitational acceleration at the solar surface
	KB = 1.38e-16;		//Boltzmann constant
	KAPPA_0 = 1e-6;		//Spitzer coefficient for thermal conduction
	MP = 1.67e-24;		//proton mass in g
	MU = 1.;			
	PI = M_PI;
	
	//Declare variables
	double L;
	double Sh;
	double Emin;
	double Emax;
	double T0;
	double n0;
	double f_thresh;
	
	int N;
	int heat_key;
	int rad_key;
	
	char filename_in[64];
	
	//Read in command line arguments
	if(argc != 3)
	{
		printf("Incorrect number of input arguments. Exiting program.\n");
		return 1;
	}
	else
	{
		L = atof(argv[1]);			//Read in loop half-length
		Sh = atof(argv[2]);			//Read in heating scale height
	}
	
	//Read in parameters from input file
	sprintf(filename_in,"hydro-loop_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error: Could not open file.\n");
		return 1;
	}
	
	fscanf(in_file,"%d\n%d\n%d\n%f\n%f\n%f\n%f\n%f\n%f\n",N,heat_key,rad_key,Emin,Emax,T0,n0,f_thresh);
	
	//Exit with no errors
	return 0;
}