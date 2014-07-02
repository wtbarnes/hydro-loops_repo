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
	
	/****Declare variables****/
	//Double
	double L;
	double Sh;
	double Emin;
	double Emax;
	double T0;
	double n0;
	double f_thresh;
	double f_test;
	double Eh;
	
	//Int
	int N;
	int heat_key;
	int rad_key;
	int Eh_length = 10;
	int i;
	int res_count = 0;
	int res_thresh = 10;
	
	//Pointers
	double *Eh_ptr;
	
	//Char
	char filename_in[64];
	
	//Struct
	struct Options inputs;
	struct hydroloops_st *loop_params;
	
	//File
	FILE *in_file;
	
	/****Read in necessary parameters****/
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
	sprintf(filename_in,"hydro-loops_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error: Could not open file.\n");
		return 1;
	}
	
	fscanf(in_file,"%d\n%d\n%d\n%le\n%le\n%le\n%le\n%le\n",&N,&heat_key,&rad_key,&Emin,&Emax,&T0,&n0,&h0,&f_thresh);
	
	//DEBUG--print the parameters read-in from the file
	printf("Print the input parameters as a test\n");
	printf("%d\t%d\t%d\t%le\t%le\t%le\t%le\t%le\n",N,heat_key,rad_key,Emin,Emax,T0,n0,f_thresh);
	
	//Add necessary inputs to input structure
	inputs.L = L;
	inputs.Sh = Sh;
	inputs.N = N;
	inputs.heat_key = heat_key;
	inputs.rad_key = rad_key;
	inputs.T0 = T0;
	inputs.n0 = n0;
	inputs.h0 = h0;
	
	/****Convergence on flux boundary condition****/
	//Set f_test to start while loop
	f_test = f_thresh + 1;
	
	//Create initial heating array
	Eh_ptr = hydroloops_linspace(Emin,Emax,Eh_length);
	
	while(fabs(f_test) > f_thresh && res_count <= res_thresh)
	{
		for(i = 0; i<Eh_length; i++)
		{
			//Set the heating from the array
			Eh = *(Eh_ptr + i);
			
			//Print the current heating rate
			fprintf("Using heating rate %le\n",Eh);
			
			//Call the flux convergence function 
			loop_params = hydroloops_fconverge(Eh,inputs);
			
			//Set the flux boundary condition
			f_test = *(loop_params->F + (N-1));
			
			//Check the value of f_test to see if it has passed zero(+ threshold flux)
			//(add in a check for complex numbers)
			if(f_test > f_thresh)
			{
				//Reset the heating array
				free(Eh_ptr);
				Eh_ptr = NULL;
				Emax = Eh;
				Emin = *(Eh_ptr + (i-1));
				Eh_ptr = hydroloops_linspace(Emin,Emax,Eh_length);
				
				//Increment the counter
				res_count += 1;
				
				//Break out of the loop 
				break;
			}
			else if(fabs(f_test) > f_thresh)
			{
				//If this condition, we break and then force the while loop to exit
				break;
			}
			
			//Clear the hydroloops_st structure and all of its members. It will be malloc'd 
			//on the next iteration.
		}
		
	}
	
	//If the reset count was exceeded, don't print any data
	if(res_count == res_thresh)
	{
		printf("Unable to converge on desired boundary conditions.\n");
		printf("F(s=L) = %f, Eh = %f for L = %f, Sh = %f\n",*(loop_params->F + (N-1)),Eh,L,Sh);
		pinttf("Trying increasing the grid size or refining the heating array.\n");
	}
	else
	{
		//Print the data to a file
	}
	
	//Clear the memory of the hydroloops structure
	
	
	//Exit with no errors
	return 0;
}