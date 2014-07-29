/***************************************************************************
File name: hydro-loops_main.c
Author: Will Barnes
Date: 23 July 2014
Description: This file holds the main function for the hydro-loops program
which solves the hydrstatic equations for a uniformly or non-uniformly heated
coronal loop. This program is based on notes from the ASTR 554 Astrophysics of 
the Sun course taught by Professor Stephen Bradshaw, Rice University, Fall 2013.
***************************************************************************/

#include "hydro-loops.h"

int main(int argc, char *argv[])
{
	//Use clock to time the entire hydro-loops program
	clock_t time_start;
	clock_t time_diff;
	double time_elapsed;
	
	//Start the timer
	time_start = clock();
	
	//Define global variables
	RSOL = 6.9550e+10; 	//Radius of the Sun in cm 
	GSOL = 27395.;		//Gravitational acceleration at the solar surface
	PI = 3.14159;
	Q_E = 4.8032e-10;	//charge of e- in stat coloumbs
	M_EL = 9.11e-28;	//mass of e- in grams
	KB = 1.38e-16;		//Boltzmann constant in erg K^-1
	
	/****Declare variables****/
	//Double
	double L;
	double Sh;
	double Emin;
	double Emax,Emax_old;
	double T0;
	double n0;
	double h0;
	double f_thresh;
	double f_test;
	double t_test;
	double Eh;
	char species[64];
	
	//Int
	int N;
	int heat_key;
	int rad_key;
	int Eh_length = 10;
	int i;
	int res_count = 0;
	int res_thresh = 20;
	
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
		L = 1e+8*L;					//Convert from Mm to cm
		Sh = atof(argv[2]);			//Read in heating scale height
		Sh = 1e+8*Sh;				//Convert from Mm to cm
	}
	
	//Read in parameters from input file
	sprintf(filename_in,"hydro-loops_parameters.txt");
	in_file = fopen(filename_in,"rt");
	if(in_file == NULL)
	{
		printf("Error: Could not open file.\n");
		return 1;
	}
	
	fscanf(in_file,"%d\n%d\n%d\n%le\n%le\n%le\n%le\n%le\n%le\n",&N,&heat_key,&rad_key,&Emin,&Emax,&T0,&n0,&h0,&f_thresh);
	
	//Add necessary inputs to input structure
	inputs.L = L;
	inputs.Sh = Sh;
	inputs.N = N;
	inputs.heat_key = heat_key;
	inputs.rad_key = rad_key;
	inputs.T0 = T0;
	inputs.n0 = n0;
	inputs.h0 = h0;
	inputs.f_thresh = f_thresh;
	
	//Calculate parameters specific to species
	hydroloops_calc_abundance();
	
	/****Print header to standard output****/
	hydroloops_print_header(inputs);
	
	/****Convergence on flux boundary condition****/
	//Set f_test to start while loop
	f_test = f_thresh + 1;
		
	while(fabs(f_test) > f_thresh && res_count <= res_thresh)
	{
		for(i = 0; i<Eh_length; i++)
		{
			//Set the heating from the array
			Eh = Emin + i*(Emax - Emin)/(Eh_length - 1.);
			
			//Print the current heating rate
			printf("Using heating rate %le\n",Eh);
			
			//Call the flux convergence function 
			loop_params = hydroloops_fconverge(Eh,inputs);
			
			//Set the flux boundary condition
			f_test = loop_params->flux_end;
			
			//Set the T condition
			t_test = loop_params->t_end;
			
			//Check the value of f_test to see if it has passed zero(+ threshold flux)
			//(add in a check for complex numbers)
			if(res_count > res_thresh)
			{
				//Make sure we have not exceeded the maximum number of resets
				//If we have, break the for loop so the while loop can be broken
				break;
				
			}
			else if(f_test > f_thresh || t_test < 0.)
			{
				//Reset the heating array
				Emax_old = Emax;
				Emax = Eh;
				Emin = Eh - (Emax_old - Emin)/(Eh_length - 1.);
				
				//Increment the counter
				res_count += 1;
				//Print the counter
				printf("Heating reset %d times\n",res_count);
				
				//Clear the hydroloops_st structure and all of its members. It will be malloc'd on the
				//next iteration.
				hydroloops_free_struct(loop_params);
				
				//Break out of the loop 
				break;
			}
			else if(fabs(f_test) < f_thresh)
			{
				//If this condition, we break and then force the while loop to exit
				//The desired boundary condition has been met to within the set tolerance
				break;
			}
			
			//Clear the hydroloops_st structure and all of its members. It will be malloc'd 
			//on the next iteration.
			hydroloops_free_struct(loop_params);
		}
		
	}
	
	//If the reset count was exceeded, don't print any data
	if(res_count > res_thresh)
	{
		printf("Unable to converge on desired boundary conditions.\n");
		printf("F(s=L) = %4.2le, Eh = %f for L = %4.2le, Sh = %4.2le\n",loop_params->flux_end,Eh,L,Sh);
		printf("Try increasing the grid size or refining the heating array.\n");
		printf("It is also possible that there is no solution for (L,Sh) = (%4.2f,%4.2f) Mm\n",L/1.e+8,Sh/1.e+8);
	}
	else
	{
		//Print the data to a file
		hydroloops_print_data(loop_params,inputs);
		//Clear the memory of the hydroloops structure
		hydroloops_free_struct(loop_params);
	}
	
	//Stop the timer
	time_diff = clock() - time_start;
	time_elapsed = time_diff*1000/CLOCKS_PER_SEC;
	
	//Time elapsed
	printf("The process took %f milliseconds to run\n",time_elapsed);
	
	//Exit with no errors
	return 0;
}