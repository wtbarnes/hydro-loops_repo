/***************************************************************************
File name: hydro-loops_functions.c
Author: Will Barnes
Date: 23 July 2014
Description: This file holds the functions for the hydro-loops program which
solves the hydrostatic equations for a uniformly or non-uniformly heated coronal
loop. A description of each function can be found below.
***************************************************************************/

//Include appropriate header file
#include "hydro-loops.h"

/***************************************************************************
Function name: hydroloops_fconverge
Description: This function solves the hydrostatic equations for some given 
heating and inputs as specified in the inputs structure. 
Inputs:
	Eh0: heating rate coefficient (erg cm^-3)
	Options inputs: structure containing all information read in from input
					file. See hydro-loops.h for layout.
Outputs:
	loop_params: structure that holds plasma parameters computed by the 
				 hydrostatic equations and updated loop coordinates.
***************************************************************************/


struct hydroloops_st *hydroloops_fconverge(double Eh0, struct Options inputs)
{
	/****Declare variables****/
	//Double
	double s;
	double r;
	double g;
	double h;
	double Fe,Fi;
	double Te,Ti;
	double n;
	double Pe,Pi;
	double delta_s;
	double dFe,dFi,dTe,dTi,dPe,dPi;
	double lambda;
	double Eh;
	double nu_ei;
	double c2e,c3e,mean_Te;
	double c2i,c3i,mean_Ti;
	
	//Int
	int i;
	int t_flag = 0;
	
	//Structures
	//Reserve memory for structure that will return loop parameters
	struct hydroloops_st *loop_params = malloc(sizeof(struct hydroloops_st));
	
	//Reserve memory for each structure member
	loop_params->s = malloc(sizeof(double[inputs.N]));
	loop_params->r = malloc(sizeof(double[inputs.N]));
	loop_params->g = malloc(sizeof(double[inputs.N]));
	loop_params->h = malloc(sizeof(double[inputs.N]));
	loop_params->Fe = malloc(sizeof(double[inputs.N]));
	loop_params->Te = malloc(sizeof(double[inputs.N]));
	loop_params->Pe = malloc(sizeof(double[inputs.N]));
	loop_params->Fi = malloc(sizeof(double[inputs.N]));
	loop_params->Ti = malloc(sizeof(double[inputs.N]));
	loop_params->Pi = malloc(sizeof(double[inputs.N]));
	loop_params->n = malloc(sizeof(double[inputs.N]));
	
	/****Initial Calculations****/
	//Calculate spatial step
	delta_s = inputs.L/inputs.N;
	
	//IC for Plasma Properties
	s = inputs.h0;
	r = RSOL + inputs.h0;
	g = GSOL*pow(RSOL/r,2.)*cos(PI*s/(2.*inputs.L));
	h = inputs.h0;
	Fe = 0.;
	Fi = Fe;
	Te = inputs.T0;
	Ti = Te;
	n = inputs.n0;
	Pe = inputs.n0*KB*inputs.T0;
	Pi = Pe;
	
	//Save the data for the first iteration
	loop_params->s[0] = s;
	loop_params->r[0] = r;
	loop_params->g[0] = g;
	loop_params->h[0] = h;
	loop_params->Fe[0] = Fe;
	loop_params->Fi[0] = Fi;
	loop_params->Te[0] = Te;
	loop_params->Ti[0] = Ti;
	loop_params->n[0] = n;
	loop_params->Pe[0] = Pe;
	loop_params->Pi[0] = Pi;
	
	/****Compute Parameters over Loop Half-Length****/
	for(i = 1; i < inputs.N; i++)
	{
		//Check if temperature is less than zero
		if(Te < 0 || Ti < 0)
		{
			//Print error to the screen
			printf("Temperature less than zero. Breaking the loop\n");
			printf("Flux: F = %f\n",F);
			
			//Reset the flux to avoid breaking early
			Fe = inputs.f_thresh + 1.;
			
			//Raise the negative temperature flag
			t_flag = 1;
			
			break;
		}
		//Update loop coordinates
		s += delta_s;
		r = RSOL + 2*inputs.L/PI*sin(PI*s/(2.*inputs.L));
		g = GSOL*pow(RSOL/r,2.)*cos(PI*s/(2.*inputs.L));
		h = r - RSOL;
		
		//Calculate the radiative loss function
		lambda = hydroloops_rad_loss(Te,inputs);
		
		//Calculate the heating
		Eh = hydroloops_heating(s,Eh0,inputs);
		
		//Calculate collisional frequency of collisions between electrons and ions
		nu_ei = hydroloops_collision_freq(Te,n);
		
		//Step through F,T,P and n parameters using energy and momentum equations
		dFe = (Eh - pow(n,2)*lambda + 3./2.*KB*n*nu_ei*(Ti - Te))*delta_s;
		dFi = (3./2.*KB*n*nu_ei*(Te - Ti))*delta_s;
		dTe = Fe/(-KAPPA_0_E*pow(Te,5./2.))*delta_s;
		dTi = Fi/(-KAPPA_0_I*pow(Ti,5./2.))*delta_s;
		dPe = -M_EL*n*g*delta_s;
		dPi = -MI*n*g*delta_s;
		
		//Update all of the parameters--FINISH UPDATING FROM HERE FOR TWO FLUID CASE!!
		Fe = dFe + Fe;
		Fi = dFi + Fi;
		Te = dTe + Te;
		Ti = dTi + Ti;
		Pe = dPe + Pe;
		Pi = dPi + Pi;
		n = Pe/(KB*Te);
		
		//Save updated parameters to data structure
		loop_params->s[i] = s;
		loop_params->r[i] = r;
		loop_params->g[i] = g;
		loop_params->h[i] = h;
		
		loop_params->Pe[i] = Pe;
		loop_params->Fe[i] = Fe;
		loop_params->Te[i] = Te;
		loop_params->n[i] = n;
		loop_params->Pi[i] = Pi;
		loop_params->Fi[i] = Fi;
		loop_params->Ti[i] = Ti;
	}
		
	//Display some results
	printf("******************************************\n");
	printf("BC: F(s=L) = %f\n",F);
	printf("T(s=L) = %f MK\n",T/1e+6);
	printf("n(s=L) = %f (10^8)\n",n/1e+8);
	printf("******************************************\n");
	
	//Check and make sure we didn't break early for T<0
	//This avoids segmentation fault if we broke early
	if(t_flag == 0)
	{
		//Calculate the EBTEL coefficients
		mean_Te = hydroloops_avg_val(loop_params->Te,inputs.N);
		mean_Ti = hydroloops_avg_val(loop_params->Ti,inputs.N);
		c2e = mean_Te/loop_params->Te[inputs.N-1];
		c3e = loop_params->Te[0]/loop_params->Te[inputs.N-1];
		c2i = mean_Ti/loop_params->Ti[inputs.N-1];
		c3i = loop_params->Ti[0]/loop_params->Ti[inputs.N-1];
	
		//Save these to the structure
		loop_params->c2e = c2e;
		loop_params->c3e = c3e;
		loop_params->c2i = c2i;
		loop_params->c3i = c3i;
	}
	
	//Set the flux end value
	loop_params->flux_end = Fe;
	
	//Set the temperature end value
	loop_params->t_end = Te;
	
	//Return structure containing plasma properties
	return loop_params;
}

/***************************************************************************
Function name: hydroloops_heating
Description: This function computes the ad-hoc volumetric heating. 
Inputs:
	s: loop coordinate (cm)
	Eh0: heating rate coefficient (erg cm^-3)
	Options inputs: structure containing all information read in from input
					file. See hydro-loops.h for layout.
Outputs:
	Eh: heating rate; for the uniform case, Eh=Eh0
***************************************************************************/

double hydroloops_heating(double s,double Eh0,struct Options inputs)
{
	//Declare variables
	double Eh;
	
	//First check if the heating is uniform or non-uniform
	if(inputs.heat_key == 0)
	{
		//Uniform heating
		Eh = Eh0;
	}
	else if(inputs.heat_key == 1)
	{
		//Non-uniform heating
		Eh = Eh0*exp(-pow((s - inputs.h0),2.)/(2.*pow(inputs.Sh,2.)));
		//Eh = Eh0*exp(-(s - inputs.h0)/inputs.Sh);
	}
	else
	{
		//Heating unspecified
		printf("Type of heating unspecified. Please choose uniform or non-uniform heating\n");
		exit(0);
	}
	
	return Eh;
}

/***************************************************************************
Function name: hydroloops_rad_loss
Description: This function computes the radiative loss function for a given
	 		 temperature.
Inputs:
	T: temperature (K)
	Options inputs: structure containing all information read in from input
					file. See hydro-loops.h for layout.
Outputs:
	lambda: radiative loss function
***************************************************************************/

double hydroloops_rad_loss(double T, struct Options inputs)
{
	//Declare some variables
	double alpha;
	double chi;
	double logT;
	double lambda;
	
	//Take log of temperature
	logT = log10(T);
	
	//Decide on full or simplified loss function
	if(inputs.rad_key == 0)
	{
		//Simplified radiative loss function
		if(logT <= 4.97)
		{
			alpha = 2.0;
			chi = 1.09e-31;
		}
		else
		{
			alpha = 0.5;
			chi = 2.19e-19;
		}
	}
	else if(inputs.rad_key == 1)
	{
		//Full radiative loss function
		if(logT >= 7.63)
		{
			alpha = 0.5;
			chi = 1.96e-27;
		}
		else if(logT >= 6.90)
		{
			alpha = -1.0;
			chi = 5.49e-16;
		}
		else if(logT >= 6.55)
		{
			alpha = 1./3.;
			chi = 3.46e-25;
		}
		else if(logT >= 6.18)
		{
			alpha = -1.5;
			chi = 3.53e-13;
		}
		else if(logT >= 5.67)
		{
			alpha = 0.;
			chi = 1.90e-22;
		}
		else if(logT >= 4.97)
		{
			alpha = -1.0;
			chi = 8.87e-17;
		}
		else
		{
			alpha = 2.0;
			chi = 1.09e-31;
		}
	}
	else
	{
		printf("Invalid input for radiative loss function. Exiting program.\n");
		exit(1);
	}
	
	//Calculate lambda
	lambda = chi*pow(T,alpha);
	
	return lambda;
}

/***************************************************************************
Function name: hydroloops_print_data
Description: This function prints the contents of the main data structure to
a specified file. 
Inputs:
	hydroloops_st loop_param: structure that holds all computed plasma 
							  parameters
	Options inputs: structure that holds input parameters
Outputs:
***************************************************************************/

void hydroloops_print_data(struct hydroloops_st *loop_param, struct Options inputs)
{
	//Declare variables
	int i;
	int i_max = inputs.N;
	char fn_out[64];
	char fn_out_coeff[64];
	FILE *out_file;
	
	//Check and see if directory 'data' exists. If it does not, then create a new one.
	struct stat st = {0};
	if(stat("data",&st) == -1){
		mkdir("data",0777);
	}
	
	//Create the filename and open the file
	sprintf(fn_out,"data/loopsdat_L_%d_Sh_%d_Eh_%d.txt",(int) (inputs.L/1e+8), (int) (inputs.Sh/1e+8), inputs.heat_key);
	out_file = fopen(fn_out,"wt");
	
	//Tell the user where the results were printed
	printf("The results were printed to the file %s\n",fn_out);
	
	//Print the data to the file
	for(i = 0; i<i_max; i++)
	{
		fprintf(out_file,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",*(loop_param->s + i),*(loop_param->r + i),*(loop_param->g + i),*(loop_param->h + i),*(loop_param->Fe + i),*(loop_param->Fi + i),*(loop_param->Te + i),*(loop_param->Ti + i),*(loop_param->Pe + i),*(loop_param->Pi + i),*(loop_param->n + i));
	}
	
	//Close the file when we are done writing to it
	fclose(out_file);
	
	/****Print Coefficient Results****/
	
	//Print relavent coefficients to separate file
	sprintf(fn_out_coeff,"data/coeffdat_L_%d_Sh_%d_Eh_%d.txt",(int) (inputs.L/1e+8), (int) (inputs.Sh/1e+8), inputs.heat_key);
	out_file = fopen(fn_out_coeff,"wt");
	
	//Tell the user where the coefficient results were printed
	printf("The coefficient results were printed to the file %s\n",fn_out_coeff);
	
	//Print the data to the file
	fprintf(out_file,"%f\n%f\n%f\n%f\n",loop_param->c2e,loop_param->c3e,loop_param->c2i,loop_param->c3i);
	
	//Close the file
	fclose(out_file);
}

/***************************************************************************
Function name: hydroloops_print_header
Description: This function prints a header to the screen each time the program
begins execution.
Inputs:
	Options inputs: structure that holds input parameters
Outputs:
***************************************************************************/

void hydroloops_print_header(struct Options inputs)
{
	//Print header to the screen
	printf("************************************************************************************\n");
	printf("                                   HydroLoops                                       \n");
	printf("                         Coronal Hydrostatics Simulation                            \n");
	printf("************************************************************************************\n");
	printf("Code based on ASTR 554, Astrophysics of the Sun\n");
	printf("Course taught by Prof. Stephen Bradshaw, Rice University, Fall 2013\n");
	printf("See also Aschwanden et al. (2001)\n");
	printf("Coded by Will Barnes, Dept. of Physics and Astronomy, Rice University\n");
	printf("************************************************************************************\n");
	printf("Inputs:\n");
	printf("Loop half-length, L = %f Mm\n",inputs.L/1e+8);
	printf("Heating scale height, Sh = %f Mm\n",inputs.Sh/1e+8);
	printf("Number of grid points,  N = %d\n",inputs.N);
	if(inputs.heat_key == 0)
	{
		printf("Using uniform heating.\n");
	}
	else
	{
		printf("Using non-uniform heating.\n");
	}
	if(inputs.rad_key == 0)
	{
		printf("Using simplified radiative loss function.\n");
	}
	else
	{
		printf("Using full radiative loss function.\n");
	}
	printf("Base temperature, T(s=s0) = %f MK\n",inputs.T0/1e+6);
	printf("Base density, n(s=s0) = %f (10^8) cm^-3\n",inputs.n0/1e+8);
	printf("Height above the solar surface, h0 = %f Mm\n",inputs.h0/1e+8);
	printf("************************************************************************************\n");
	
}

/***************************************************************************
Function name: hydroloops_free_struct
Description: This function frees memory used by the main structure to store
plasma parameters computed by the hydrostatic equations.
Inputs:
	hydroloops_st loop_param: structure that holds all computed plasma 
							  parameters
Outputs:
***************************************************************************/

void hydroloops_free_struct(struct hydroloops_st *loop_param)
{
	//Free all structure members
	
	//First check that the pointer structure is valid	
	assert(loop_param != NULL);
	
	//Now free individual members
	free(loop_param->Fe);
	loop_param->Fe = NULL;
	free(loop_param->Te);
	loop_param->Te = NULL;
	free(loop_param->Pe);
	loop_param->Pe = NULL;
	free(loop_param->Fi);
	loop_param->Fi = NULL;
	free(loop_param->Ti);
	loop_param->Ti = NULL;
	free(loop_param->Pi);
	loop_param->Pi = NULL;
	
	free(loop_param->n);
	loop_param->n = NULL;
	free(loop_param->s);
	loop_param->s = NULL;
	free(loop_param->r);
	loop_param->r = NULL;
	free(loop_param->g);
	loop_param->g = NULL;
	free(loop_param->h);
	loop_param->h = NULL;
	
	//Free the structure itself
	free(loop_param);
}

/***************************************************************************
Function name: hydroloops_linspace
Description: This function constructs a vector beginning at a and ending at 
			 b with n entries; meant to emulate equivalen MATLAB function.
Inputs:
	a: first entry of array
	b: last entry of array
	n: number of entries in array
Outputs:
***************************************************************************/

double * hydroloops_linspace( double a, double b, int n)
{
	//Declare necessary variables
	int i;	
	double *linspace = malloc(sizeof(double[n]));
	double N = n;
	double interval = (b - a)/(N-1);	//spacing between points
	
	//Make the array
	linspace[0] = a;
	for(i = 1; i<n; i++)
	{
		linspace[i] = interval + *(linspace + (i-1));
	}

	return linspace;
}

/***************************************************************************
Function name: hydroloops_calc_abundance
Description: This function calculates the mean molecular weight given some abundance
and modifies expressions for KB and MP to account for presence of additional
species. 
Inputs:
	species: specifies whether run is for electrons or ions
Outputs:
***************************************************************************/

void hydroloops_calc_abundance(void)
{
	double m_p = 1.67e-24;
	double m_p_old = m_p;
	
	//Calculate average ion mass
    double n_he_n_p = 0.075;   //He/p abundance.
    double m_fact = (1.0 + n_he_n_p*4.0)/(2.0 + 3.0*n_he_n_p); //Include Helium
    //double m_fact = (1 + n_he_n_p*4.)/2.; //For Hydrad comparison
    Z_AVG = (1.0 + 2.0*n_he_n_p)/(1.0 + n_he_n_p); //Include Helium
    //Z_AVG = 1.; //For Hydrad comparison.
    KB_FACT = 0.5*(1.0+1.0/Z_AVG);
	
	KAPPA_0_E = 7.8e-7;		//Spitzer coefficient for thermal conduction (electrons)
	KAPPA_0_I = 3.2e-8;		//Spitzer coefficient for thermal conduciton (ions)		

    MI = m_p*m_fact*(1.0 + Z_AVG)/Z_AVG; 	//Average ion mass
}

/***************************************************************************
Function name: hydroloops_avg_val
Description: This function computes the average value for some array.
Inputs:
	numbers: array of numbers whose mean is calculated
	length: length of the numbers array
Outputs:
	mean: mean of the numbers array
***************************************************************************/

double hydroloops_avg_val (double numbers[], int length)
{
	//Declare some variables
	int i;
	double mean;
	double sum = 0;
	
	//Calculate the sum
	for (i = 0; i<length; i++)
	{
		sum = sum + numbers[i];	
	}
	
	//Calculate the mean
	mean = sum/length;
	
	//Return the mean
	return mean;
}

/***********************************************************************************

FUNCTION NAME: hydroloops_collision_freq

FUNCTION_DESCRIPTION: This function calculates the coulomb collision frequency between
the ions and electrons. The coulomb logarithm is calculated using the procedure found 
in the 2009 edition of the NRL Plasma Formulary.

INPUTS:
	T_e--electron temperature (K)
	T_i--ion temperature (K)
	n--number density (cm^-3)
	
OUTPUTS:
	nu_ei--electron-ion collision frequency

***********************************************************************************/

double ebtel_collision_freq(double T_e, double n)
{
	//Declare variables
	double ln_lambda;
	double nu_ei;
	double beta_1 = 1.0e+13;
	double beta_2 = 1.602*1e-9;
	
	//Expression for the Coulomb logarithm from Physics of the Solar Corona by M.J. Aschwanden
	ln_lambda = 23 - log(sqrt(n/beta_1)*pow(KB*T_e/beta_2,-3./2.));	
		
	//Calculate collision frequency
	nu_ei = 16./3.*sqrt(PI)*pow(Q_E,4.)/(M_EL*M_P)*pow(2*KB*T_e/M_EL,-3./2.)*n*ln_lambda;
	
	return nu_ei;
}
