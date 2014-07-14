#include "hydro-loops.h"

struct hydroloops_st *hydroloops_fconverge(double Eh0, struct Options inputs)
{
	/****Declare variables****/
	//Double
	double s;
	double r;
	double g;
	double h;
	double F;
	double T;
	double n;
	double P;
	double delta_s;
	double dF,dT,dP;
	double lambda;
	double Eh;
	double c2,c3,mean_T;
	
	//Int
	int i;
	
	//Structures
	//Reserve memory for structure that will return loop parameters
	struct hydroloops_st *loop_params = malloc(sizeof(struct hydroloops_st));
	
	//Reserve memory for each structure member
	loop_params->s = malloc(sizeof(double[inputs.N]));
	loop_params->r = malloc(sizeof(double[inputs.N]));
	loop_params->g = malloc(sizeof(double[inputs.N]));
	loop_params->h = malloc(sizeof(double[inputs.N]));
	loop_params->F = malloc(sizeof(double[inputs.N]));
	loop_params->T = malloc(sizeof(double[inputs.N]));
	loop_params->n = malloc(sizeof(double[inputs.N]));
	loop_params->P = malloc(sizeof(double[inputs.N]));
	
	
	/****Initial Calculations****/
	//Calculate spatial step
	delta_s = inputs.L/inputs.N;
	
	//IC for Plasma Properties
	s = inputs.h0;
	r = RSOL + inputs.h0;
	g = GSOL*pow(RSOL/r,2.)*cos(PI*s/(2.*inputs.L));
	h = inputs.h0;
	F = 0.;
	T = inputs.T0;
	n = inputs.n0;
	P = 2*inputs.n0*KB*inputs.T0;
	
	//Save the data for the first iteration
	loop_params->s[0] = s;
	loop_params->r[0] = r;
	loop_params->g[0] = g;
	loop_params->h[0] = h;
	loop_params->F[0] = F;
	loop_params->T[0] = T;
	loop_params->n[0] = n;
	loop_params->P[0] = P;
	
	/****Compute Parameters over Loop Half-Length****/
	for(i = 1; i < inputs.N; i++)
	{
		//Check if temperature is less than zero
		if(T < 0)
		{
			printf("Temperature less than zero. Breaking the loop\n");
			printf("Flux: F = %f\n",F);
			break;
		}
		//Update loop coordinates
		s += delta_s;
		r = RSOL + 2*inputs.L/PI*sin(PI*s/(2.*inputs.L));
		g = GSOL*pow(RSOL/r,2.)*cos(PI*s/(2.*inputs.L));
		h = r - RSOL;
		
		//Calculate the radiative loss function
		lambda = hydroloops_rad_loss(T,inputs);
		
		//Calculate the heating
		Eh = hydroloops_heating(s,Eh0,inputs);
		
		//Step through F,T,P and n parameters using energy and momentum equations
		dF = (Eh - pow(n,2)*lambda)*delta_s;
		dT = F/(-KAPPA_0*pow(T,5./2.))*delta_s;
		dP = -MI*n*g*delta_s;
		
		//Update all of the parameters
		F = dF + F;
		T = dT + T;
		P = dP + P;
		n = P/(2.*KB*T);
		
		//Save updated parameters to data structure
		loop_params->s[i] = s;
		loop_params->r[i] = r;
		loop_params->g[i] = g;
		loop_params->h[i] = h;
		loop_params->F[i] = F;
		loop_params->T[i] = T;
		loop_params->n[i] = n;
		loop_params->P[i] = P;
	}
	
	//Display some results
	printf("******************************************\n");
	printf("BC: F(s=L) = %f\n",F);
	printf("T(s=L) = %f MK\n",T/1e+6);
	printf("n(s=L) = %f (10^8)\n",n/1e+8);
	printf("******************************************\n");
	
	//Calculate the EBTEL coefficients
	mean_T = hydroloops_avg_val(loop_params->T,inputs.N);
	c2 = mean_T/loop_params->T[inputs.N-1];
	c3 = loop_params->T[0]/loop_params->T[inputs.N-1];
	
	//Save these to the structure
	loop_params->c2 = c2;
	loop_params->c3 = c3;
	
	//Return structure containing plasma properties
	return loop_params;
}

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
		Eh = Eh0*exp(-s/inputs.Sh);
	}
	else
	{
		//Heating unspecified
		printf("Type of heating unspecified. Please choose uniform or non-uniform heating\n");
		exit(0);
	}
	
	return Eh;
}

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
		fprintf(out_file,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",*(loop_param->s + i),*(loop_param->r + i),*(loop_param->g + i),*(loop_param->h + i),*(loop_param->F + i),*(loop_param->T + i),*(loop_param->P + i),*(loop_param->n + i));
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
	fprintf(out_file,"%f\n%f\n",loop_param->c2,loop_param->c3);
	
	//Close the file
	fclose(out_file);
}

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


void hydroloops_free_struct(struct hydroloops_st *loop_param)
{
	//Free all structure members
	
	//First check that the pointer structure is valid	
	assert(loop_param != NULL);
	
	//Now free individual members
	free(loop_param->F);
	loop_param->F = NULL;
	free(loop_param->T);
	loop_param->T = NULL;
	free(loop_param->P);
	loop_param->P = NULL;
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

void hydroloops_calc_abundance(char *species)
{
	double k_b = 1.38e-16;
	double m_p = 1.67e-24;
	double m_p_old = m_p;
	
	//Calculate average ion mass
    double n_he_n_p = 0.075;   //He/p abundance.
    Z_AVG = (1.0 + 2.0*n_he_n_p)/(1.0 + n_he_n_p); //Include Helium
    //Z_AVG = 1.; //For Hydrad comparison.
    double kb_fact = 0.5*(1.0+1.0/Z_AVG);
    double m_fact = (1.0 + n_he_n_p*4.0)/(2.0 + 3.0*n_he_n_p); //Include Helium
    //double m_fact = (1 + n_he_n_p*4.)/2.; //For Hydrad comparison
	
	//Set global variables
	if(strcmp(species,"electron") == 0)
	{
		MI = 9.11e-28;
		KB = k_b;
		MU = 1;
		KAPPA_0 = 7.8e-7;		//Spitzer coefficient for thermal conduction (electrons)
	}
	else if(strcmp(species,"ion") == 0)
	{
	    MI = m_p*m_fact*(1.0 + Z_AVG)/Z_AVG; 	//Average ion mass
		KB = k_b*kb_fact; 						//Modify equation of state for non-e-p plasma
		MU = MI/m_p_old;						//Mean molecular weight
		KAPPA_0 = 3.2e-8;						//Spitzer coefficient for thermal conduciton (ions)		
	}
	else
	{
		printf("Incorrect species input. Choose either electron or ion. Exiting the program\n");
		exit(1);
	}
}

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
