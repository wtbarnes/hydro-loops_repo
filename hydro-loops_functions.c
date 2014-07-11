#include "hydro-loops.h"

double * hydroloops_linspace( double a, double b, double n)
{
	//Declare necessary variables
	int i;	
	double *linspace = malloc(sizeof(double[n]));
	double interval = (b - a)/(n-1);	//spacing between points
	
	//Make the array
	linspace[0] = a;
	for(i = 1; i<n; i++)
	{
		linspace[i] = interval + *(linspace + (i-1));
	}

	return linspace;
}

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
	s = h0;
	r = RSOL + h0;
	g = GSOL*pow(RSOL/r,2.)*cos(PI*s/(2.*inputs.L));
	h = h0;
	F = 0.;
	T = T0;
	n = n0;
	P = 2*n0*KB*T0;
	
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
		//Calculate some inital loop coordinates
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
		dT = F/(-KAPPA_O*pow(T,5./2.))*delta_s;
		dP = -MU*MP*n*g*delta_s;
		
		//Update all of the parameters
		F = dF + F;
		T = dT + T;
		P = dP + P;
		n = P/(2.*K_B*T);
		
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
	printf("Using heating rate E_H = %f\n",Eh0);
	printf("BC: F(s=L) = %f\n",F);
	printf("T(s=L) = %f MK\n",T/1e+6);
	printf("n(s=L) = %f (10^8)\n",n/1e+8);
	printf("******************************************\n");
	
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


void hydroloops_struct_free(struct hydroloops_st *loop_param)
{
	//Free all structure members
	
	//First check that the pointer structure is valid	
	
}