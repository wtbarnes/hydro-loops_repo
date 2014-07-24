#ifndef HYDROLOOPS_H
#define HYDROLOOPS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>

//Declare some global variables
double RSOL;
double GSOL;
double KB;
double KAPPA_0;
double MU;
double Z_AVG;
double MI;
double PI;

//Declare input structure
struct Options {
	double L;
	double Sh;
	double T0;
	double n0;
	double h0;
	int N;
	int rad_key;
	int heat_key;
};

//Declare structure to return loop data
struct hydroloops_st {
	double *F;
	double *T;
	double *P;
	double *n;
	double *s;
	double *r;
	double *g;
	double *h;
	double c2;
	double c3;
	double flux_end;
	
};

//Function Prototypes

//Declare function hydroloops_linspace of type double
double * hydroloops_linspace(double, double, int);

//Declare function hydroloops_fconverge of type struct hydroloops_st
struct hydroloops_st *hydroloops_fconverge(double, struct Options);

//Declare function hydroloops_heating of type double
double hydroloops_heating(double, double, struct Options);

//Declare function hydroloops_rad_loss of type double
double hydroloops_rad_loss(double, struct Options);

//Declare function hydroloops_free_struct of type void
void hydroloops_free_struct(struct hydroloops_st *);

//Declare function hydroloops_print_data of type void
void hydroloops_print_data(struct hydroloops_st *,struct Options);

//Declare function hydroloops_print_header of type void
void hydroloops_print_header(struct Options);

//Declare function hydroloops_calc_abundance of type void
void hydroloops_calc_abundance(char *);

//Declare function hydroloops_avg_val of type double
double hydroloops_avg_val(double[], int);

#endif