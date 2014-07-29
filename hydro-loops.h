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
double KB, KB_FACT;
double KAPPA_0_E, KAPPA_0_I;
double Z_AVG;
double MI;
double M_EL;
double Q_E;
double PI;

//Declare input structure
struct Options {
	double L;
	double Sh;
	double T0;
	double n0;
	double h0;
	double f_thresh;
	int N;
	int rad_key;
	int heat_key;
};

//Declare structure to return loop data
struct hydroloops_st {
	double *Fe;
	double *Fi;
	double *Te;
	double *Ti;
	double *Pe;
	double *Pi;
	double *n;
	double *s;
	double *r;
	double *g;
	double *h;
	double c2;
	double c3;
	double flux_end;
	double t_end;
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
void hydroloops_calc_abundance(void);

//Declare function hydroloops_avg_val of type double
double hydroloops_avg_val(double[], int);

//Declare function hydroloops_collision_freq of type double
double ebtel_collision_freq(double, double);

#endif