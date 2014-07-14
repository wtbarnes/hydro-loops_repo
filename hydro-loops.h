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
	
};

//Function Prototypes
double * hydroloops_linspace(double, double, int);

struct hydroloops_st *hydroloops_fconverge(double, struct Options);

double hydroloops_heating(double, double, struct Options);

double hydroloops_rad_loss(double, struct Options);

void hydroloops_free_struct(struct hydroloops_st *);

void hydroloops_print_data(struct hydroloops_st *,struct Options);

void hydroloops_print_header(struct Options);

void hydroloops_calc_abundance(char *);

double hydroloops_avg_val(double[], int);

#endif