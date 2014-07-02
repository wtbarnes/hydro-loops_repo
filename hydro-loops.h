#ifndef HYDROLOOPS_H
#define HYDROLOOPS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <float.h>

//Declare some global variables
double RSOL;
double GSOL;
double KB;
double KAPPA_0;
double MP;
double MU;
double PI;

//Declare input structure
struct Options {
	double L;
	double Sh;
	double T0;
	double n0;
	int N;
	int rad_key;
};

//Declare structure to return loop data
struct hydroloops_st {
	double *F;
	double *T;
	double *P;
	double *n;
	double *s;
	
};

//Function Prototypes
double * hydroloops_linspace(double, double, double);

struct hydroloops_st *hydroloops_fconverge(double, struct Options);

double hydroloops_heating(double, double, struct Options);

#endif