#ifndef __FVM_H__
#define __FVM_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kse.h"

#define PI 3.14159265
//#define MAX(a, b) ((a) > (b)) ? (a) : (b)


void init_state(state*, double, double, double);
void execute(state*, prev_layer*, double);
void evalVar(state*, prev_layer*, int, int, int, double, double *, double *);
void evalQ(state*, prev_layer*, int, int, int, double, double*, double*);

extern double gl_rho;
extern double gl_rhoU;
extern double gl_rhoV;
extern double gl_rhoE;

#endif
