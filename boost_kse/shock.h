#ifndef __SHOCK_H__
#define __SHOCK_H__

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "field.h"


extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;
extern double gamma_g;

extern int process_id_g;

struct Shock {
	double shock_start;
	double mul;
	int startx;
	int top_osc;
	int bottom_osc;
	double rho0;
	double rho1;
	double p0;
	double p1;
	double u0;
	double u1;
	double eps0;
	double eps1;
	double (*tracef)(double);

public:
	Shock(int nx, double gamma, double rho, double p, double u, double v, double rho_sh, double start_sh, double mul);

	void
	condition(Field &st, double gamma, int partx, int party);
};



#endif
