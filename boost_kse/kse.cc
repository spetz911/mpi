#include "fvm.h"
#include "kse.h"

state *create_state(int bx, int by)
{
	int i;
	state *st = (state*) malloc(sizeof(state));
	
	st->bx = bx;
	st->by = by;
	
	st->rho  = (double**) malloc(sizeof(double*) * (bx + 2));
	st->rhoU = (double**) malloc(sizeof(double*) * (bx + 2));
	st->rhoV = (double**) malloc(sizeof(double*) * (bx + 2));
	st->rhoE = (double**) malloc(sizeof(double*) * (bx + 2));

	for (i = 0; i < bx + 2; ++i) {
		st->rho[i]  = (double*) malloc(sizeof(double) * (by + 2));
		st->rhoU[i] = (double*) malloc(sizeof(double) * (by + 2));
		st->rhoV[i] = (double*) malloc(sizeof(double) * (by + 2));
		st->rhoE[i] = (double*) malloc(sizeof(double) * (by + 2));
	}

	//printf("1:rho_init-%d\n", sizeof(st->rho[0]));
	return st;
}

void clear_state(state *st)
{
	int i;
	
	for (i = 0; i < st->bx + 2; ++i) {
		free(st->rho[i]);
		free(st->rhoU[i]);
		free(st->rhoV[i]);
		free(st->rhoE[i]);
	}

	free(st->rho);
	free(st->rhoU);
	free(st->rhoV);
	free(st->rhoE);
	free(st);
}

prev_layer *create_prev_layer(int by)
{
	prev_layer *pl = (prev_layer*) malloc(sizeof(prev_layer));
	
	pl->by = by;
	pl->rho_prev_row  = (double*) malloc(sizeof(double) * by);
	pl->rhoU_prev_row = (double*) malloc(sizeof(double) * by);
	pl->rhoV_prev_row = (double*) malloc(sizeof(double) * by);
	pl->rhoE_prev_row = (double*) malloc(sizeof(double) * by);
	
	return pl;
}

void init_prev_layer(prev_layer *pl)
{
	
}

void clear_prev_layer(prev_layer *pl)
{
	free(pl->rho_prev_row);
	free(pl->rhoU_prev_row);
	free(pl->rhoV_prev_row);
	free(pl->rhoE_prev_row);
	
	free(pl);
}


bound *create_bound(int bx, int by)
{
	bound *b = (bound*) malloc(sizeof(bound));
	
	b->bx = bx;
	b->by = by;
	
	b->top_bound    = (double *) malloc(sizeof(double) * by);
	b->bottom_bound = (double *) malloc(sizeof(double) * by);
	b->tmpY         = (double *) malloc(sizeof(double) * by);
	
	b->left_bound   = (double *) malloc(sizeof(double) * bx);
	b->right_bound  = (double *) malloc(sizeof(double) * bx);
	b->tmpX         = (double *) malloc(sizeof(double) * bx);
	
	return b;
}


void init_bound(bound *b, double **val)
{
	int i;

	for (i = 1; i < b->by + 1; ++i) {
		b->top_bound[i-1]    = val[1][i];
		b->bottom_bound[i-1] = val[b->bx][i];
	}
	
	for (i = 1; i < b->bx + 1; ++i) {
		b->left_bound[i-1]   = val[i][1];
		b->right_bound[i-1]  = val[i][b->by];
	}
}

void copy_bound(bound *b, double **val)
{
	int i;

	for (i = 1; i < b->bx + 1; ++i) {
		val[i][0]         = b->left_bound[i-1];
		val[i][b->by + 1] = b->right_bound[i-1];
	}
	
	for (i = 1; i < b->by + 1; ++i) {
		val[0][i]         = b->top_bound[i-1];
		val[b->bx + 1][i] = b->bottom_bound[i-1];
	}
}

void clear_bound(bound *b)
{
	free(b->bottom_bound);
	free(b->left_bound);
	free(b->right_bound);
	free(b->tmpX);
	free(b->tmpY);
	
	free(b);
}

void set_partition(int *partx, int *party, int bx, int by, int num_section)
/*
 * set partition of region
 * example:
 *    partx = 2
 *    party = 2
 *    *------------*------------*
 *    |   fisrt    |   second   |
 *    |   block    |   block    |
 *    *------------*------------|
 *    |   third    |   fourth   |
 *    |   block    |   block    |
 *    *------------*------------*
 */
{
	int max_sect = num_section;

	while(max_sect % 2 == 0) {
		if(bx / *partx >= by / *party)
			*partx <<= 1;
		else
			*party <<= 1;

		max_sect >>= 1;

		continue;
	}

	if(bx / *partx >= by / *party)
		*partx *= max_sect;
	else
		*party *= max_sect;

	//error is not possible
	if(*partx > bx || *party > by)
		exit(0);
}


void set_step_partition(int x_part, int y_part, int *bx, int *by, int nx, int ny, int rank)
{
  //number of points for each process
  if((rank + 1) % y_part == 0)
    *by = ny / y_part + ny % y_part;
  else
    *by = ny / y_part;
  
  if(rank >= y_part * (x_part - 1))
    *bx = nx / x_part + nx % x_part;
  else
    *bx = nx / x_part;
}

double top_boundary_condition(int l, double x, int size)
{
	switch (l) {
	case 0:
		return rho_g;
	case 1:
		return rhoU_g;
	case 2:
		return rhoV_g;
	case 3:
		return rhoE_g;
	};
}

double bottom_boundary_condition(int l, double x, int size)
{
	switch (l) {
	case 0:
		return rho_g;
	case 1:
		return rhoU_g;
	case 2:
		return rhoV_g;
	case 3:
		return rhoE_g;
	};
}

double left_boundary_condition(int l, double y, int size)
{
	switch (l) {
	case 0:
		return rho_g;
	case 1:
		return rhoU_g;
	case 2:
		return rhoV_g;
	case 3:
		return rhoE_g;
	};
}

double right_boundary_condition(int l, double y, int size)
{
	switch (l) {
	case 0:
		return rho_g;
	case 1:
		return rhoU_g;
	case 2:
		return rhoV_g;
	case 3:
		return rhoE_g;
	};
}

void show_result(int bx, int by, double **val, const char *fname, const char *var)
{
	int i, j;
	FILE *F = fopen(fname, "w");

	fprintf(F, "%s variable %d %d\n-------------------------------------\n", var, bx, by);

	for (i = 1; i <= bx; ++i) {
		for (j = 1; j <= by; ++j)
			fprintf(F, "%f ", val[i][j]);

		fprintf(F, "\n");
		//fprintf("\n-----   line %d    --------\n", i);
	}

	fclose(F);
}

/*
void print(double *u, int partx, int party, int bx, int by)
{
	
  int i, j, k, l, inc = bx * by;
	FILE *fp;
	fp = fopen("result","w");
	
	for (l = 0; l < party; ++l)
	  for (i = 0; i < by; ++i)
	    for (j = 0; j < partx; ++j)
	      for (k = 0; k < bx; ++k)
		fprintf(fp, "%f", u[k + j * inc + i * bx + l * inc * partx]);
	
	fclose(fp);
}*/

double psqrt(double x) 
{
	int i = 0;
	double sum = 1.0;
	double mul = x;
	double coeff = 2.0;
	double sig = 1.0;
	
	for (i = 0; i < 10; ++i) {
		sum   += sig * mul / coeff;
		mul   *= x;
		coeff *= 2.0; 
		sig = -sig;
	}
	
	return sum - 1.0;
}
 
double psin(double x) 
{
	int i = 0;
	double sum = 1.0;
	double mul = x;
	double coeff = 1.0;
	double sig = 1.0;
	
	for (i = 0; i < 10; ++i) {
		sum   += sig * mul / coeff;
		mul   *= x * x;
		coeff *= (i + 2) * (i + 3); 
		sig = -sig;
	}
	
	return sum - 1.0;
}

double pline(double x)
{
  return 1.0;
}
