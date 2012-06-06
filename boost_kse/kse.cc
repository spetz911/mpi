#include "fvm.h"
#include "kse.h"



prev_layer::prev_layer(size_t by)
{
	this->by = by;
	prev_row.resize(by);
}


bounds::bounds(size_t bx, size_t by)
{
	this->bx = bx;
	this->by = by;
	
	top_bound.resize(by);
	bottom_bound.resize(by);
	tmpY.resize(by);
	left_bound.resize(bx);
	right_bound.resize(bx);
	tmpX.resize(bx);
}


void
bounds::init(const std::vector< std::vector<double> > &val)
{
	for (int i = 1; i < by + 1; ++i) {
		top_bound[i-1]    = val[1][i];
		bottom_bound[i-1] = val[bx][i];
	}
	
	for (int i = 1; i < bx + 1; ++i) {
		left_bound[i-1]   = val[i][1];
		right_bound[i-1]  = val[i][by];
	}
}

void
bounds::copy(std::vector< std::vector<double> > &val) const
{
	for (int i = 1; i < bx + 1; ++i) {
		val[i][0]      = left_bound[i-1];
		val[i][by + 1] = right_bound[i-1];
	}
	
	for (int i = 1; i < by + 1; ++i) {
		val[0][i]      = top_bound[i-1];
		val[bx + 1][i] = bottom_bound[i-1];
	}
}



void set_partition(int *partx, int *party, int bx, int by)
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
	int max_sect = process_id_g;

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


void set_step_partition(int x_part, int y_part, int *bx, int *by, int nx, int ny)
{
	int rank = process_id_g;

  //number of points for each process
  if ((rank + 1) % y_part == 0)
    *by = ny / y_part + ny % y_part;
  else
    *by = ny / y_part;
  
  if (rank >= y_part * (x_part - 1))
    *bx = nx / x_part + nx % x_part;
  else
    *bx = nx / x_part;
}



/*
void print_res(double *u, int partx, int party, int bx, int by)
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


