#ifndef __KSE_H__
#define __KSE_H__

#include <stdio.h>
#include <stdlib.h>

extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;
extern double gamma_g;

// eval shock wave in the node
const double mul_g = ;


struct state {
	double p;  // rho
	double u; // rhoU
	double v; // rhoV
	double e; // rhoE

public:

	const state
	operator/ (double x) const
	{
	    state tmp = *this;
	    
	    tmp.rho  /= x;
	    tmp.rhoU /= x;
	    tmp.rhoV /= x;
	    tmp.rhoE /= x;
	    
	    return tmp;
	}

	double
	calc_eps() const
	{
		return (2 * rhoE) / (rhoU * rhoU + rhoV * rhoV);
	}

	static
	const state
	update_state(const state &st, const double &gamma)
	{
		state tmp = st;
		tmp.rhoU /= st.rho;
		tmp.rhoV /= st.rho;
		tmp.rhoE /= st.rho;
		tmp.rho = (gamma_g - 1.0) * st.rho * tmp.calc_eps();
		return tmp;
	}

	double
	calc_velocity() const
	{
		double tmp;
		
		tmp  = sqrt(gamma_g * (gamma_g - 1.0) * this->calc_eps());
		tmp += sqrt(rhoU * rhoU + rhoV * rhoV);
		return tmp;
	}


	

};

/*
 * field state of the flat gas dynamic system
 */
class field {
	state **data;
	short unsigned int bx;
	short unsigned int by;
	double hx;
	double hy;
	double ht;
	
public:

	double *&
	operator[](size_t i)
	{
		return data[i];
	}

};


/* /ru/
 * Для расчета значения характеристики    2
 * в узле необходимы значения в         1 * 3
 * четырех окружающих точках              4
 * При этом значения точек 1 и 4 уже рассчитаны
 * Чтобы получить значения точек 1 и 4 на предыдущем 
 * слое, необходимо их предварительно сохранить.
 * 
 * Значение <var>_old      соответствует текущей точке
 * Значение <var>_prev_old соответствует точке 1
 * Значения <var>_prev_row соответствуют токам 4, т.е. нижней строки 
 * by - длина строки
 */

typedef struct {
	state old;
	state prev_old;
	state *prev_row;
	unsigned short int by;
} prev_layer;


prev_layer

typedef struct {
	double *top_bound;
	double *bottom_bound;
	double *left_bound;
	double *right_bound;
	double *tmpX;
	double *tmpY;
	//number of points for each process
	unsigned short int bx;
	unsigned short int by;
} bound;

typedef struct {
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
} shock;

field 	   *create_field(int, int);
void   		clear_field(field*);

prev_layer *create_prev_layer(int);
void        init_prev_layer(prev_layer*);
void        clear_prev_layer(prev_layer*);

bound      *create_bound(int, int);
void
init_bound(bound*, double**);
void   		copy_bound(bound*, double**);
void        clear_bound(bound*);

void   		set_partition(int*, int*, int, int, int);
void   		set_step_partition(int, int, int*, int*, int, int, int);
double 		top_boundary_condition(int, double, int);
double 		bottom_boundary_condition(int, double, int);
double 		left_boundary_condition(int, double, int);
double 		right_boundary_condition(int, double, int);
void 		show_result(int bx, int by, double**, const char*, const char*);


double psqrt(double);
double psin(double);
double pline(double);

#endif