#ifndef __KSE_H__
#define __KSE_H__

#include <stdio.h>
#include <stdlib.h>

/*
 * /ru/ состояние системы плоской газовой динамики
 * /en/ state of the flat gas dynamic system
 */
typedef struct {
	double **rho;
	double **rhoU;
	double **rhoV;
	double **rhoE;
	short unsigned int bx;
	short unsigned int by;
	double hx;
	double hy;
	double ht;
} state;


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
	double rho_old;
	double rhoU_old;
	double rhoV_old;
	double rhoE_old;
	double rho_prev_old;
	double rhoU_prev_old;
	double rhoV_prev_old;
	double rhoE_prev_old;
	double *rho_prev_row;
	double *rhoU_prev_row;
	double *rhoV_prev_row;
	double *rhoE_prev_row;
	unsigned short int by;
} prev_layer;

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

state 	   *create_state(int, int);
void   		clear_state(state*);

prev_layer *create_prev_layer(int);
void        init_prev_layer(prev_layer*);
void        clear_prev_layer(prev_layer*);

bound      *create_bound(int, int);
void     	init_bound(bound*, double**);
void   		copy_bound(bound*, double**);
void        clear_bound(bound*);

void   		set_partition(int*, int*, int, int, int);
void   		set_step_partition(int, int, int*, int*, int, int, int);
double 		top_boundary_condition(int, double, int);
double 		bottom_boundary_condition(int, double, int);
double 		left_boundary_condition(int, double, int);
double 		right_boundary_condition(int, double, int);
void 		show_result(int bx, int by, double**, const char*, const char*);

extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;

double psqrt(double);
double psin(double);
double pline(double);

#endif
