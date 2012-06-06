#ifndef __KSE_H__
#define __KSE_H__

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "state.h"
#include "field.h"

extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;
extern double gamma_g;

extern int process_id_g;

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

struct prev_layer {
	State old;
	State prev_old;
	std::vector<State> prev_row;
	size_t by;
public:
	prev_layer(size_t by);

};


struct bounds {
	std::vector<double> top_bound;
	std::vector<double> bottom_bound;
	std::vector<double> left_bound;
	std::vector<double> right_bound;
	std::vector<double> tmpX;
	std::vector<double> tmpY;
	//number of points for each process
	size_t bx;
	size_t by;

public:
	bounds(size_t bx, size_t by);
	void
	init(const std::vector< std::vector<double> > &val);
	void
	copy(std::vector< std::vector<double> > &val) const;

};



void
set_partition(int*, int*, int, int);
void
set_step_partition(int, int, int*, int*, int, int);


double psqrt(double);
double psin(double);

#endif
