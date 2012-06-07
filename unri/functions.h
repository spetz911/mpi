#ifndef SOLVE_FUNCTION
#define SOLVE_FUNCTION
#include "structures.h"

#include <math.h>

void solve(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, int flag);
void firstVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt,int print);
void secondVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt, int print);
void thirdVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt);
void fourthVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt);


#endif

