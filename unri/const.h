#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__
#include <math.h>
//const double ro = 1.0;	// plotnost' 
const double gamma1 = 5.0 / 3.0;
//const double pressure = 1.0;


const double mkm = 0.000001;

const double ro1 = 50; // plotnost' in layer 1
const double ro2 = 1000; // plotnost' in layer 2
const double ro3 = 1; // plotnost' in layer 3

const double h1 = 10 * mkm;  // layer 1
const double h2 = 30 * mkm;  // layer 2
const double h3 = 105 * mkm;  // layer 3

const double pressure1 = 1.0E10; // pressure in layer 1;
const double pressure2 = 1.0E10; // pressure in layer 2;
const double pressure3 = 1.0E10; // pressure in layer 3;

const double h = 150 * mkm; // max Y
const double l = 5 * mkm; // max X
const double a0 = 1 * mkm;

const double T = 1E-12 * 1300 * 10;
const int K = 1000 * 10; // num steps
const double M = 2; // Mah num


const int printstep = 50;

const double PI = atan(1.0) * 4.0; 


#endif
