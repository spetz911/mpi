#include "functions.h"
#include "const.h"
#include <math.h>
#include <stdio.h>

double fmax(double a, double b)
{
	if (a > b)
		return a;
	else return b;
}

void solve(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, int flag)
{
	double dt = T / ( K +1.0 );

  //	printf ("start_i = %d, end_i = %d, start_j = %d, end_j = %d \n", start_i, end_i, start_j, end_j);
	// define border cell or not!
	int print = 0;
	for (int i = 1; i < cell.bsX-1; i++)
	{
		for (int j = 1; j < cell.bsY-1; j++)
		{
			cell.i = i;
			cell.j = j;


			if ((flag == 1) && (i == 7) && (j == 6))
			{
				print = 1;
			}
			else
			{
				print = 0;
			}

			firstVector(prev_layer, cur_layer, cell, dx, dy, dt, 0);
			secondVector(prev_layer, cur_layer, cell, dx, dy, dt, print );
			thirdVector(prev_layer, cur_layer, cell, dx, dy, dt);
			fourthVector(prev_layer, cur_layer, cell, dx, dy, dt);
		
		}	
	}


}
/*

void firstVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt, int print)
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;

	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);

	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i+1].speedU * prev_layer[j*bsX + i+1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].ro		// F(U_i) - F(U_i+1)
		- ciphj * (prev_layer[j*bsX + i+1].ro - prev_layer[j*bsX + i].ro) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].ro		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].ro - prev_layer[j*bsX + i-1].ro) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[(j+1)*bsX + i].speedV * prev_layer[(j+1)*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j+1)*bsX + i].ro - prev_layer[j*bsX + i].ro) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( prev_layer[(j-1)*bsX + i].speedV * prev_layer[(j-1)*bsX + i].ro + prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijmh * (prev_layer[j*bsX + i].ro - prev_layer[(j-1)*bsX + i].ro) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value1 = prev_layer[j*bsX + i].value1 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;

	if (print)
	{
		printf("\n\n prev_layer[j*bsX + i].ro = %f", prev_layer[j*bsX + i].ro );
		printf("\n\n prev_layer[(j-1)*bsX + i].ro = %f", prev_layer[(j-1)*bsX + i].ro );
		printf("\n\n SV - ciphj = %f", ciphj);
		printf("\n\n SV - Fiphj = %f", Fiphj);
		printf("\n\n SV - cimhj = %f", cimhj);
		printf("\n\n SV - Fimhj = %f", Fimhj);
		printf("\n\n SV - cijph = %f", cijph);
		printf("\n\n SV - Gijph = %f", Gijph);
		printf("\n\n SV - cijmh = %f", cijmh);
		printf("\n\n SV - Gijmh = %f", Gijmh);
		printf("\n\n SV - E = %f", prev_layer[j*bsX + i].value4);
		printf("\n j = %d  i = %d\n", j, i);
	}
}


void secondVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt, int print)
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;
	/// Second element of Vector
		// 1)	C_{i+1/2, j}
		double ciphj;	// c i plus half, j
		ciphj = fmax (
			( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
			( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
			);
		
		// 2)	F_{i+1/2, j}
		double Fiphj;	// F i plus half, j
		Fiphj = 0.5 * ( pow(prev_layer[j*bsX + i].speedU, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure + pow(prev_layer[j*bsX + i+1].speedU, 2.0) * prev_layer[j*bsX + i+1].ro + prev_layer[j*bsX + i+1].pressure		// F(U_i) - F(U_i+1)
			- ciphj*(prev_layer[j*bsX + i+1].value2 - prev_layer[j*bsX + i].value2) );
			
		// 3)	C_{i-1/2, j}
		double cimhj;	// c i minus half, j
		cimhj = fmax (
			( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
			( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
			);
		
		// 4)	F_{i-1/2, j}
		double Fimhj;	// F i minus half, j
		Fimhj = 0.5 * ( pow(prev_layer[j*bsX + i-1].speedU, 2.0) * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i-1].pressure + pow(prev_layer[j*bsX + i].speedU, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure		// F(U_i-1) - F(U_i)
			- cimhj*(prev_layer[j*bsX + i].value2 - prev_layer[j*bsX + i-1].value2) );
		
		// 5)	C_{i, j+1/2}
		double cijph;	// c i, j plus half
		cijph = fmax (
			( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
			( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
			);
	
		// 6)	G_{i, j+1/2}
		double Gijph;	// G i, j plus half
		Gijph = 0.5 * ( prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[(j+1)*bsX + i].speedU * prev_layer[(j+1)*bsX + i].speedV * prev_layer[(j+1)*bsX + i].ro		// G(U_i-1) - G(U_i)
			- cijph*(prev_layer[(j+1)*bsX + i].value2 - prev_layer[j*bsX + i].value2) );
		
		// 7)	C_{i, j-1/2}
		double cijmh;	// c i, j minus half
		cijmh = fmax (
			( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
			( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
			);
		
		// 8)	G_{i, j-1/2}
		double Gijmh;	// G i, j minus half
		Gijmh = 0.5 * ( prev_layer[(j-1)*bsX + i].speedU * prev_layer[(j-1)*bsX + i].speedV * prev_layer[(j-1)*bsX + i].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// G(U_i-1) - G(U_i)
			- cijmh*(prev_layer[j*bsX + i].value2 - prev_layer[(j-1)*bsX + i].value2) );
		
		// 9)	U_{i,j}^{k+1}
		cur_layer[j*bsX + i].value2 = prev_layer[j*bsX + i].value2 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;



		if (print)
		{
			printf("\n\n SV - ciphj = %f", ciphj);
			printf("\n\n SV - Fiphj = %f", Fiphj);
			printf("\n\n SV - cimhj = %f", cimhj);
			printf("\n\n SV - Fimhj = %f", Fimhj);
			printf("\n\n SV - cijph = %f", cijph);
			printf("\n\n SV - Gijph = %f", Gijph);
			printf("\n\n SV - cijmh = %f", cijmh);
			printf("\n\n SV - Gijmh = %f", Gijmh);
			printf("\n\n SV - E = %f", prev_layer[j*bsX + i].value4);
			printf("\n j = %d  i = %d\n", j, i);
		}
}

void thirdVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt)// get Third element of Vector
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;
	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);
	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i+1].speedU * prev_layer[j*bsX + i+1].speedV * prev_layer[j*bsX + i+1].ro		// F(U_i) - F(U_i+1)
		- ciphj*(prev_layer[j*bsX + i+1].value3 - prev_layer[j*bsX + i].value3) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * prev_layer[j*bsX + i-1].speedV * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].value3 - prev_layer[j*bsX + i-1].value3) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( pow(prev_layer[j*bsX + i].speedV, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) * prev_layer[(j+1)*bsX + i].ro + prev_layer[(j+1)*bsX + i].pressure		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j+1)*bsX + i].value3 - prev_layer[j*bsX + i].value3) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) * prev_layer[(j-1)*bsX + i].ro + prev_layer[(j-1)*bsX + i].pressure + pow(prev_layer[j*bsX + i].speedV, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure		// G(U_i-1) - G(U_i)
		- cijmh*(prev_layer[j*bsX + i].value3 - prev_layer[(j-1)*bsX + i].value3) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value3 = prev_layer[j*bsX + i].value3 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;

}


void fourthVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt)
/// Fourth element of Vector
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;

	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);

	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i].speedU * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure) + prev_layer[j*bsX + i+1].speedU * (prev_layer[j*bsX + i+1].value4 + prev_layer[j*bsX + i+1].pressure)		// F(U_i) - F(U_i+1)
		- ciphj*(prev_layer[j*bsX + i+1].value4 - prev_layer[j*bsX + i].value4) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * (prev_layer[j*bsX + i-1].value4 + prev_layer[j*bsX + i-1].pressure) + prev_layer[j*bsX + i].speedU * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure)		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].value4 - prev_layer[j*bsX + i-1].value4) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( prev_layer[j*bsX + i].speedV * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure) + prev_layer[(j+1)*bsX + i].speedV * (prev_layer[(j+1)*bsX + i].value4 + prev_layer[(j+1)*bsX + i].pressure)		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j+1)*bsX + i].value4 - prev_layer[j*bsX + i].value4) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( prev_layer[(j-1)*bsX + i].speedV * (prev_layer[(j-1)*bsX + i].value4 + prev_layer[(j-1)*bsX + i].pressure) + prev_layer[j*bsX + i].speedV * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure)		// G(U_i-1) - G(U_i)
		- cijmh*(prev_layer[j*bsX + i].value4 - prev_layer[(j-1)*bsX + i].value4) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value4 = prev_layer[j*bsX + i].value4 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;
	
	
}
*/


void firstVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt, int print)
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;

	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);

	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i+1].speedU * prev_layer[j*bsX + i+1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].ro		// F(U_i) - F(U_i+1)
		- ciphj * (prev_layer[j*bsX + i+1].ro - prev_layer[j*bsX + i].ro) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].ro		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].ro - prev_layer[j*bsX + i-1].ro) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[(j-1)*bsX + i].speedV * prev_layer[(j-1)*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j-1)*bsX + i].ro - prev_layer[j*bsX + i].ro) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( prev_layer[(j+1)*bsX + i].speedV * prev_layer[(j+1)*bsX + i].ro + prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijmh * (prev_layer[j*bsX + i].ro - prev_layer[(j+1)*bsX + i].ro) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value1 = prev_layer[j*bsX + i].value1 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;

	if (print)
	{
		printf("\n\n prev_layer[j*bsX + i].ro = %f", prev_layer[j*bsX + i].ro );
		printf("\n\n prev_layer[(j-1)*bsX + i].ro = %f", prev_layer[(j-1)*bsX + i].ro );
		printf("\n\n SV - ciphj = %f", ciphj);
		printf("\n\n SV - Fiphj = %f", Fiphj);
		printf("\n\n SV - cimhj = %f", cimhj);
		printf("\n\n SV - Fimhj = %f", Fimhj);
		printf("\n\n SV - cijph = %f", cijph);
		printf("\n\n SV - Gijph = %f", Gijph);
		printf("\n\n SV - cijmh = %f", cijmh);
		printf("\n\n SV - Gijmh = %f", Gijmh);
		printf("\n\n SV - E = %f", prev_layer[j*bsX + i].value4);
		printf("\n j = %d  i = %d\n", j, i);
	}
}


void secondVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt, int print)
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;
	/// Second element of Vector
	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);

	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( pow(prev_layer[j*bsX + i].speedU, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure + pow(prev_layer[j*bsX + i+1].speedU, 2.0) * prev_layer[j*bsX + i+1].ro + prev_layer[j*bsX + i+1].pressure		// F(U_i) - F(U_i+1)
		- ciphj*(prev_layer[j*bsX + i+1].value2 - prev_layer[j*bsX + i].value2) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( pow(prev_layer[j*bsX + i-1].speedU, 2.0) * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i-1].pressure + pow(prev_layer[j*bsX + i].speedU, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].value2 - prev_layer[j*bsX + i-1].value2) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[(j-1)*bsX + i].speedU * prev_layer[(j-1)*bsX + i].speedV * prev_layer[(j-1)*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j-1)*bsX + i].value2 - prev_layer[j*bsX + i].value2) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( prev_layer[(j+1)*bsX + i].speedU * prev_layer[(j+1)*bsX + i].speedV * prev_layer[(j+1)*bsX + i].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// G(U_i-1) - G(U_i)
		- cijmh*(prev_layer[j*bsX + i].value2 - prev_layer[(j+1)*bsX + i].value2) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value2 = prev_layer[j*bsX + i].value2 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;



	if (print)
	{
		printf("\n\n prev_layer[j*bsX + i].ro = %f", prev_layer[j*bsX + i].ro );
		printf("\n\n prev_layer[(j-1)*bsX + i].ro = %f", prev_layer[(j+1)*bsX + i].ro );
		printf("\n\n SV - ciphj = %f", ciphj);
		printf("\n\n SV - Fiphj = %f", Fiphj);
		printf("\n\n SV - cimhj = %f", cimhj);
		printf("\n\n SV - Fimhj = %f", Fimhj);
		printf("\n\n SV - cijph = %f", cijph);
		printf("\n\n SV - Gijph = %f", Gijph);
		printf("\n\n SV - cijmh = %f", cijmh);
		printf("\n\n SV - Gijmh = %f", Gijmh);
		printf("\n\n SV - E = %f", prev_layer[j*bsX + i].value4);
		printf("\n j = %d  i = %d\n", j, i);
	}
}

void thirdVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt)// get Third element of Vector
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;
	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);
	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i+1].speedU * prev_layer[j*bsX + i+1].speedV * prev_layer[j*bsX + i+1].ro		// F(U_i) - F(U_i+1)
		- ciphj*(prev_layer[j*bsX + i+1].value3 - prev_layer[j*bsX + i].value3) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * prev_layer[j*bsX + i-1].speedV * prev_layer[j*bsX + i-1].ro + prev_layer[j*bsX + i].speedU * prev_layer[j*bsX + i].speedV * prev_layer[j*bsX + i].ro		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].value3 - prev_layer[j*bsX + i-1].value3) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( pow(prev_layer[j*bsX + i].speedV, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) * prev_layer[(j-1)*bsX + i].ro + prev_layer[(j-1)*bsX + i].pressure		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j-1)*bsX + i].value3 - prev_layer[j*bsX + i].value3) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) * prev_layer[(j+1)*bsX + i].ro + prev_layer[(j+1)*bsX + i].pressure + pow(prev_layer[j*bsX + i].speedV, 2.0) * prev_layer[j*bsX + i].ro + prev_layer[j*bsX + i].pressure		// G(U_i-1) - G(U_i)
		- cijmh*(prev_layer[j*bsX + i].value3 - prev_layer[(j+1)*bsX + i].value3) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value3 = prev_layer[j*bsX + i].value3 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;

}


void fourthVector(TCell *prev_layer, TCell *cur_layer, Cell cell, double dx, double dy, double dt)
/// Fourth element of Vector
{
	int i = cell.i;
	int j = cell.j;
	int bsX = cell.bsX;

	// 1)	C_{i+1/2, j}
	double ciphj;	// c i plus half, j
	ciphj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0 ) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i+1].speedU, 2.0) + pow(prev_layer[j*bsX + i+1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i+1].eps))
		);

	// 2)	F_{i+1/2, j}
	double Fiphj;	// F i plus half, j
	Fiphj = 0.5 * ( prev_layer[j*bsX + i].speedU * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure) + prev_layer[j*bsX + i+1].speedU * (prev_layer[j*bsX + i+1].value4 + prev_layer[j*bsX + i+1].pressure)		// F(U_i) - F(U_i+1)
		- ciphj*(prev_layer[j*bsX + i+1].value4 - prev_layer[j*bsX + i].value4) );

	// 3)	C_{i-1/2, j}
	double cimhj;	// c i minus half, j
	cimhj = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0) * prev_layer[j*bsX + i].eps)),
		( sqrt( pow(prev_layer[j*bsX + i-1].speedU, 2.0) + pow(prev_layer[j*bsX + i-1].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i-1].eps))
		);

	// 4)	F_{i-1/2, j}
	double Fimhj;	// F i minus half, j
	Fimhj = 0.5 * ( prev_layer[j*bsX + i-1].speedU * (prev_layer[j*bsX + i-1].value4 + prev_layer[j*bsX + i-1].pressure) + prev_layer[j*bsX + i].speedU * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure)		// F(U_i-1) - F(U_i)
		- cimhj*(prev_layer[j*bsX + i].value4 - prev_layer[j*bsX + i-1].value4) );

	// 5)	C_{i, j+1/2}
	double cijph;	// c i, j plus half
	cijph = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j-1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j-1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j-1)*bsX + i].eps))
		);

	// 6)	G_{i, j+1/2}
	double Gijph;	// G i, j plus half
	Gijph = 0.5 * ( prev_layer[j*bsX + i].speedV * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure) + prev_layer[(j-1)*bsX + i].speedV * (prev_layer[(j-1)*bsX + i].value4 + prev_layer[(j-1)*bsX + i].pressure)		// G(U_i-1) - G(U_i)
		- cijph*(prev_layer[(j-1)*bsX + i].value4 - prev_layer[j*bsX + i].value4) );

	// 7)	C_{i, j-1/2}
	double cijmh;	// c i, j minus half
	cijmh = fmax (
		( sqrt( pow(prev_layer[j*bsX + i].speedU, 2.0) + pow(prev_layer[j*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[j*bsX + i].eps) ),
		( sqrt( pow(prev_layer[(j+1)*bsX + i].speedU, 2.0) + pow(prev_layer[(j+1)*bsX + i].speedV, 2.0) ) + sqrt (gamma1*(gamma1 - 1.0)*prev_layer[(j+1)*bsX + i].eps))
		);

	// 8)	G_{i, j-1/2}
	double Gijmh;	// G i, j minus half
	Gijmh = 0.5 * ( prev_layer[(j+1)*bsX + i].speedV * (prev_layer[(j+1)*bsX + i].value4 + prev_layer[(j+1)*bsX + i].pressure) + prev_layer[j*bsX + i].speedV * (prev_layer[j*bsX + i].value4 + prev_layer[j*bsX + i].pressure)		// G(U_i-1) - G(U_i)
		- cijmh*(prev_layer[j*bsX + i].value4 - prev_layer[(j+1)*bsX + i].value4) );

	// 9)	U_{i,j}^{k+1}
	cur_layer[j*bsX + i].value4 = prev_layer[j*bsX + i].value4 - dt * (Fiphj - Fimhj) / dx - dt * (Gijph - Gijmh) / dy;


}

