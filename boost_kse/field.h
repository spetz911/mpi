#ifndef __FIELD_H__
#define __FIELD_H__

extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;
extern double gamma_g;

#include <vector>
#include <string>

#include "state.h"

/*
 * field state of the flat gas dynamic system
 */
struct Field {
	std::vector< std::vector<State> > data;
	size_t bx;
	size_t by;
	double hx;
	double hy;
	double ht;
	
public:
	Field(int bx, int by);
	void
	init(double hx, double hy, double ht);

	std::vector<State>&
	operator[](size_t i)
	{
		return data[i];
	}

	const std::vector<State>&
	operator[](size_t i) const
	{
		return data[i];
	}
	
	void
	show_result(const char *fname);
};


#endif
