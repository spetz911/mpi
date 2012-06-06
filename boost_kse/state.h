#ifndef __STATE_H__
#define __STATE_H__

extern double rho_g;
extern double rhoU_g;
extern double rhoV_g;
extern double rhoE_g;
extern double gamma_g;

#include <math.h>

struct State {
	double p; // rho
	double u; // rhoU
	double v; // rhoV
	double e; // rhoE

public:
	State()
		: p(rho_g),
		  u(rhoU_g),
		  v(rhoV_g),
		  e(rhoE_g)
	{}
	
	
	double
	calc_eps() const
	{
		return 2.0 * e / (u*u + v*v);
	}

	static
	const State
	update(const State &st)
	{
		State tmp;
		tmp.u = st.u / st.p;
		tmp.v = st.v / st.p;
		tmp.e = st.e / st.p;
		tmp.p = (gamma_g - 1.0) * st.p * tmp.calc_eps();
		return tmp;
	}

	double
	calc_velocity() const
	{
		return sqrt(gamma_g * (gamma_g - 1.0) * this->calc_eps()) + sqrt(u * u + v * v);
	}


};

#endif
