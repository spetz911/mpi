#include "fvm.h"

#include <algorithm>

void init_state(state *st, double hx, double hy, double ht)
{
	st->hx = hx;
	st->hy = hy;
	st->ht = ht;
	
	for (int i = 0; i < st->bx + 2; ++i)
		for (int j = 0; j < st->by + 2; ++j) {
			st[i][j]->init();
		}
}


void
evalVar(state *st, prev_layer *pl, int i, int j)
{
	// eval Q
	
	//en/ for each point it's coumpute separately
	//ru/ следующие параметры рассчитываются для каждой точки отдельно
	//  #
	// $ * #
	//  $
	//  |
	//en/ get value for all variables v
	//ru/ получим значения всех переменных в top
	//ru/ 5-ти ближайших точках: текущей left cc right
	//ru/ и 4-х окружающих bottom
	
	double Q[] = { 0, 0, 0, 0 };
	double D[] = { 0, 0, 0, 0 };
	double P[] = { 0, 0, 0, 0 };
	double E[] = { 0, 0, 0, 0 };
	
	
	state old = pl->old;
	state old_prev = pl->old_prev;
	state *prev_row = pl->prev_row;
	
	// FIXME p == p??
	state cc	= update_state(old); // center
	// top FIXME error p[i][j+1] vs u[i+1][j]
	state top	= update_state(st[i+1][j]);
	state bottom	= update_state(prev_row[j-1]);
	state right	= update_state(st[i][j+1]);
	state left	= update_state(prev_old);

	state xc	= old; // ex-center
	state xtop	= st[i+1][j];
	state xbottom	= prev_row[j-1];
	state xright	= st[i][j+1];
	state xleft	= prev_old;

	
	// value of flow velocity in the central point
	double central_velocity = cc.calc_velocity();

	//--------------------------------------------------------------------------------
	D[0] = std::max(central_velocity, top.calc_velocity());	// top
	D[1] = std::max(central_velocity, bottom.calc_velocity());	// bottom
	D[2] = std::max(central_velocity, right.calc_velocity());	// right
	D[3] = std::max(central_velocity, left.calc_velocity());	// left
	
	P[0] = D[0] * (xtop.p - xc.p);
	P[1] = D[1] * (xc.p - xbottom.p);
	P[2] = D[2] * (xright.p - xc.p);
	P[3] = D[3] * (xc.p - xleft.p);
	
	E[0] = D[0] * (xtop.e - xc.e);
	E[1] = D[1] * (xc.e - xbottom.e);
	E[2] = D[2] * (xright.e - xc.e);
	E[3] = D[3] * (xc.e - xleft.e);
	
	double hx = st.hx;
	double hy = st.hy;
	double ht = st.ht;

	// mech = 0
	Q[0] = (xtop.u + xc.u)    - P[0];
	Q[1] = (xc.u + xbottom.u) - P[1];
	Q[2] = ((xright.v + xc.v) - P[2];
	Q[3] = ((xc.v + xleft.v)  - P[3];

	st[i][j].p = xc.p - 0.5 * ((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht;

	// mech = 1
	Q[0] = ((xtop.u * top.u + top.p) + (xc.u * cc.u + cc.p))          - P[0];
	Q[1] = ((cc.u * cc.u + cc.p) + (xbottom.u * bottom.u + bottom.p)) - P[1];
	Q[2] = ((xright.u * right.v + right.p) + (xc.u * cc.v + cc.p))    - P[2];
	Q[3] = ((xleft.u * left.v + left.p) + (xc.u * cc.v + cc.p))       - P[3];
	
	st[i][j].u = xc.u - 0.5 * ((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht;
	
	// mech = 2
	Q[0] = (xtop.v * top.v + top.p + xc.v * cc.v + cc.p)          - P[0];
	Q[1] = (xbottom.v * bottom.v + bottom.p + xc.v * cc.v + cc.p) - P[1];
	Q[2] = (xright.v * right.u + right.p + xc.u * cc.u + cc.p)    - P[2];
	Q[3] = (xleft.v * left.u + left.p + xc.u * cc.u + cc.p)       - P[3];

	st[i][j].v = xc.v - 0.5 * ((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht;
	
	// mech = 3
	Q[0] = (xtop.e * top.u + top.p * top.u + xleft.e * cc.u + cc.p * cc.u)             - E[0];
	Q[1] = (xbottom.e * bottom.u + bottom.p * bottom.u + xleft.e * cc.u + cc.p * cc.u) - E[1];
	Q[2] = (xright.e * right.v + right.p * right.p + xleft.e * cc.v + cc.p * cc.v)     - E[2];
	Q[3] = (xleft.e * left.v + left.p * left.v + xleft.e * cc.v + cc.p * cc.v)         - E[3];
	
	st[i][j].e = xc.e - 0.5 * ((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht;


}

void
execute(state *st, prev_layer *pl, double gamma)
{
	// velocity of the shock wave 
	
	pl->old = st[1][0];

	for (int k = 1; k < st->by + 1; ++k)
		pl->prev_row[k-1] = st[0][k];
	
	//--------------------------------------------------------------------------------
	printf("start execute\n");
	
	for (int i = 1; i < st->bx + 1; ++i) {
		for (int j = 1; j < st->by + 1; ++j) {
			pl->prev_old = pl->old;
			pl->old = st[i][j];

			// eval conservative variables
			evalVar(st, pl, i, j);
			pl->prev_row[j-1] = pl->prev_old;
		}
	}
	
	printf("end execute\n");
}



