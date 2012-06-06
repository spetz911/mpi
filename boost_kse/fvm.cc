#include "fvm.h"

double MAX(double a, double b)
{
//	if (a < 0.0)
//		a = -a;
//	
//	if (b < 0.0)
//		b = -b;
	
	return a > b ?
			a : b;
}

void init_state(state *st, double hx, double hy, double ht)
{
	st->hx = hx;
	st->hy = hy;
	st->ht = ht;
	
	for (int i = 0; i < st->bx + 2; ++i)
		for (int j = 0; j < st->by + 2; ++j) {
			st->p[i][j] = p_g;
			st->u[i][j] = u_g;
			st->v[i][j] = v_g;
			st->e[i][j] = e_g;
		}
}

void
evalVar(state *st, prev_layer *pl, int i, int j, double gamma, double *Q, double *D);


void
execute(state *st, prev_layer *pl, double gamma)
{
	int i, j, k, meth;
	//en/ velocity of the shock wave 
	//ru/ скорость ударной волны
	
	pl->old = st->data[1][0];

	for (k = 1; k < st->by + 1; ++k) {
		pl->prev_row[k-1] = st->data[0][k];
	}
	
	//--------------------------------------------------------------------------------
	printf("start execute\n");
	
	for (i = 1; i < st->bx + 1; ++i) {
		for (j = 1; j < st->by + 1; ++j) {
			//printf("(%d, %d) of [%d, %d] iter\n", i, j, bx, by);

			pl->prev_old = pl->old;
			pl->old = st->data[i][j];

			// eval conservative variables
			evalVar(st, pl, i, j);
			pl->prev_row[j-1] = pl->prev_old;
		}
	}
	
	printf("end execute\n");
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
	D[0] = MAX(central_velocity, top.calc_velocity());
	D[1] = MAX(central_velocity, bottom.calc_velocity());
	D[2] = MAX(central_velocity, right.calc_velocity());
	D[3] = MAX(central_velocity, left.calc_velocity());
	
	double hx = st.hx;
	double hy = st.hy;
	double ht = st.ht;

	// mech = 0
	Q[0] = (xtop.u + xc.u) / 2.0 - 
	       D[0]*(xtop.p - xc.p) / 2.0;
	
	Q[1] = (xbottom.u + xc.u) / 2.0 - 
	       D[1]*(xc.p - xbottom.p) / 2.0;
	
	Q[2] = (xright.v + xc.v) / 2.0 - 
	       D[2]*(xright.p - xc.p) / 2.0;
	
	Q[3] = (xleft.v + xc.v) / 2.0 - 
	       D[3]*(xc.p - xleft.p) / 2.0;

	st[i][j].p = -((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht + xc.p;

	// mech = 1
	Q[0] = (xtop.u * top.u + top.p + xc.u * cc.u + cc.p) / 2.0 -
	       D[0]*(xtop.p - xc.p) / 2.0;
	
	Q[1] = (xbottom.u * bottom.u + bottom.p + cc.u * cc.u + cc.p) / 2.0 -
	       D[1]*(xc.p - xbottom.p) / 2.0;
	
	Q[2] = (xright.u * right.v + right.p + xc.u * cc.v + cc.p) / 2.0 -
	       D[2]*(xright.p - xc.p) / 2.0;
	
	Q[3] = (xleft.u * left.v + left.p + cc.u * cc.v + cc.p) / 2.0 -
	       D[3]*(xc.p - xleft.p) / 2.0;
	
	st[i][j].u = -((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht + xc.u;
	
	// mech = 2
	Q[0] = (xtop.v * top.v + top.p + xc.v * cc.v + cc.p) / 2.0 -
	       D[0]*(xtop.p - xc.p) / 2.0;
	
	Q[1] = (xbottom.v * bottom.v + bottom.p + xc.v * cc.v + cc.p) / 2.0 -
	       D[1]*(xc.p - xbottom.p) / 2.0;
	
	Q[2] = (xright.v * right.u + right.p + xc.u * cc.u + cc.p) / 2.0 -
	       D[2]*(xright.p - xc.p) / 2.0;
	
	Q[3] = (xleft.v * left.u + left.p + xc.u * cc.u + cc.p) / 2.0 -
	       D[3]*(xc.p - xleft.p) / 2.0;

	st[i][j].v = -((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht + xc.v;
	
	// mech = 3
	Q[0] = (xtop.e * top.u + top.p * top.u + xleft.e * cc.u + cc.p * cc.u) / 2.0 -
	       D[0]*(xtop.e - xc.e) / 2.0;
	
	Q[1] = (xbottom.e * bottom.u + bottom.p * bottom.u + xleft.e * cc.u + cc.p * cc.u) / 2.0 -
	       D[1]*(xc.e - xbottom.e) / 2.0;
	
	Q[2] = (xright.e * right.v + right.p * right.p + xleft.e * cc.v + cc.p * cc.v) / 2.0 -
	       D[2]*(xright.e - xc.e) / 2.0;
	
	Q[3] = (xleft.e * left.v + left.p * left.v + xleft.e * cc.v + cc.p * cc.v) / 2.0 -
	       D[3]*(xc.e - xleft.e) / 2.0;
	
	st[i][j].e = -((Q[0] - Q[1]) / hx + (Q[2] - Q[3]) / hy) * ht + xc.e;


}


