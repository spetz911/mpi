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
	int i, j;
	
	st->hx = hx;
	st->hy = hy;
	st->ht = ht;
	
	for (i = 0; i < st->bx + 2; ++i)
		for (j = 0; j < st->by + 2; ++j) {
			st->rho[i][j]  = rho_g;
			st->rhoU[i][j] = rhoU_g;
			st->rhoV[i][j] = rhoV_g;
			st->rhoE[i][j] = rhoE_g;
		}
}

void execute(state *st, prev_layer *pl, double gamma)
{
	int i, j, k, meth;
	//en/ velocity of the shock wave 
	//ru/ скорость ударной волны
	double D[] = { 0, 0, 0, 0 };
	double Q[] = { 0, 0, 0, 0 };
	
	pl->rho_old  = st->rho[1][0];
	pl->rhoU_old = st->rhoV[1][0];
	pl->rhoV_old = st->rhoU[1][0];
	pl->rhoE_old = st->rhoE[1][0];

	for (k = 1; k < st->by + 1; ++k) {
		pl->rho_prev_row[k-1]  = st->rho[0][k];
		pl->rhoU_prev_row[k-1] = st->rhoU[0][k];
		pl->rhoV_prev_row[k-1] = st->rhoV[0][k];
		pl->rhoE_prev_row[k-1] = st->rhoE[0][k];
	}
		
	//--------------------------------------------------------------------------------
	printf("start execute\n"); 
	
	for (i = 1; i < st->bx + 1; ++i) {
		for (j = 1; j < st->by + 1; ++j) {
			//printf("(%d, %d) of [%d, %d] iter\n", i, j, bx, by);

			pl->rho_prev_old  = pl->rho_old;
			pl->rhoU_prev_old = pl->rhoU_old;
			pl->rhoV_prev_old = pl->rhoV_old;
			pl->rhoE_prev_old = pl->rhoE_old;
			
			pl->rho_old  = st->rho[i][j];
			pl->rhoU_old = st->rhoV[i][j];
			pl->rhoV_old = st->rhoU[i][j];
			pl->rhoE_old = st->rhoE[i][j];

			//en/ eval caonservative variables
			//ru/ рассчит. консервативные переменные
			for (meth = 0; meth < 4; ++meth)
				evalVar(st, pl, i, j, meth, gamma, Q, D);
			
			pl->rho_prev_row[j-1] = pl->rho_prev_old;
			pl->rhoU_prev_row[j-1] = pl->rhoU_prev_old;
			pl->rhoV_prev_row[j-1] = pl->rhoV_prev_old;
			pl->rhoE_prev_row[j-1] = pl->rhoE_prev_old;
		}
	}
	
	printf("end   execute\n");
}

void evalVar(state *st, prev_layer *pl, int i, int j, int meth, double gamma, double *Q, double *D)
{
	//--------------------------------------------------------------------------------
	evalQ(st, pl, i, j, meth, gamma, Q, D);

	switch(meth) {
	case 0:
		st->rho[i][j]  =-((Q[0] - Q[1]) / st->hx + (Q[2] - Q[3]) / st->hy) * st->ht + pl->rho_old;
		break;
	case 1:
		st->rhoU[i][j] =-((Q[0] - Q[1]) / st->hx + (Q[2] - Q[3]) / st->hy) * st->ht + pl->rhoU_old;
		break;
	case 2:
		st->rhoV[i][j] =-((Q[0] - Q[1]) / st->hx + (Q[2] - Q[3]) / st->hy) * st->ht + pl->rhoV_old;
		break;
	case 3:
		st->rhoE[i][j] =-((Q[0] - Q[1]) / st->hx + (Q[2] - Q[3]) / st->hy) * st->ht + pl->rhoE_old;
		break;
	default:
		break;
	}
}

void evalQ(state *st, prev_layer *pl, int i, int j, int meth, double gamma, double *Q, double *D)
{
	//en/ for each point it's coumpute separately
	//ru/ следующие параметры рассчитываются для каждой точки отдельно
	//                                               #
	//                                             $ * #
	//                                               $
	//                                               |
	//en/ get value for all variables                v
	//ru/ получим значения всех переменных в        top
	//ru/ 5-ти ближайших точках: текущей      left center right
	//ru/ и 4-х окружающих                         bottom
	
	//center
	double u_center   = pl->rhoU_old / pl->rho_old;
	double v_center   = pl->rhoV_old / pl->rho_old;
	double E_center   = pl->rhoE_old / pl->rho_old;
	double eps_center = (2 * E_center) / (u_center * u_center + v_center * v_center);
	double p_center   = (gamma - 1.0) * pl->rho_old * eps_center;
	
	//center
	double u_top      = st->rhoU[i+1][j] / st->rho[i+1][j];
	double v_top      = st->rhoV[i+1][j] / st->rho[i+1][j];
	double E_top      = st->rhoE[i+1][j] / st->rho[i+1][j];
	double eps_top    = (2 * E_top) / (u_top * u_top + v_top * v_top);
	double p_top      = (gamma - 1.0) * st->rho[i][j+1] * eps_top;
	
	//bottom
	double u_bottom   = pl->rhoU_prev_row[j-1] / pl->rho_prev_row[j-1];
	double v_bottom   = pl->rhoV_prev_row[j-1] / pl->rho_prev_row[j-1];
	double E_bottom   = pl->rhoE_prev_row[j-1] / pl->rho_prev_row[j-1];
	double eps_bottom = (2 * E_bottom) / (u_bottom * u_bottom + v_bottom * v_bottom);
	double p_bottom   = (gamma - 1.0) * pl->rho_prev_row[j-1] * eps_bottom;
	
	//right
	double u_right    = st->rhoU[i][j+1] / st->rho[i][j+1];
	double v_right    = st->rhoV[i][j+1] / st->rho[i][j+1];
	double E_right    = st->rhoE[i][j+1] / st->rho[i][j+1];
	double eps_right  = (2 * E_right) / (u_right * u_right + v_right * v_right);
	double p_right    = (gamma - 1.0) * st->rho[i][j+1] * eps_right;
	
	//left
	double u_left     = pl->rhoU_prev_old / pl->rho_prev_old;
	double v_left     = pl->rhoV_prev_old / pl->rho_prev_old;
	double E_left     = pl->rhoE_prev_old / pl->rho_prev_old;
	double eps_left   = (2 * E_left) / (u_left * u_left + v_left * v_left);
	double p_left     = (gamma - 1.0) * pl->rho_prev_old * eps_left;
	
	//en/ eval shock wave in the node
	//ru/ рассчит. скорость ударной волны в узле
	double mul     = gamma * (gamma - 1.0);
	
	//en/ value of flow velocity in the central point
	//ru/ значение скорости потока в центральной точке
	double central = sqrt(mul * eps_center) + sqrt(u_center * u_center + v_center * v_center);

	//--------------------------------------------------------------------------------
	D[0] = MAX(
		central,
		sqrt(mul * eps_top) + sqrt(u_top * u_top + v_top * v_top)
		);
	D[1] = MAX(
		central,
		sqrt(mul * eps_bottom) + sqrt(u_bottom  * u_bottom  + v_bottom  * v_bottom)
		);
	D[2] = MAX(
		central,
		sqrt(mul * eps_right) + sqrt(u_right * u_right + v_right * v_right)
		);
	D[3] = MAX(
		central,
		sqrt(mul * eps_left) + sqrt(u_left * u_left + v_left * v_left)
		);
	
	
	switch(meth) {
	case 0:
		Q[0] = (st->rhoU[i+1][j]       + pl->rhoU_old) / 2.0 - 
				D[0]*(st->rho[i+1][j] - pl->rho_old)           / 2.0;
		
		Q[1] = (pl->rhoU_prev_row[j-1] + pl->rhoU_old) / 2.0 - 
				D[1]*(pl->rho_old     - pl->rho_prev_row[j-1]) / 2.0;
		
		Q[2] = (st->rhoV[i][j+1]       + pl->rhoV_old) / 2.0 - 
				D[2]*(st->rho[i][j+1] - pl->rho_old)   		   / 2.0;
		
		Q[3] = (pl->rhoV_prev_old      + pl->rhoV_old) / 2.0 - 
				D[3]*(pl->rho_old    - pl->rho_prev_old) 	   / 2.0;
		break;
	case 1:
		Q[0] = (st->rhoU[i+1][j]       * u_top    + p_top    + pl->rhoU_old * u_center + p_center) / 2.0 -
				D[0]*(st->rho[i+1][j] - pl->rho_old)           / 2.0;
		
		Q[1] = (pl->rhoU_prev_row[j-1] * u_bottom + p_bottom + pl->rhoU_old * u_center + p_center) / 2.0 -
				D[1]*(pl->rho_old     - pl->rho_prev_row[j-1]) / 2.0;
		
		Q[2] = (st->rhoU[i][j+1]       * v_right  + p_right  + pl->rhoU_old * v_center + p_center) / 2.0 -
				D[2]*(st->rho[i][j+1] - pl->rho_old)           / 2.0;
		
		Q[3] = (pl->rhoU_prev_old      * v_left   + p_left   + pl->rhoU_old * v_center + p_center) / 2.0 -
				D[3]*(pl->rho_old - pl->rho_prev_old)          / 2.0;
		break;
	case 2:
		Q[0] = (st->rhoV[i+1][j]       * v_top    + p_top    + pl->rhoV_old * v_center + p_center) / 2.0 -
				D[0]*(st->rho[i+1][j] - pl->rho_old)           / 2.0;
		
		Q[1] = (pl->rhoV_prev_row[j-1] * v_bottom + p_bottom + pl->rhoV_old * v_center + p_center) / 2.0 -
				D[1]*(pl->rho_old     - pl->rho_prev_row[j-1]) / 2.0;
		
		Q[2] = (st->rhoV[i][j+1]       * u_right  + p_right  + pl->rhoU_old * u_center + p_center) / 2.0 -
				D[2]*(st->rho[i][j+1] - pl->rho_old)           / 2.0;
		
		Q[3] = (pl->rhoV_prev_old      * u_left   + p_left   + pl->rhoU_old * u_center + p_center) / 2.0 -
				D[3]*(pl->rho_old - pl->rho_prev_old)          / 2.0;
		break;
	case 3:
		Q[0] = (st->rhoE[i+1][j]       * u_top    + p_top    * u_top    + pl->rhoE_prev_old * u_center + p_center * u_center) / 2.0 -
				D[0]*(st->rhoE[i+1][j] - pl->rhoE_old)           / 2.0;
		
		Q[1] = (pl->rhoE_prev_row[j-1] * u_bottom + p_bottom * u_bottom + pl->rhoE_prev_old * u_center + p_center * u_center) / 2.0 -
				D[1]*(pl->rhoE_old     - pl->rhoE_prev_row[j-1]) / 2.0;
		
		Q[2] = (st->rhoE[i][j+1]       * v_right + p_right   * v_right  + pl->rhoE_prev_old * v_center + p_center * v_center) / 2.0 -
				D[2]*(st->rhoE[i][j+1] - pl->rhoE_old)           / 2.0;
		
		Q[3] = (pl->rhoE_prev_old      * v_left  + p_left    * v_left   + pl->rhoE_prev_old * v_center + p_center * v_center) / 2.0 -
				D[3]*(pl->rhoE_old - pl->rhoE_prev_old)          / 2.0;
		
		break;
	}
}

//void evalD(double gamma, double eps, double u, double v, double *D)
//{
//
//}
