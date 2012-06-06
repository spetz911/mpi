

#include "shock.h"

double pline(double x)
{
	return 1.0;
}

static void
addf(double *p, double *D, double *u)
{
	//*p = *p + 1.5;
	*u = *u + 0.5;
}

Shock::Shock(int nx, double gamma_g, double rho, double p, double u, double v, double rho_sh, double start_sh, double mul)
{
	double c02 = gamma_g * p * (1.0 / rho); //ru/ скорость звука покоя
	double D;
	double c12;
	
	this->shock_start = start_sh;
	this->mul = mul;
	this->startx = nx * start_sh;
	this->top_osc = bottom_osc = nx * mul;
	this->tracef = &pline;
	
	this->rho0 = rho;
	this->p0   = p;
	this->u0   = u;
	this->eps0 = p / ((gamma_g - 1.0) * rho);
	
	this->rho1 = rho_sh;
	//p1   = p * (((u0 * u0) / c02) * 2.0 * gamma_g - (gamma_g - 1.0)) / (gamma_g + 1.0);
	this->p1 = p * (
		      ((gamma_g + 1.0) * (1.0 / rho) - (gamma_g - 1.0) * (1.0 / rho_sh)) 
		      /
		      ((gamma_g + 1.0) * (1.0 / rho_sh) - (gamma_g - 1.0) * (1.0 / rho)) 
        );
	double D1 = 1.5 * sqrt((gamma_g  * (p0 / ( rho))));
	double D2 = 1.5 * sqrt(fabs(gamma_g  * (p1 / ( rho_sh)))); 
	//printf("p1 = %f\n", p1);
	c12 =  gamma_g * this->p1 * (1 / rho1);
	printf("c12 = %f\n", c12);
	u1   = (c12 / (2 * gamma_g)) * ((gamma_g - 1.0) + (gamma_g + 1.0) * (p0 / p1));
	addf(&p1, &D2, &u1);
	eps1 = p1 /((gamma_g - 1.0) * rho_sh);
	printf("before\n");
	printf("rho = %f, u = %f, p = %f, D = %f\n", rho0, u0, p0, D1);//evalD(1.67, gamma_g, p1 / ((gamma_g - 1) * rho1)));
	printf("after\n");
	printf("rho = %f, u = %f, p = %f, D = %f\n", rho_sh, u1, fabs(p1), D2);

}


void
Shock::condition(Field &st, double gamma_g, int partx, int party)
{
	double bx = st.bx;
	double by = st.by;
	double hx = st.hx;
	double hy = st.hy;
	
	
	int f = true;
	int topx     = bx * (int) (process_id_g / party);
	int bottomx  = topx +  bx;
	int lefty    = by * (int) (process_id_g % party);
	int righty   = lefty + by;
	
	if (topx > startx + bottom_osc) {
	//ru/ мы находимся за фронтом ударной волны
	//ru/ ничего не делаем, т.к. начальные условия
	//ru/ не изменились
	} else {
	//ru/ рассчитываем ударные условия адиабаты Гюгонио
	
		int *cords = (int*) malloc(sizeof(int) * by);
		bool is_inner = false;
		
		for (int i = 0; i < by; ++i) {
			cords[i] = startx - (int) mul * tracef(hy * i);
			
			if (cords[i] >= topx && cords[i] <= bottomx)
				is_inner = true;
		}

	
		if (startx - top_osc > bottomx || is_inner) {
		//ru/ рассчит. ударные условия для всех точек
	
			for (int i = 0; i < bx; ++i)
				for (int j = 0; j < by; ++j) {
					st[i][j].p  = rho1;
					st[i][j].u = rho1 * u1;
					st[i][j].v = 0.0;
					st[i][j].e = (rho1 + u1 * u1 / 2.0) * p1 / ((gamma_g - 1.0) * rho1);
					f = false;
				}
		}
		else {
		  //printf("point 5, bx = %d, by = %d\n", bx, by);  
			for (int i = 0; i < bx; ++i)
				for (int j = 0; j < by; ++j) {
				        if (cords[i] > topx + j) {
					//ru/ находимся над волной
						st[i][j].p  = rho1;
						st[i][j].u = rho1 * u1;
						st[i][j].v = 0.0;
						st[i][j].e = rho1 * p1 / ((gamma_g - 1.0) * rho1);
					}
					else if (cords[i] == topx + j) {
					//ru/ находимя на волне - смешанные условия
						double S = 0.0, Stop = 0.0, Sbottom =0.0;
						
						S = bx * hx * by * hy;
						
						for (int k = 0; k < by; ++k) {
							Sbottom = cords[k] * hx;
						}
						
						Stop = S - Sbottom;
						
						st[i][j].p  = (rho0 * Sbottom + rho1 * Stop) / (Sbottom + Stop);
						st[i][j].u = st[i][j].p * 
								(u0 * Sbottom + u1 * Stop) / (Sbottom + Stop);
						st[i][j].v = 0.0;
						st[i][j].e = st[i][j].p *
								(u0 * Sbottom + 
								(rho1 * p1 / ((gamma_g - 1.0) * rho1)) * Stop) 
								/ 
								(Sbottom + Stop);
					}
				}

		}
	}
}


