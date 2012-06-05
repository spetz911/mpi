void addf(double *p, double *D, double *u)
{
  //*p = *p + 1.5;
  *u = *u + 0.5;
}

shock *create_shock(int nx, double gamma, double rho, double p, double u, double v, double rho_sh, double start_sh, double mul)
{
	double c02 = gamma * p * (1.0 / rho); //ru/ скорость звука покоя
	double D;
	double c12;
	
	shock *sh = (shock*) malloc(sizeof(shock));
	
	sh->shock_start = start_sh;
	sh->mul = mul;
	sh->startx = nx * start_sh;
	sh->top_osc = sh->bottom_osc = nx * mul;
	sh->tracef = &pline;
	
	sh->rho0 = rho;
	sh->p0   = p;
	sh->u0   = u;
	sh->eps0 = p /((gamma - 1.0) * rho);
	
	sh->rho1 = rho_sh;
	//sh->p1   = p * (((sh->u0 * sh->u0) / c02) * 2.0 * gamma - (gamma - 1.0)) / (gamma + 1.0);
	sh->p1 = p * (
		      ((gamma + 1.0) * (1.0 / rho) - (gamma - 1.0) * (1 / rho_sh)) 
		      /
		      ((gamma + 1.0) * (1.0 / rho_sh) - (gamma - 1.0) * (1 / rho)) 
        );
	double D1 = 1.5 * sqrt((gamma  * (sh->p0 / ( rho))));
	double D2 = 1.5 * sqrt(fabs(gamma  * (sh->p1 / ( rho_sh)))); 
	//printf("p1 = %f\n", sh->p1);
	c12 =  gamma * sh->p1 * (1 / sh->rho1);
	printf("c12 = %f\n", c12);
	sh->u1   = (c12 / (2 * gamma)) * ((gamma - 1.0) + (gamma + 1.0) * (sh->p0 / sh->p1)); addf(&sh->p1, &D2, &sh->u1);
	sh->eps1 = sh->p1 /((gamma - 1.0) * rho_sh);
	printf("before\n");
	printf("rho = %f, u = %f, p = %f, D = %f\n", sh->rho0, sh->u0, sh->p0, D1);//evalD(1.67, gamma, sh->p1 / ((gamma - 1) * sh->rho1)));
	printf("after\n");
	printf("rho = %f, u = %f, p = %f, D = %f\n", rho_sh, sh->u1, fabs(sh->p1), D2);
	
	
	return sh;
}


void shock_condition(state *st, shock *sh, double gamma, int partx, int party)
{
	int i, j, k;
	int f = true;
	int topx     = st->bx * (int) (process_id_g / party);
	int bottomx  = topx +  st->bx;
	int lefty    = st->by * (int) (process_id_g % party);
	int righty   = lefty + st->by;
	
	if (topx > sh->startx + sh->bottom_osc) {
	//ru/ мы находимся за фронтом ударной волны
	//ru/ ничего не делаем, т.к. начальные условия
	//ru/ не изменились
	}
	else {
	//ru/ рассчитываем ударные условия адиабаты Гюгонио
	
		int *cords = (int*) malloc(sizeof(int) * st->by);
		bool is_inner = false;
		
		for (i = 0; i < st->by; ++i) {
			cords[i] = sh->startx - (int) sh->mul * sh->tracef(st->hy * i);
			
			if (cords[i] >= topx && cords[i] <= bottomx)
				is_inner = true;
		}

	
		if (sh->startx - sh->top_osc > bottomx || is_inner) {
		//ru/ рассчит. ударные условия для всех точек
	
			for (i = 0; i < st->bx; ++i)
				for (j = 0; j < st->by; ++j) {
					st->rho[i][j]  = sh->rho1;
					st->rhoU[i][j] = sh->rho1 * sh->u1;
					st->rhoV[i][j] = 0.0;
					st->rhoE[i][j] = (sh->rho1 + sh->u1 * sh->u1 / 2.0) * sh->p1 / ((gamma - 1.0) * sh->rho1);
					f = false;
				}
		}
		else {
		  //printf("point 5, st->bx = %d, st->by = %d\n", st->bx, st->by);  
			for (i = 0; i < st->bx; ++i)
				for (j = 0; j < st->by; ++j) {
				        if (cords[i] > topx + j) {
					//ru/ находимся над волной
						st->rho[i][j]  = sh->rho1;
						st->rhoU[i][j] = sh->rho1 * sh->u1;
						st->rhoV[i][j] = 0.0;
						st->rhoE[i][j] = sh->rho1 * sh->p1 / ((gamma - 1.0) * sh->rho1);
					}
					else if (cords[i] == topx + j) {
					//ru/ находимя на волне - смешанные условия
						double S = 0.0, Stop = 0.0, Sbottom =0.0;
						
						S = st->bx * st->hx * st->by * st->hy;
						
						for (k = 0; k < st->by; ++k) {
							Sbottom = cords[k] * st->hx;
						}
						
						Stop = S - Sbottom;
						
						st->rho[i][j]  = (sh->rho0 * Sbottom + sh->rho1 * Stop) / (Sbottom + Stop);
						st->rhoU[i][j] = st->rho[i][j] * 
								(sh->u0 * Sbottom + sh->u1 * Stop) / (Sbottom + Stop);
						st->rhoV[i][j] = 0.0;
						st->rhoE[i][j] = st->rho[i][j] *
								(sh->u0 * Sbottom + 
								(sh->rho1 * sh->p1 / ((gamma - 1.0) * sh->rho1)) * Stop) 
								/ 
								(Sbottom + Stop);
					}
				}

		}
	}
}
