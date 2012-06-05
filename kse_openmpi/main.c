#include "mpi.h"
#include "fvm.h"
#include "kse.h"

int        gl_rank;
int        gl_num_procs;
char      *gl_processor_name;
MPI_Status gl_status;
double gl_rho;
double gl_rhoU;
double gl_rhoV;
double gl_rhoE;
double M = 1.0;

void addf(double *p, double *D, double *u)
{
  //*p = *p + 1.5;
  *u = *u + 0.5;
}
void prmsg(const char *str)
{
	fprintf(stdout,"Process %d of %d on %s: %s\n",
		  gl_rank, gl_num_procs, gl_processor_name, str);
}

void echange_bound(state *st, bound *b, int partx, int party)
{
	char str[100];
	int i, l;
	double **block_val;

	for (l = 0; l < 4; ++l) {

		switch (l) {
		case 0:
			block_val = st->rho;
			break;
		case 1:
			block_val = st->rhoU;
			break;
		case 2:
			block_val = st->rhoV;
			break;
		case 3:
			block_val = st->rhoE;
			break;
		}

		prmsg("control point");
		init_bound(b, block_val);
	 	prmsg("copy bounds complete");

		//for top    --> bottom layer
		if (gl_rank >= party) {
//			sprintf(str, "send to %d start", gl_rank - party);
//			prmsg(str);
						
			for (i = 1; i < b->by + 1; ++i)
				b->tmpY[i-1] = block_val[1][i];

			MPI_Send(b->tmpY, b->by, MPI_DOUBLE, gl_rank - party, gl_rank - party,
					MPI_COMM_WORLD);
//			sprintf(str, "send to %d complete", gl_rank - party);
//			prmsg(str);
		}

		if (gl_rank < party * (partx - 1))
			MPI_Recv(b->bottom_bound, b->by, MPI_DOUBLE, gl_rank + party,
					MPI_ANY_TAG, MPI_COMM_WORLD, &gl_status);
		else
			for (i = 0; i < b->by; ++i)
				b->bottom_bound[i] = top_boundary_condition(l, i * st->hy, b->by);

		prmsg("stage 1");
//		MPI_Barrier( MPI_COMM_WORLD);
		
		//for bottom --> top     layer
		if (gl_rank < party * (partx - 1)) {
			for (i = 1; i < b->by + 1; ++i)
				b->tmpY[i-1] = block_val[b->bx][i];

			MPI_Send(b->tmpY, b->by, MPI_DOUBLE, gl_rank + party, gl_rank + party,
					MPI_COMM_WORLD);
		}

		if (gl_rank >= party)
			MPI_Recv(b->top_bound, b->by, MPI_DOUBLE, gl_rank - party,
					MPI_ANY_TAG, MPI_COMM_WORLD, &gl_status);
		else
			for (i = 0; i < b->by; ++i)
				b->top_bound[i] = bottom_boundary_condition(l, i * st->hy, b->by);

		prmsg("stage 2");

		//for left   --> right   layer
		if ((gl_rank + 1) % party != 0) {
			for (i = 1; i < b->bx + 1; ++i)
				b->tmpX[i-1] = block_val[i][b->by];

			MPI_Send(b->tmpX, b->bx, MPI_DOUBLE, gl_rank + 1, gl_rank + 1,
					MPI_COMM_WORLD);
		}

		if (gl_rank % party > 0)
			MPI_Recv(b->left_bound, b->bx, MPI_DOUBLE, gl_rank - 1, MPI_ANY_TAG,
					MPI_COMM_WORLD, &gl_status);
		else
			for (i = 0; i < b->bx; ++i)
				b->left_bound[i] = left_boundary_condition(l, i * st->hx, b->bx);

		prmsg("stage 3");

		//for right  --> left    layer
		if (gl_rank % party > 0) {
			for (i = 1; i < b->bx + 1; ++i)
				b->tmpX[i-1] = block_val[i][1];

			MPI_Send(b->tmpX, b->bx, MPI_DOUBLE, gl_rank - 1, gl_rank - 1,
					MPI_COMM_WORLD);
		}

		if ((gl_rank + 1) % party != 0)
			MPI_Recv(b->right_bound, b->bx, MPI_DOUBLE, gl_rank + 1, MPI_ANY_TAG,
					MPI_COMM_WORLD, &gl_status);
		else
			for (i = 0; i < b->bx; ++i)
				b->right_bound[i] = right_boundary_condition(l, i * st->hx, b->bx);

		copy_bound(b, block_val);
		MPI_Barrier( MPI_COMM_WORLD);

		prmsg("stage 4");
	}	
}

double evalD(double M, double gamma, double eps)
{
  return M * sqrt(gamma * (gamma - 1.0) * eps);
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
	int topx     = st->bx * (int) (gl_rank / party);
	int bottomx  = topx +  st->bx;
	int lefty    = st->by * (int) (gl_rank % party);
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

int main(int argc, char *argv[])
{
	int rank;                     //process id
	int num_procs;                //number of process     
	double startwtime;            //start time
	double endwtime;              //end   time
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int  namelen;
	MPI_Comm comm_gr;
	char str[100];

	int nx, ny, nt;
	int partx = 1, party = 1;     //partition along X and Y axis
	int bx, by;                   //number of points for each process
	int i, j, k, l, stt;          //stat;
	double hx, hy, ht;
	double gamma;
	double *result;              //finally result
	bound *b;                    //four bounds
	state *st;                   //values
	prev_layer *pl;              //stored values
	shock *sh;
	
	double c;                      //ru/ скорость звука
	double rho, p, u, v;
	double rho_sh, sh_start, mul;  //ru/поверхность разрыва
		
	//-------------


	FILE *F = fopen("conf/init-config", "r");
	
	fscanf(F, "%d ",  &nx);
	fscanf(F, "%d ",  &ny);
	fscanf(F, "%d ",  &nt);
	fscanf(F, "%lf ", &hx);
	fscanf(F, "%lf ", &hy);
	fscanf(F, "%lf ", &ht);
	fscanf(F, "%lf ", &gamma);
	fscanf(F, "%lf ", &rho);
	fscanf(F, "%lf ", &p);
	fscanf(F, "%lf ", &u);
	fscanf(F, "%lf ", &v);
	fscanf(F, "%lf ", &rho_sh);
	fscanf(F, "%lf ", &sh_start);
	fscanf(F, "%lf ", &mul);
	
	fclose(F);
	
	//en/ Curant-Friedrichs-Levi condition
	c = sqrt((gamma * p) / rho);
	
	while (ht > 1.0 / ((v + c) * (sqrt(hx * hx + hy * hy) / (hx * hy)))) {
		ht /= 10.0;
	}


	hx = ((1.0 / nx) < hx) ?
			1.0 / nx : hx;
	hy = ((1.0 / ny) < hy) ?
			1.0 / ny : hy;
	ht = hx * hy;
	//------------------------------------------------------------
	gl_rho = rho;
	gl_rhoU = rho * u;
	gl_rhoV = rho * v;
	gl_rhoE = rho * (p / ((gamma - 1.0) * rho)) * (u * u + v * v) / 2.0;
	
	//------------------------------------------------------------
	srand(time(NULL));
	printf("nx = %d ny = %d\nhx = %f hy = %f ht = %f\n", nx, ny, hx, hy, ht);
	
	//init mpi
	//------------------------------------------------------------
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	set_partition(&partx, &party, nx, ny, num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

	gl_num_procs      = num_procs;
	gl_rank           = rank;
	gl_processor_name = processor_name;

	prmsg("start");

	//set steps partition

	set_step_partition(partx, party, &bx, &by, nx, ny, rank);

	if (rank == 0) {
		sprintf(str, "partx = %d party = %d", partx, party);
		prmsg(str);
	}

	sprintf(str, "bx = %d by = %d rank = %d", nx, ny, bx, by, rank);
	prmsg(str);
	
	//init structures
	st = create_state(bx, by);
	init_state(st, hx, hy, ht);
	b  = create_bound(bx, by);
	pl = create_prev_layer(st->by);
	sh = create_shock(nx, gamma, rho, p, u, v, rho_sh, sh_start, mul);
		
	//------------------------------------------------------------
	//while ((stat = MPI_Barrier(MPI_COMM_WORLD)) == 0) ;
	sprintf(str, "barrier val: %d", MPI_Barrier(MPI_COMM_WORLD));
	prmsg(str);
	prmsg("initialize complete");
    

	//------------------------------------------------------------
	startwtime = MPI_Wtime();

 	sprintf(str, "nt = %d", nt);
	prmsg(str);
 	
	//Huhonio
	shock_condition(st, sh, gamma, partx, party);
	//***************************************
//	MPI_Barrier(MPI_COMM_WORLD);
//	exit(0);
		
	/* for (stt = 0; stt < nt; ++stt) { */
	/* 	prmsg("exchange"); */
	/* 	echange_bound(st, b, partx, party); */
	/* 	prmsg("evaluate"); */
	/* 	execute(st, pl, gamma); */
	/* } */

	//eval_block(block_val, bx, by, partx, party);

	prmsg("stage 5");
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(str, "out/rho_%d",  rank);
	print(bx, by, st->rho,  str, "Rho");
	sprintf(str, "out/rhoU_%d", rank);
	print(bx, by, st->rhoU, str, "RhoU");
	sprintf(str, "out/rhoV_%d", rank);
	print(bx, by, st->rhoV, str, "RhoV");
	sprintf(str, "out/rhoE_%d", rank);
	print(bx, by, st->rhoE, str, "RhoE");

	
	prmsg("stage 6");
	//***************************************
	MPI_Barrier(MPI_COMM_WORLD);
//	exit(0);
	clear_state(st);
	clear_bound(b);
	clear_prev_layer(pl);

	prmsg("stage 7");
	//stop mpi

	if (rank == 0) {
		//end timer
		endwtime = MPI_Wtime();
		sprintf(str, "wall clock time = %.10f",
				endwtime-startwtime);
		prmsg(str);
	}


	MPI_Finalize();

	return 0;
}


/*
  if(rank >= 1) {
	//MPI_Send(&blocks[myid],sizeof(blockPartition),MPI_CHAR,0,0,MPI_COMM_WORLD);
	//MPI_Send(b->top_bound, partx,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	MPI_Send(block_val, bx * by, MPI_INT, 0, 0, MPI_COMM_WORLD);
	//MPI_Send(b->b->bottom_bound, partx,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	//MPI_Send(b->left_bound, partx,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	//MPI_Send(b->right_bound, partx,MPI_DOUBLE,0,0,MPI_COMM_WORLD);	
    }



    if(rank == 0) {
      u = (double *) malloc(sizeof(double) * bx * by * num_procs);
      tmp = (double *) malloc(sizeof(double) * bx * by);

      j = bx * by * num_procs;

      for (i = 0; i < j; ++i)
    	  u[i] = 0.0;
      
      j = bx * by;

      for (i = 0; i < j; ++i)
    	  tmp[i] = 0.0;
      
      for (k = 0; k < j; ++k)
    	  u[k] = block_val[i];

      for(i = 1; i < num_procs; ++i) {
    	  MPI_Recv(tmp, bx * by, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &gl_status);
    	  j = bx * by;

    	  for (k = 0; k < j; ++k)
    		  u[k + i * j] = tmp[k];
		  }

      print(u, partx, party, bx, by);
    }


    
      if (rank == 0) {
    	  prmsg("gather start");
      result = (double *) malloc(sizeof(double) * bx * by * num_procs);
      j = bx * by * num_procs;
      
      for (i = 0; i < j; ++i)
    	  result[i] = 0.0;
      
      MPI_Gather(block_val, bx * by, MPI_DOUBLE, result, bx * by, MPI_DOUBLE, rank, MPI_COMM_WORLD);
      
      print(result, partx, party, bx, by);
      
      prmsg("gather complete");
    }
*/
