#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/timer.hpp>

#include <iostream>

#include "kse.h"
#include "fvm.h"

std::string	processor_name_g;

int	process_id_g;
int	process_num_g;

double rho_g;
double rhoU_g;
double rhoV_g;
double rhoE_g;

// FIXME
int nx = 100;
int ny = 100;
int nt = 10;

//number of points for each process
int bx;
int by;

double hx = 0.0001;
double hy = 0.0001;
double ht = 0.0000001;

double gamma_g = 1.66666666667;

#include "stuff.cc"


namespace mpi = boost::mpi;

class KSE {
	static const int nx = 100;
	static const int ny = 100;
	static const int nt = 10;

};

double *result;              //finally result
bound *b;                    //four bounds
state *st;                   //values
prev_layer *pl;              //stored values
shock *sh;



void
log_(const std::string &str)
{
	fprintf(stdout,"Process %d of %d on %s: %s\n",
		process_id_g, process_num_g, processor_name_g.c_str(), str.c_str());
}

void
init_globals()
{
	double rho	= 1.0;
	double p	= 1.0;
	double u	= 1.0;
	double v	= 0.0;
	
	// discontinuity surface
	double rho_sh	= 0.5;
	double sh_start	= 0.333;
	double mul	= 0.2;
	
	// Curant-Friedrichs-Levi condition
	double c_speed = sqrt((gamma_g * p) / rho);

	rho_g = rho;
	rhoU_g = rho * u;
	rhoV_g = rho * v;
	rhoE_g = rho * (p / ((gamma_g - 1.0) * rho)) * (u * u + v * v) / 2.0;

	while (ht > 1.0 / ((v + c_speed) * (sqrt(hx * hx + hy * hy) / (hx * hy)))) {
		ht /= 10.0;
	}

	hx = ((1.0 / nx) < hx) ?
			1.0 / nx : hx;
	hy = ((1.0 / ny) < hy) ?
			1.0 / ny : hy;
	ht = hx * hy;


	//init structures
	st = create_state(bx, by);
	init_state(st, hx, hy, ht);
	b  = create_bound(bx, by);
	pl = create_prev_layer(st->by);
	sh = create_shock(nx, gamma_g, rho, p, u, v, rho_sh, sh_start, mul);
	

}


int
main(int argc, char* argv[])
{
	//partition along X and Y axis
	int partx = 1;
	int party = 1;
	
	

	int i, j, k, l, stt;          //stat;

	
	

	//------------------------------------------------------------


	mpi::environment env(argc, argv);
	mpi::communicator world;
	mpi::timer timer;
	process_num_g = world.size();
	set_partition(&partx, &party, nx, ny, process_num_g);
	process_id_g = world.rank();
	processor_name_g = env.processor_name();
	
	set_step_partition(partx, party, &bx, &by, nx, ny, process_id_g);

	if (process_id_g == 0)
		printf("partx = %d party = %d", partx, party);
	
	init_globals();

	printf("bx = %d by = %d process_id_g = %d", nx, ny, bx, by, process_id_g);
	
	
	//------------------------------------------------------------
	//while ((stat = MPI_Barrier(MPI_COMM_WORLD)) == 0) ;
	world.barrier();
	printf("initialize complete");
	
	//------------------------------------------------------------
	timer.restart();

 	printf("nt = %d", nt);
	
	//Huhonio
	shock_condition(st, sh, gamma_g, partx, party);

	//***************************************
//	MPI_Barrier(MPI_COMM_WORLD);
//	exit(0);
		
	/* for (stt = 0; stt < nt; ++stt) { */
	/* 	prmsg("exchange"); */
	/* 	echange_bound(st, b, partx, party); */
	/* 	prmsg("evaluate"); */
	/* 	execute(st, pl, gamma_g); */
	/* } */

	//eval_block(block_val, bx, by, partx, party);

	printf("stage 5");

	// Wait for all processes within a communicator to reach the barrier.
	// This routine is a collective operation that blocks each process until all processes
	// have entered it, then releases all of the processes "simultaneously".
	world.barrier();
	
	char str[100];
	
	sprintf(str, "out/rho_%d",  process_id_g);
	show_result(bx, by, st->rho,  str, "Rho");
	sprintf(str, "out/rhoU_%d", process_id_g);
	show_result(bx, by, st->rhoU, str, "RhoU");
	sprintf(str, "out/rhoV_%d", process_id_g);
	show_result(bx, by, st->rhoV, str, "RhoV");
	sprintf(str, "out/rhoE_%d", process_id_g);
	show_result(bx, by, st->rhoE, str, "RhoE");

	
	printf("stage 6");
	//***************************************
	world.barrier();

//	exit(0);
	clear_state(st);
	clear_bound(b);
	clear_prev_layer(pl);

	printf("stage 7");
	//stop mpi

	if (process_id_g == 0)
		printf(str, "wall clock time = %.10f", timer.elapsed());

	//	MPI_Finalize();
	// decstructor ~env();

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


