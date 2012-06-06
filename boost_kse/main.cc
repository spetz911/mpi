#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/timer.hpp>

#include <iostream>

#include "kse.h"
#include "fvm.h"
#include "field.h"
#include "shock.h"

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


namespace mpi = boost::mpi;

class KSE {
	static const int nx = 100;
	static const int ny = 100;
	static const int nt = 10;

};




void
log_(const std::string &str)
{
	fprintf(stdout,"Process %d of %d on %s: %s\n",
		process_id_g, process_num_g, processor_name_g.c_str(), str.c_str());
}



int
main(int argc, char* argv[])
{
	//partition along X and Y axis
	int partx = 1;
	int party = 1;
	
	double *result;              //finally result
	
	//------------------------------------------------------------


	mpi::environment env(argc, argv);
	mpi::communicator world;
	mpi::timer timer;
	process_num_g = world.size();

	set_partition(&partx, &party, nx, ny);

	process_id_g = world.rank();
	processor_name_g = env.processor_name();
	
	set_step_partition(partx, party, &bx, &by, nx, ny);

	if (process_id_g == 0)
		printf("partx = %d party = %d", partx, party);
	
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


	printf("bx = %d by = %d process_id_g = %d", nx, ny, bx, by, process_id_g);

	//init structures
	Field st(bx, by);
	st.init(hx, hy, ht);
	bounds b(bx, by);
	
	prev_layer pl(by);
	Shock sh(nx, gamma_g, rho, p, u, v, rho_sh, sh_start, mul);
		
	
	//------------------------------------------------------------
	//while ((stat = MPI_Barrier(MPI_COMM_WORLD)) == 0) ;
	world.barrier();
	printf("initialize complete");
	
	//------------------------------------------------------------
	timer.restart();

 	printf("nt = %d", nt);
	
	//Huhonio
	sh.condition(st, gamma_g, partx, party);

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
	
	st.show_result("out/");

	
	printf("stage 6");
	//***************************************
	world.barrier();

//	exit(0);

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


