//#include <iostream>
//#include <omp.h>

//using namespace std;

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structures.h"
#include "work_withCell.h"
#include "functions.h"
#include "const.h"


/*
double fmax(double a, double b)
{
	if (a > b)
		return a;
	else return b;
}*/


int main(int argc,char *argv[])
{
	int np;				// number of processes (summary x, y)	//by input as a parameter in exe
	//int nX = atoi(argv[1]), nY = atoi(argv[2]);			//number X, Y points					// by input in program
	int nX = 10;
	int nY = 100;
	int wave_index = 3;
	//int bsX, bsY;		//blocksize of X, Y
	int npX, npY;		//number processes X, Y
	int idx;
	
	double dx = h / nX;
	double dy = l / nY;
	double startwtime, endwtime;
	startwtime = endwtime = 0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &idx);

	npY = separate(np, nX, nY); // separate our grid
	npX = (np - 1) / npY;	

	int master_idx = 0;

	if (idx == master_idx)	// master
	{
		double dx = l / nX;
		double dy = h / nY;
		//printf("MASTER %d \n\n", idx);

		startwtime = MPI_Wtime();

		int i, j;
		int bsX;
		int bsY;

		double speedU, speedV;
		double a, b, c;
		double D;
		TCell cur;

		speedU = 0.0;
		speedV = 0.0;
		double eps3 = (pressure3 / (gamma1 - 1.0) / ro3);

		D = -M * sqrt(gamma1* (gamma1 - 1.0) * eps3);
		a = ro3 * (D - speedV);
		b = ro3 * pow((D - speedV), 2.0) + pressure3;
		c = eps3 + (pressure3 / ro3 ) +  pow((D - speedV), 2.0) / 2.0;

		cur.speedU = speedU;
		cur.ro = ( b * gamma1 + sqrt(b*b * gamma1*gamma1 - 2.0* a * a * c * (gamma1 * gamma1 - 1.0)) ) / ( 2.0 * c *( gamma1 - 1.0));
		cur.speedV = D - a / cur.ro;
		cur.pressure = b - a*a / cur.ro;
		cur.eps = cur.pressure / ( (gamma1 - 1.0 )* cur.ro );
		
		double curx = 0.0;
		double cury = h;
		double prevy = h;
		TMessage message;
		message.gcell = cur;

		printf("npX = %d, npY = %d\n\n", npX, npY);
		for (i = 0; i < npX; ++i)
		{
			for (j = 0; j < npY; ++j)
			{
				int id_process = j * npX + i + 1;
				bsX = nX / npX;		
				bsY = nY / npY;

				if (id_process % npX == 0)
				{
					bsX += (nX % npX);
				}

				if (id_process >=  np - npX)
				{
					bsY += (nY % npY);
				}

				cury -= bsY * dy;

				message.start_x = curx;
				message.start_y = cury;

				//printf("! idx = %d, bsX = %d, bsY = %d\n", id_process, bsX, bsY);
				//printf("st_x = %f, st_y = %f\n", message.start_x * 1000000, message.start_y * 1000000);
				//printf("process id = %d\n", j * npX + i + 1);

				MPI_Send(&message, sizeof(struct TMessage), MPI_CHAR, id_process, 2, MPI_COMM_WORLD );
			}

			curx += bsX * dx;
			cury = prevy;
		}
	}

	if (idx != master_idx)		// usual process
	{
		int work_id = 1;

		double dx = l / nX;
		double dy = h / nY;

		Cell cell = cellInit(cell, nX, nY, idx, npX, npY, np);

		TMessage mes1;
		// Receive from master
		MPI_Status status0;
		MPI_Recv(&mes1, sizeof(struct TMessage), MPI_CHAR, master_idx, 2, MPI_COMM_WORLD, &status0);
		
		cell.start_x = mes1.start_x;
		cell.start_y = mes1.start_y;
		
		TCell* layer1 = (TCell*) malloc( cell.bsX * cell.bsY * sizeof(TCell));
		TCell* layer2 = (TCell*) malloc( cell.bsX * cell.bsY * sizeof(TCell));
		//TCell* layer3 = (TCell*) malloc( cell.bsX * cell.bsY * sizeof(TCell));

		TCell* prev_layer = layer1;
		TCell* cur_layer = layer2;
		TCell* temp_layer;		
		

		layerInit(prev_layer, mes1.gcell, cell, dx, dy, idx);

		TCell* buf = (TCell*) malloc( (cell.bsX-2) * sizeof(TCell));

		for (int t = 0; t < K; ++t)
		{
			sendReceive(prev_layer, cell, npX);
			solve(prev_layer, cur_layer, cell, dx, dy, 0);		// u^{*} = dt * L(u^{k})
	
			updateLayer(cur_layer,cell);
			setBorder(cur_layer, prev_layer, cell);
			
			temp_layer = prev_layer;
			prev_layer = cur_layer;
			cur_layer = temp_layer;



			if ((t % printstep) == 0)
			{
				for (int i = 1; i < cell.bsY - 1; ++i)
				{
					for (int j = 1; j < cell.bsX - 1; ++j)
					{
						buf[j - 1] = prev_layer[i*cell.bsX + j];
					}
					MPI_Send(buf, (cell.bsX - 2) * sizeof(struct TCell), MPI_CHAR, master_idx, 2, MPI_COMM_WORLD );
				}
			}
		}

		//TCell* buf = (TCell*) malloc( (cell.bsX-2) * sizeof(TCell));
	
		for (int i = 1; i < cell.bsY - 1; ++i)
		{
			for (int j = 1; j < cell.bsX - 1; ++j)
				{
					buf[j - 1] = prev_layer[i*cell.bsX + j];
				}
				MPI_Send(buf, (cell.bsX - 2) * sizeof(struct TCell), MPI_CHAR, master_idx, 2, MPI_COMM_WORLD );
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (idx == master_idx)	// master
	{
		MPI_Status status0;
		int i, j;
		int max_x;

		int bsX;
		int bsY;

		FILE* fout_u;
		FILE* fout_v;
		FILE* fout_ro;
		FILE* fout_p;
		FILE* fout_eps;
		fout_u = fopen("output_speedU.txt", "w");
		fout_v = fopen("output_speedV.txt", "w");
		fout_ro = fopen("output_ro.txt", "w");
		fout_p = fopen("output_pressure.txt", "w");
		fout_eps = fopen("output_eps.txt", "w");

		fprintf(fout_u, "Manipulate[ListPlot3D[{\n");
		fprintf(fout_v, "Manipulate[ListPlot3D[{\n");
		fprintf(fout_ro, "Manipulate[ListPlot3D[{\n");
		fprintf(fout_eps, "Manipulate[ListPlot3D[{\n");
		fprintf(fout_p, "Manipulate[ListPlot3D[{\n");

		for (int iprint = 0; iprint <= (K / printstep); ++iprint)
		{

		//fprintf(fout_u, "\n%d\n\n", iprint * printstep);
		//fprintf(fout_v, "\n%d\n\n", iprint * printstep);
		//fprintf(fout_ro, "\n%d\n\n", iprint * printstep);
		//fprintf(fout_eps, "\n%d\n\n", iprint * printstep);
		//fprintf(fout_p, "\n%d\n\n", iprint * printstep);

		fprintf(fout_u, "{ ");
		fprintf(fout_v, "{ ");
		fprintf(fout_ro, "{ ");
		fprintf(fout_eps, "{ ");
		fprintf(fout_p, "{ ");

		max_x = nX / npX;
		max_x += (nX % npX);

		bsY = nY / npY;

		if ((master_idx + 1) >=  np - npX)
		{
			bsY += (nY % npY);
		}

		TCell* buf = (TCell*) malloc (max_x * sizeof(TCell));
		double curx = 0.0;
		double cury = h;
		double prevy = h;
		dx = l / nX;
		dy = h / nY;

		for (i = 0; i < npY; ++i)
		{
			for (int vlines = 0; vlines < bsY; vlines++)
			{
				// send all horizontal lines
				for (j = 0; j < npX; ++j)
				{
					int id_process = i * npX + j + 1;
					bsX = nX / npX;
					bsY = nY / npY;

					if (id_process % npX == 0)
					{
						bsX += (nX % npX);
					}

					if (id_process >=  np - npX)
					{
						bsY += (nY % npY);
					}

					MPI_Recv(buf, bsX * sizeof(struct TCell), MPI_CHAR, id_process, 2, MPI_COMM_WORLD, &status0);

					for (int x_size = 0; x_size < bsX; ++x_size)
					{
						if ( (vlines % 5 == 0) || (vlines == bsY - 1))

						if ((x_size == bsX - 1) && (vlines == bsY - 1) && (id_process == np - 1))
						{
							double tmpx = (curx + x_size * dx)*1000000;
							double tmpy = (cury - (vlines+1) * dy)*1000000;
							
							fprintf(fout_u,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].speedU);
							fprintf(fout_v,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].speedV);
							fprintf(fout_ro,  "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].ro);
							fprintf(fout_eps, "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].eps);
							fprintf(fout_p,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].pressure);
						}
						else
						{
							double tmpx = (curx + x_size * dx)*1000000;
							double tmpy = (cury - (vlines+1) * dy)*1000000;
							
							fprintf(fout_u,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].speedU);
							fprintf(fout_v,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].speedV);
							fprintf(fout_ro,  "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].ro);
							fprintf(fout_eps, "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].eps);
							fprintf(fout_p,   "{%f, %f, %f} ", tmpx, tmpy, buf[x_size].pressure);
					
						}
					}
					curx = 0.0;
				}
				fprintf(fout_u, "\n");
				fprintf(fout_v, "\n");
				fprintf(fout_ro, "\n");
				fprintf(fout_eps, "\n");
				fprintf(fout_p, "\n");
				// new line
			}
			curx = 0.0;
			cury -= bsY * dy;
			// new process line
		}

		fprintf(fout_u, " },\n\n");
		fprintf(fout_v, " },\n\n");
		fprintf(fout_ro, " },\n\n");
		fprintf(fout_eps, " },\n\n");
		fprintf(fout_p, " },\n\n");

		}
		const char str[] = "}[[t]], Mesh -> None, PlotRange -> Full, ImageSize -> 600, ColorFunction -> \"TemperatureMap\"], {t, 1, %d, 1}]\n";
		fprintf(fout_u, str, (K / printstep) + 1);
		fprintf(fout_v, str, (K / printstep) + 1);
		fprintf(fout_ro, str, (K / printstep) + 1);
		fprintf(fout_eps, str, (K / printstep) + 1);
		fprintf(fout_p, str, (K / printstep) + 1);

		fclose(fout_u);
		fclose(fout_v);
		fclose(fout_ro);
		fclose(fout_eps);
		fclose(fout_p);
		endwtime = MPI_Wtime();
		printf("time = %lf\n", endwtime - startwtime);
	}

	MPI_Finalize();

	return 0;
}
