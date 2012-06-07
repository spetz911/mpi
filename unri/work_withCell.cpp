#include <stdlib.h>
#include <stdio.h>
#include "work_withCell.h"
#include "const.h"


void layerInit(TCell *cur, TCell gugonio, Cell cell, double dx, double dy, int id_x)
{
	double x , y;
	for (int i = 1 ; i < cell.bsY-1; i++)
	{
		for (int j = 1 ; j < cell.bsX - 1; j++)
		{
			x = cell.start_x + (j-1) * dx;
			y = cell.start_y + (i - 1) * dy;
			cur[(cell.bsY - 1 - i) * cell.bsX + j] = Psi(cur[i * cell.bsX + j], gugonio, x, y);
		}

	}
	setBorder(cur, cur, cell);
}




double getH2(double x, double y)
{
	return (h2 + a0* cos(PI * x / l ) );
}

void make_layer(TCell* cur_layer, TCell* layer1, Cell cell)
{
	int i, j;
	for (j = 0; j < cell.bsX; ++j)
	{
		for (i = 0; i < cell.bsY; ++i)
		{
			cur_layer[i*cell.bsX + j].eps = 0.5 * cur_layer[i*cell.bsX + j].eps + 0.5 * layer1[i*cell.bsX + j].eps;
			cur_layer[i*cell.bsX + j].speedU = 0.5 * cur_layer[i*cell.bsX + j].speedU + 0.5 * layer1[i*cell.bsX + j].speedU;
			cur_layer[i*cell.bsX + j].speedV = 0.5 * cur_layer[i*cell.bsX + j].speedV + 0.5 * layer1[i*cell.bsX + j].speedV;
			cur_layer[i*cell.bsX + j].value1 = 0.5 * cur_layer[i*cell.bsX + j].value1 + 0.5 * layer1[i*cell.bsX + j].value1;
			cur_layer[i*cell.bsX + j].value2 = 0.5 * cur_layer[i*cell.bsX + j].value2 + 0.5 * layer1[i*cell.bsX + j].value2;
			cur_layer[i*cell.bsX + j].value3 = 0.5 * cur_layer[i*cell.bsX + j].value3 + 0.5 * layer1[i*cell.bsX + j].value3;
			cur_layer[i*cell.bsX + j].value4 = 0.5 * cur_layer[i*cell.bsX + j].value4 + 0.5 * layer1[i*cell.bsX + j].value4;
		}
	}
}



TCell Psi(TCell work_cell, TCell gugonio, double x, double y) // work_cell - cell we are working with
{

	if ( y < h1) // we are in 1 layer
	{
		work_cell.pressure = pressure1;
		work_cell.ro = ro1;
		work_cell.speedU = 0.0 ;
		work_cell.speedV = 0.0;
		work_cell.eps =  work_cell.pressure / ( (gamma1 - 1) * work_cell.ro );
		work_cell.value1 = work_cell.ro; 
		work_cell.value2 = 0.0;
		work_cell.value3 = 0.0;
		work_cell.value4 = work_cell.ro * work_cell.eps;
	}

	else if ( y < getH2(x,y) )
		{
			work_cell.pressure = pressure2;
			work_cell.ro = ro2;
			work_cell.speedU = 0.0 ;
			work_cell.speedV = 0.0;
			work_cell.eps =  work_cell.pressure / ( (gamma1 - 1) * work_cell.ro );
			work_cell.value1 = work_cell.ro; 
			work_cell.value2 = 0.0;
			work_cell.value3 = 0.0;
			work_cell.value4 = work_cell.ro * work_cell.eps;
	
		}
	
	else if ( y < h3 )
	{
		work_cell.pressure = pressure3;
		work_cell.ro = ro3;
		work_cell.speedU = 0.0 ;
		work_cell.speedV = 0.0;
		work_cell.eps =  work_cell.pressure / ( (gamma1 - 1) * work_cell.ro );
		work_cell.value1 = work_cell.ro; 
		work_cell.value2 = 0.0;
		work_cell.value3 = 0.0;
		work_cell.value4 = work_cell.ro * work_cell.eps;
	}
	else
	{
		work_cell = gugonio;
		work_cell.value1 = gugonio.ro; 
		work_cell.value2 = gugonio.ro * gugonio.speedU;
		work_cell.value3 = gugonio.ro * gugonio.speedV;
		work_cell.value4 = gugonio.ro * (gugonio.eps + (gugonio.speedU * gugonio.speedU + gugonio.speedV * gugonio.speedV) / 2.0);
	}
	return work_cell;
}




int countExchang(int Nx, int Ny, int npX, int npY)
{
	int count = (npX - 1) * Ny + ( npY - 1 ) * Nx;
	return count;
}


int separate(int np, int Nx, int Ny)
{
	int npY, npX, temp_npY, temp_npX;
	int i;
	int th = (int)(sqrt((double)(np - 1)));

	for (i = th; i > 0; --i)
	{
		if ((np-1) % i == 0)
		{
			npY = i;
			npX = (np - 1) / npY;
			break;
		}
	}

	int count = countExchang( Nx, Ny, npX, npY);
	

	//printf ("Rus_count = %d\n", count);
	temp_npY = 1;                // horizontal line
	temp_npX = np - 1;

	int temp_count = countExchang(Nx, Ny, temp_npX, temp_npY);
	//printf ("horiz_count = %d\n", temp_count);
	if ( temp_count <= count )
	{
	//	printf("horizontal_line\n");
		npY = temp_npY;
	}
	
	temp_npY = np - 1;   // vertical line
	temp_npX = 1;
	
	int temp_count2 = countExchang(Nx, Ny, temp_npX, temp_npY);
	//printf ("vertical_count = %d\n", temp_count2);
	if ( temp_count2 <= count && temp_count2 <= temp_count )
	{
	//	printf("vertical_line\n");
		npY = temp_npY;
	}

	return npY;
	
}



void setBorder(TCell * cur, TCell * prev, Cell cell)
{
	int border = cell.border;

	if (border >= 8)// down
	{
		for (int i = 0; i < cell.bsX; i++)
		{
			cur[(cell.bsY - 1 ) * cell.bsX + i] = cur[(cell.bsY - 2 ) * cell.bsX + i];
			cur[(cell.bsY - 1 ) * cell.bsX + i].speedV = -cur[(cell.bsY - 2) * cell.bsX + i].speedV; 
			cur[(cell.bsY - 1 ) * cell.bsX + i].value3  = -cur[(cell.bsY - 2) * cell.bsX + i].value3;
		}
		border -= 8;

	}

	if (border >= 4) // right
	{
		for (int j = 0 ; j < cell.bsY; j++)
		{
			cur[( j * cell.bsX) + (cell.bsX-1)] = cur[( j * cell.bsX) + (cell.bsX-2)];
			cur[( j * cell.bsX) + (cell.bsX-1)].speedU = -cur[( j * cell.bsX) + (cell.bsX-2)].speedU;
			cur[( j * cell.bsX) + (cell.bsX-1)].value2 = -cur[( j * cell.bsX) + (cell.bsX-2)].value2;
		}
		border -= 4;
	}

	if (border >= 2) // up
	{
		for (int i = 0; i < cell.bsX; i++)
		{
			cur[ ( 0 * cell.bsX ) + i] = cur[ ( 1 * cell.bsX ) + i];
			cur[ ( 0 * cell.bsX ) + i].speedV = -cur[(1 * cell.bsX) + i].speedV; 
			cur[ ( 0 * cell.bsX ) + i].value3 = -cur[ ( 0 * cell.bsX ) + i].value3;
		}
		border -= 2;
	}

	if (border >= 1) // left
	{
		for (int j = 0 ; j < cell.bsY; j++)
		{
			cur[ j * cell.bsX + 0] = cur[ (j * cell.bsX) + 1];
			cur[ j * cell.bsX + 0].speedU = -cur[( j * cell.bsX) + 1].speedU;
			cur[ j * cell.bsX + 0].value2 = -cur[( j * cell.bsX) + 1].value2;

		}
		border -= 1;
	}
}



void printResult(TCell *cur, Cell cell)
{
	int start_i, start_j, end_i, end_j; 
	int border = cell.border;
	start_j = start_i = 1;
	end_i = cell.bsX-1;  // т.к. cell.bsX - 1 - мы попадаем в границу, а нам, необходимо начать раньше.
	end_j = cell.bsY-1;

	if (border >= 7)
	{
		end_j--;
		border -= 7;

	}

	if (border >= 4)
	{
		end_i--;
		border -= 4;
	}

	if (border >= 2)
	{
		start_j++;
		border -= 2;
	}

	if (border >= 1)
	{
		start_i++;
		border -= 1;
	}


}



/*
void layerInit(TCell *cur, TCell gugonio, Cell cell, double dx, double dy)
{
		double x, y;
		for (int i = 1; i < cell.bsY - 1; i++)
		{
			for (int j = 1; j < cell.bsX - 1; j++)
			{
				x = cell.start_x + (j-1) * dx;
				y = cell.start_y + (i-1) * dy;
				cur[i* cell.bsX + j] = Psi(cur, gugonio,x,y);
			}
		}
		// заполним граничные условия!
		setBorder(cur, cell);


}*/


Cell cellInit(Cell cell, int nX, int nY, int idx, int npX, int npY, int np)
{
	// define Cell size;
	int bsX = nX / npX;		
	int bsY = nY / npY;
		
	cell.border = 0;
	
	if (idx % npX == 0)
	{
		bsX += (nX % npX);
	}

	if (idx >=  np - npX)
	{
		bsY += (nY % npY);
	}
	bsX += 2;
	bsY += 2;

	// Has no left neighbour
	if (idx % npX == 1 || npX == 1)
	{
		cell.border += 1;
	}

	
	// Has no Right neighbour
	if (idx % npX == 0)
	{
		cell.border += 4;
	}

	// Has no Up neighbour
	if (idx - npX <= 0)
	{
		cell.border += 2;
	}

	// Has no Down neighbour
	if (idx + npX >= np)	// down
	{
		cell.border += 8;
	}
	
	cell.bsX = bsX;
	cell.bsY = bsY;
	return cell;
}


void copyCol(TCell *cur, TCell * buff, Cell cell, int j_to_copy) // j - отвечает за ту колонку, которую хотим скопировать (левую или правую)
{
	int i;
	for ( i = 0; i < cell.bsY ;  i++)
	{
		buff[i] = cur[i*cell.bsX + j_to_copy];
	}
}

void copyRow(TCell *cur, TCell * buff, Cell cell, int i_to_copy) // j - отвечает за ту строку, которую хотим скопировать (левую или правую)
{
	int j;
	for ( j = 0; j < cell.bsX ;  j++)
	{
		buff[j] = cur[i_to_copy* cell.bsX + j];
	}
}


void insertCol(TCell *cur, TCell * buff, Cell cell, int j_to_insert)
{
	int i;
	for ( i = 0; i < cell.bsY ;  i++)
	{
		cur[i*cell.bsX + j_to_insert] = buff[i];
	}
}

void insertRow(TCell *cur, TCell * buff, Cell cell, int i_to_insert)
{
	int j;
	for ( j = 0; j < cell.bsX ;  j++)
	{
		cur[i_to_insert* cell.bsX + j] = buff[j];
	}
}

void sendReceive(TCell *cur, Cell cell, int npX)
{
	int np, idx;
	cell.border = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &idx);


	// Send to left
	if (idx % npX != 1 && npX != 1 )
	{
		//printf("%d   %d   send to left\n", idx, cell.bsY * sizeof( struct TCell) );
		TCell* buf_StL = (TCell*) malloc (cell.bsY * sizeof( struct TCell));
		copyCol(cur, buf_StL, cell, 1); // заполнили буфер;
		MPI_Send( buf_StL, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx - 1, 0, MPI_COMM_WORLD );
	}

	//Send to down
	if (idx + npX < np)	// down
	{
		//printf("%d   %d   send to down\n", idx, cell.bsX * sizeof( struct TCell));
		TCell* buf_ds = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		copyRow(cur, buf_ds, cell, cell.bsY - 2);

		MPI_Send( buf_ds, cell.bsX * sizeof(struct TCell), MPI_CHAR, idx + npX, 3, MPI_COMM_WORLD );

	}

	//Receive from right
	if (idx % npX != 0)
	{
		//printf("%d  %d   --Receive from right\n ",idx, cell.bsY * sizeof(struct TCell)  );
		TCell* buf_rr = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		MPI_Status status1;
		MPI_Recv(buf_rr, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx + 1, 0, MPI_COMM_WORLD, &status1);
		insertCol(cur, buf_rr, cell, cell.bsX - 1);
	}

	//Receive from up
	if (idx - npX > 0)	//up
	{
		//printf("%d    %d    --Receive from up\n ", idx, cell.bsX * sizeof(struct TCell));
		TCell* buf_ur = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		MPI_Status status4;
		MPI_Recv(buf_ur, cell.bsX * sizeof( struct TCell), MPI_CHAR, idx - npX, 3, MPI_COMM_WORLD, &status4);
		insertRow(cur, buf_ur, cell, 0);
	}

	//Send to Up
	if (idx - npX > 0)
	{
		//printf("%d    %d   send to up\n", idx, cell.bsX * sizeof(struct TCell));
		TCell* buf_StU = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		copyRow(cur, buf_StU, cell, 1);
		MPI_Send( buf_StU, cell.bsX * sizeof(struct TCell), MPI_CHAR, idx - npX, 2, MPI_COMM_WORLD );
	}

	//Send to Right
	if (idx % npX != 0)
	{
		//printf("%d    %d   send to right\n", idx, cell.bsY * sizeof(struct TCell));
		TCell* buf_StR = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		copyCol(cur, buf_StR, cell, cell.bsX - 2); // заполнили буфер
		MPI_Send( buf_StR, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx + 1, 1, MPI_COMM_WORLD );
	}

	//Receive from down
	if (idx + npX < np)
	{
		//printf("%d   %d   --Receive from down\n ", idx, cell.bsX * sizeof(struct TCell));
		TCell* buf_dr = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		MPI_Status status3;
		MPI_Recv(buf_dr, cell.bsX * sizeof (struct TCell), MPI_CHAR, idx + npX, 2, MPI_COMM_WORLD, &status3);
		insertRow(cur, buf_dr, cell, cell.bsY - 1 );
	}

	//Receive from left
	if (idx % npX != 1 && npX != 1)
	{
		//printf("%d   %d   --Receive from left\n ", idx, cell.bsY * sizeof(struct TCell) );
		TCell* buf_lr = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		MPI_Status status2;
		MPI_Recv(buf_lr, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx - 1, 1, MPI_COMM_WORLD, &status2);
		insertCol(cur, buf_lr, cell,0);
	}
}



void updateLayer(TCell *cur, Cell cell)
{
	for (int i = 1; i < cell.bsY - 1; i++)
	{
		for (int j = 1; j < cell.bsX - 1; j++)
		{
			cur[i*cell.bsX + j ].ro = cur[i*cell.bsX + j ].value1;
			cur[i*cell.bsX + j ].speedU = cur[i*cell.bsX + j ].value2 / cur[i*cell.bsX + j ].ro;
			cur[i*cell.bsX + j ].speedV = cur[i*cell.bsX + j ].value3 / cur[i*cell.bsX + j ].ro;
			cur[i*cell.bsX + j ].eps = cur[i*cell.bsX + j ].value4 / cur[i*cell.bsX + j ].ro - ( cur[i*cell.bsX + j ].speedU * cur[i*cell.bsX + j ].speedU + cur[i*cell.bsX + j ].speedV * cur[i*cell.bsX + j ].speedV) / 2.0;
			cur[i*cell.bsX + j ].pressure = (gamma1 - 1.0) * cur[i*cell.bsX + j ].eps * cur[i*cell.bsX + j ].ro;
		}
	}
}


void printAllCell(TCell * cur, Cell cell)
{
	for (int i = 0; i < cell.bsY; i++)
	{
		printf ("%d ", i);
		for (int j = 0; j < cell.bsX; j++)
		{
			printf ("%d  u = %f, v = %f ||", j, cur[i*cell.bsX + j ].speedU, cur[i*cell.bsX + j ].speedV ) ;
		}
		printf("\n");
	}
}

TCell* Gugonio(TCell * cur, Cell cell, int wave_front_index)
{
	int i = wave_front_index;			// номер строки которую надо заполнить
	int index, index2;
	double a, b, c;
	double D;
	TCell *buf = new TCell[cell.bsX];
	for (int t = i+1 ; t < i; t++)          // Присвоение всем ячейкам от волны до конца блока НАЧАЛЬНЫХ значений
	{
		for (int j = 0; j < cell.bsX; j++)
		{
			cur[i*cell.bsX + j].pressure = 1.0;
			cur[i*cell.bsX + j].ro = 1.0;
			cur[i*cell.bsX + j].speedU = 0;
			cur[i*cell.bsX + j].speedV = 0;
			cur[i*cell.bsX + j].value1 = 0;
			cur[i*cell.bsX + j].value2 = 0;
			cur[i*cell.bsX + j].value3 = 0;
			cur[i*cell.bsX + j].value4 = 0;
			cur[i*cell.bsX + j].eps = cur[i*cell.bsX + j].pressure / (gamma1 - 1) / cur[i*cell.bsX + j].ro;
		}
	}
//	setBorder(cur, cell);

	//-------------------------------------Гюгонио------------------
	for ( int j = 0 ; j < cell.bsX; j++  )
	{
		index = i*cell.bsX + j ;
		index2 = (i+1) * cell.bsX + j;
		
		D = M * sqrt(gamma1* (gamma1 - 1) * cur[index2].eps);
		a = cur[index2].ro * (D - cur[index2].speedV);
		b = cur[index2].ro * pow((D - cur[index2].speedV),2.0) + cur[index2].pressure;
		c = cur[index2].eps + (cur[index2].pressure / cur[index2].ro ) +  pow((D - cur[index2].speedV),2.0) / 2.0;
		
		cur[index].speedU = cur[index2].speedU;
		cur[index].ro = ( b * gamma1 + sqrt(b*b * gamma1*gamma1 - 2.0* a * c * (gamma1 * gamma1 - 1.0)) ) / ( 2.0 * c *( gamma1 - 1.0));
		cur[index].speedV = D - a / cur[index].ro;
		cur[index].pressure = b - a*a / cur[index].ro;
		cur[index].eps = cur[index].pressure / ( (gamma1 - 1.0 )* cur[index].ro );
		
		buf[j] = cur[index];


//-----------------------------------------------------------------------------

		for (int t = 0 ; t < i; t++)          // Присвоение всем ячейкам от волны до начала блока Значения ГЮГОНИО
		{
			for (int j = 0; j < cell.bsX; j++)
			{
				index = i*cell.bsX + j ;
				cur[t*cell.bsX + j] = cur[index] ; // Так можем присваивать - вроде да со структрурами такое возможно.
			}
		}
	}	
	return buf;
}
