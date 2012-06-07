/*
void sendAll(TCell *cur, Cell cell, int npX)
{
	int np, idx;
	cell.border = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &idx);
	
	// Send to left
	if (idx % npX != 1)
	{
		cell.border += 1;
		TCell* buf_StL = (TCell*) malloc (cell.bsY * sizeof( struct TCell));
		copyCol(cur, buf_StL, cell, 1); // заполнили буфер;
		MPI_Send( buf_StL, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx - 1, 0, MPI_COMM_WORLD );
	}

	//Send to Right
	if (idx % npX != 0)
	{
		cell.border += 4;
		TCell* buf_StR = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		copyCol(cur, buf_StR, cell, cell.bsX - 2); // заполнили буфер
		MPI_Send( buf_StR, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx + 1, 1, MPI_COMM_WORLD );
	}

	//Send to Up
	if (idx - npX > 0)
	{
		cell.border += 2;
		TCell* buf_StU = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		copyRow(cur, buf_StU, cell, 1);
		MPI_Send( buf_StU, cell.bsX * sizeof(struct TCell), MPI_DOUBLE, idx - npX, 2, MPI_COMM_WORLD );
	}

	//Send to down
	if (idx + npX < np)	// down
	{
		cell.border += 7;
		TCell* buf_ds = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		copyRow(cur, buf_ds, cell, cell.bsY - 2);
		MPI_Send( buf_ds, cell.bsX * sizeof(struct TCell), MPI_CHAR, idx + npX, 3, MPI_COMM_WORLD );
	}
}


void receiveAll(TCell *cur, Cell cell, int npX)
{
	int np, idx;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &idx);

	//Receive from left
	if (idx % npX != 1)
	{
		TCell* buf_lr = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		MPI_Status status2;
		MPI_Recv(buf_lr, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx - 1, 1, MPI_COMM_WORLD, &status2);
		insertCol(cur, buf_lr, cell,0);
	}

	//Receive from right
	if (idx % npX != 0)
	{
		TCell* buf_rr = (TCell*) malloc (cell.bsY * sizeof(struct TCell));
		MPI_Status status1;
		MPI_Recv(buf_rr, cell.bsY * sizeof(struct TCell), MPI_CHAR, idx + 1, 0, MPI_COMM_WORLD, &status1);
		insertCol(cur, buf_rr, cell, cell.bsX - 1);
	}

	//Receive from up
	if (idx - npX > 0)	//up
	{
		TCell* buf_ur = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		MPI_Status status4;
		MPI_Recv(buf_ur, cell.bsX * sizeof( struct TCell), MPI_CHAR, idx - npX, 3, MPI_COMM_WORLD, &status4);
		insertRow(cur, buf_ur, cell, 0);
	}

	//Receive from down
	if (idx + npX < np)
	{
		TCell* buf_dr = (TCell*) malloc (cell.bsX * sizeof(struct TCell));
		MPI_Status status3;
		MPI_Recv(buf_dr, cell.bsX * sizeof (struct TCell), MPI_CHAR, idx + npX, 2, MPI_COMM_WORLD, &status3);
		insertRow(cur, buf_dr, cell, cell.bsY - 1 );
	}

}
*/