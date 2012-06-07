#ifndef WORK_WITHCELL
#define WORK_WITHCELL
#include "structures.h"
#include <math.h>
#include <mpi.h>

int separate(int np, int Nx, int Ny); // return number process Y (nY)
int countExchang(int Nx, int Ny, int npX, int npY);
void make_layer(TCell* cur_layer, TCell* layer1, Cell cell);
void layerInit(TCell *cur, TCell gugonio, Cell cell, double dx, double dy, int id_x);
Cell cellInit(Cell cell, int nX, int nY, int idx, int npX, int npY, int np);
void copyCol(TCell *cur, TCell * buff, Cell cell, int j_to_copy);
void copyRow(TCell *cur, TCell * buff, Cell cell, int i_to_copy);
void sendReceive(TCell *cur, Cell cell, int npX);
void sendAll(TCell *cur, Cell cell, int npX);
void receiveAll(TCell *cur, Cell cell, int npX);
void insertCol(TCell *cur, TCell * buff, Cell cell, int j_to_insert);
void insertRow(TCell *cur, TCell * buff, Cell cell, int i_to_insert);
void updateLayer(TCell *cur, Cell cell);
void printResult(TCell *cur, Cell cell);
void printAllCell(TCell *cur, Cell cell);
void setBorder(TCell * cur, TCell *prev, Cell cell);
TCell* Gugonio(TCell * cur, Cell cell, int wave_front_index);
TCell Psi(TCell work_cell, TCell gugonio, double x, double y);



#endif
