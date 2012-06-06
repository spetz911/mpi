#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

//#define MSMPI_NO_DEPRECATE_20

#define gamma 5.0/3.0

#define ERR_BADORDER    255
#define TAG_INIT      31337
#define TAG_RESULT       42

struct Node;

double calcEps(double E, double u, double v);
double calcP(double rho, double eps);

double getP(Node node);
double getU(Node node);
double getV(Node node);
double getRHO(Node node);

double max(double a, double b){return a>=b?a:b;}


struct Node{	//узел сетки
	double rho;
    double urho;
    double vrho;
    double rhoE;
};

Node createNode(double rho, double u, double v, double p)
{
	Node res;
    res.rho = rho;
    res.urho = rho*u;
    res.vrho = rho*v;
    res.rhoE = rho*( (p/(rho*(gamma-1))) + (u*u + v*v)/2.0);

    return res;
}

void build_derived_type(Node*, MPI_Datatype*);


//==========================================================================================
class SoDE{
public:
	int idMainProc;	//id of main process
	int idCurrProc;	//id of curr process
	int usedProcCount;	//number of nodes, processes matrix
	int maxProcCount;	//number of available nodes

	int K; //max iteration

	int Nx, Ny;	//размеры кластерной сетки
	double dx, dy, dt; //Шаги сетки

	double edgeXLeft, edgeXRight, edgeYBefore, edgeYAfter; //границы области

	// размеры сетки глобальные и локальные
	int sizeGlobX;
	int sizeGlobY;
	int sizeLocX;
	int sizeLocY;
	Node** currLayer; //node[sizeLocX][sizeLocY]
	Node** prevLayer; //node[sizeLocX][sizeLocY]

//	borders 
	Node* borderLeft;
	Node* borderRight;
	Node* borderBefore;
	Node* borderAfter;

	Node prevLayerGetNode(int x, int y);	//Gets currLayer[x][y] if needed node is in or get value from borders
	void synch();	//synchronize borders

	bool isMain();
	
	void setStartCond();
	void setStartCond_test();
	void setBorderCond();

	void currL2prevL();
	
	int locX2globX(int locX);
	int locY2globY(int locY);
	int globX2locX(int globX);
	int globY2locY(int globY);
	int getProcID(int globX, int globY);

	void build_derived_type(Node* indata, MPI_Datatype* message_type_ptr);
	
	// return true if element has according neighbor
	bool hasLeft();
	bool hasRight();
	bool hasBefore();
	bool hasAfter();

	/*
	double getP(Node node);
	double getU(Node node);
	double getV(Node node);
	double getRHO(Node node);*/

	void collectOnMain(); // collect data from all nodes on main process
	void print2disck(Node** res, int sizeX, int sizeY);

	//------------------------------------------- Method -------------------------------------
	//double calcEps(double E, double u, double v);
	//double calcP(double rho, double eps);

	double Dx(double i, int j);
	double Dy(int i, double j);
	double Q1_rho_u(double i, int j);
	double Q1_rho_v(int i, double j);
	double Q2_rho_u2_p(double i, int j);
	double Q2_rho_u_v(int i, double j);
	double Q3_rho_u_v(double i, int j);
	double Q3_rho_v2_p(int i, double j);
	double Q4_rho_u_E_p_u(double i, int j);
	double Q4_rho_v_E_p_v(int i, double j);
	//------------------------------------------- END Method ---------------------------------

public:
	SoDE();			//count sizes, initialize layers
	void solve();
	void print();

};
//==========================================================================================

void SoDE::build_derived_type(Node* indata, MPI_Datatype* message_type_ptr)
{
  int block_lengths[4];
  MPI_Datatype typelist[4];
  MPI_Aint displacements[4];
  MPI_Aint addresses[5];

  /* Создает производный тип данных, содержащий два элемента float и один int */
  /* Сначала нужно определить типы элементов */
  
  typelist[0] = MPI_DOUBLE;
  typelist[1] = MPI_DOUBLE;
  typelist[2] = MPI_DOUBLE;
  typelist[3] = MPI_DOUBLE;

  /* Определить количество элементов каждого типа */

  block_lengths[0] = block_lengths[1] = block_lengths[2] = block_lengths[3] = 1;

	//double rho;
	//double urho;
	//double vrho;
	//double rhoE;

  /* Вычислить смещения элементов относительно indata */

  MPI_Address(indata, &addresses[0]);
  MPI_Address(&(indata->rho), &addresses[1]);
  MPI_Address(&(indata->urho), &addresses[2]);
  MPI_Address(&(indata->vrho), &addresses[3]);
  MPI_Address(&(indata->rhoE), &addresses[4]);

  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];

  /* Создать производный тип */
  MPI_Type_struct(4, block_lengths, displacements, typelist, message_type_ptr);

  /* Зарегистрировать его для использования */
  MPI_Type_commit(message_type_ptr);

} /* Build_derived_type */ 

double getP(Node node)
{
    double e = calcEps(node.rhoE/node.rho,getU(node),getV(node));
    return calcP(node.rho,e);
}
double getU(Node node){ return node.urho/node.rho; }
double getV(Node node){ return node.vrho/node.rho; }
double getRHO(Node node){ return node.rho; }


bool SoDE::hasLeft(){ return idCurrProc%Nx != 0;}
bool SoDE::hasRight(){ return (idCurrProc+1)%Nx !=0;}
bool SoDE::hasBefore(){ return idCurrProc >= Nx;}
bool SoDE::hasAfter(){ return idCurrProc < usedProcCount-Nx;}

bool SoDE::isMain(){return true;}

int SoDE::getProcID(int globX, int globY)
{
	int yPos = (globY/(sizeGlobY/Ny));
	if(yPos==Ny)
		yPos--;

	
	int xPos = (globX/(sizeGlobX/Nx));
	if(xPos==Nx)
		xPos--;
	return yPos*Nx + xPos;
}

void SoDE::print()
{
	printf("Process %d:\nNx = %d, Ny = %d, sizeLocX=%d, sizeLocY=%d\n", idCurrProc, Nx, Ny, sizeLocX, sizeLocY);

	if(idCurrProc < usedProcCount)
	{
		//printf("%d\\%d", idCurrProc%Nx!=0?1:0, idCurrProc>=Nx?1:0);
		if(idCurrProc>=Nx)
		{
			printf("\t");
			for(int i=0;i<sizeLocX;i++)
				printf("%.2lf ", borderBefore[i].rho);
		}
		printf("\n");

		for(int j=0;j<sizeLocY;j++)
		{
			if(idCurrProc%Nx !=0)
				printf("%.2lf|", borderLeft[j].rho);
			printf("\t");
			for(int i=0;i<sizeLocX;i++)
			{
				printf("%.2lf ", prevLayer[i][j].rho);
			}
			if( (idCurrProc+1)%Nx !=0 )
				printf("\t|%.2lf", borderRight[j].rho);
			printf("\n");
		}
	
		if(idCurrProc < usedProcCount-Nx)
		{
			printf("\t");
			for(int i=0;i<sizeLocX;i++)
				printf("%.2lf ", borderAfter[i].rho);
			printf("\n");
		}
	}
	else
	{
		printf("Not used\n");
	}
}

void SoDE::setStartCond_test()
{
	for(int j=0;j<sizeLocY;j++)
	{
		for(int i=0;i<sizeLocX;i++)
		{
			currLayer[i][j].rho = locX2globX(i)*100 + locY2globY(j);
		}
	}
}

int SoDE::locX2globX(int locX)
{
	return (idCurrProc%Nx)*(sizeGlobX/Nx) + locX;
	return (idCurrProc*sizeLocX)%sizeGlobX + locX;
}
int SoDE::locY2globY(int locY)
{
	return (idCurrProc/Nx)*(sizeGlobY/Ny) + locY;
}

int SoDE::globX2locX(int globX){ return globX%(sizeGlobX/Nx); }
int SoDE::globY2locY(int globY){ return globY%(sizeGlobY/Ny); }

Node SoDE::prevLayerGetNode(int x, int y)
{
	if(x>=0 && x<sizeLocX)
	{
		if(y>=0 && y<sizeLocY)
			return prevLayer[x][y];
		if(y==-1)
			return borderBefore[x];
		if(y==sizeLocY)
			return borderAfter[x];
	}
	if(x==-1)
		return borderLeft[y];
	else
		return borderRight[y];
}

void SoDE::currL2prevL()
{
	for(int i=0;i<sizeLocX;i++)
	{
		for(int j=0;j<sizeLocY;j++)
		{
			prevLayer[i][j] = currLayer[i][j];
		}
	}
	return;
}

void SoDE::synch(){
	MPI_Barrier(MPI_COMM_WORLD);
	if(idCurrProc < usedProcCount)
	{
		MPI_Datatype MY_NODE;
		build_derived_type(prevLayer[0], &MY_NODE);
		if(idCurrProc >= Nx)
		{
			for(int i=0; i<sizeLocX;i++)
				MPI_Send(&prevLayer[i][0], 1, MY_NODE, idCurrProc-Nx, i, MPI_COMM_WORLD);
			for(int i=0; i<sizeLocX;i++)
				MPI_Recv(&borderBefore[i], 1, MY_NODE, idCurrProc-Nx, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	
		if(idCurrProc < usedProcCount-Nx)
		{
			for(int i=0; i<sizeLocX;i++)
				MPI_Send(&prevLayer[i][sizeLocY-1], 1, MY_NODE, idCurrProc+Nx, i, MPI_COMM_WORLD);
			for(int i=0; i<sizeLocX;i++)
				MPI_Recv(&borderAfter[i], 1, MY_NODE, idCurrProc+Nx, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if(idCurrProc%Nx != 0)
		{
			for(int i=0; i<sizeLocY;i++)
				MPI_Send(&prevLayer[0][i], 1, MY_NODE, idCurrProc-1, i, MPI_COMM_WORLD);
			for(int i=0; i<sizeLocY;i++)
				MPI_Recv(&borderLeft[i], 1, MY_NODE, idCurrProc-1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//MPI_Send(&prevLayer[0], sizeLocY, MY_NODE, idCurrProc-1, 0, MPI_COMM_WORLD);
			//MPI_Recv(&borderLeft[0], sizeLocY, MY_NODE, idCurrProc-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if( (idCurrProc+1)%Nx !=0 )
		{
			for(int i=0; i<sizeLocY;i++)
				MPI_Send(&prevLayer[0][i], 1, MY_NODE, idCurrProc+1, i, MPI_COMM_WORLD);
			for(int i=0; i<sizeLocY;i++)
				MPI_Recv(&borderRight[i], 1, MY_NODE, idCurrProc+1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//MPI_Send(&prevLayer[sizeLocX-1], sizeLocY, MY_NODE, idCurrProc+1, 0, MPI_COMM_WORLD);
			//MPI_Recv(&borderRight[0], sizeLocY, MY_NODE, idCurrProc+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
}


//-------------------------------------------------------------------------
//========================== METHOD =======================================
//-------------------------------------------------------------------------


double calcEps(double E, double u, double v){ return E - (u*u + v*v)/2.0; }
double calcP(double rho, double eps){ return (gamma - 1.0)*rho*eps; }

double SoDE::Dx(double i, int j)
{
    int ind = (int) i;

    double u1 = prevLayerGetNode(ind, j).urho/prevLayerGetNode(ind, j).rho;
    double v1 = prevLayerGetNode(ind, j).vrho/prevLayerGetNode(ind, j).rho;
    double E1 = prevLayerGetNode(ind, j).rhoE/prevLayerGetNode(ind, j).rho;
    double e1 = calcEps(E1,u1,v1);

    double u2 = prevLayerGetNode(ind+1, j).urho/prevLayerGetNode(ind+1, j).rho;
    double v2 = prevLayerGetNode(ind+1, j).vrho/prevLayerGetNode(ind+1, j).rho;
    double E2 = prevLayerGetNode(ind+1, j).rhoE/prevLayerGetNode(ind+1, j).rho;
    double e2 = calcEps(E2,u2,v2);

    double a = pow(gamma*(gamma-1)*e1,0.5) + pow(u1*u1+v1*v1,0.5);
    double b = pow(gamma*(gamma-1)*e2,0.5) + pow(u2*u2+v2*v2,0.5);
    return max(a,b);
}

double SoDE::Dy(int i, double j)
{
    int ind = (int) j;

    double u1 = prevLayerGetNode(i, ind).urho/prevLayerGetNode(i, ind).rho;
    double v1 = prevLayerGetNode(i, ind).vrho/prevLayerGetNode(i, ind).rho;
    double E1 = prevLayerGetNode(i, ind).rhoE/prevLayerGetNode(i, ind).rho;
    double e1 = calcEps(E1,u1,v1);

    double u2 = prevLayerGetNode(i, ind+1).urho/prevLayerGetNode(i, ind+1).rho;
    double v2 = prevLayerGetNode(i, ind+1).vrho/prevLayerGetNode(i, ind+1).rho;
    double E2 = prevLayerGetNode(i, ind+1).rhoE/prevLayerGetNode(i, ind+1).rho;
    double e2 = calcEps(E2,u2,v2);

    double a = pow(gamma*(gamma-1)*e1,0.5) + pow(u1*u1+v1*v1,0.5);
    double b = pow(gamma*(gamma-1)*e2,0.5) + pow(u2*u2+v2*v2,0.5);
    return max(a,b);
}

double SoDE::Q1_rho_u(double i, int j)
{
    int ind = (int) i;
    return (prevLayerGetNode(ind, j).urho + prevLayerGetNode(ind+1, j).urho)/2.0 - Dx(i,j)*(prevLayerGetNode(ind+1, j).rho - prevLayerGetNode(ind, j).rho)/2.0;
}

double SoDE::Q1_rho_v(int i, double j)
{
    int ind = (int) j;
    return (prevLayerGetNode(i, ind).vrho + prevLayerGetNode(i, ind+1).vrho)/2.0 - Dy(i,j)*(prevLayerGetNode(i, ind+1).rho - prevLayerGetNode(i, ind).rho)/2.0;
}

double SoDE::Q2_rho_u2_p(double i, int j)
{
    int ind = (int) i;
    double u1 = prevLayerGetNode(ind, j).urho/prevLayerGetNode(ind, j).rho;
    double v1 = prevLayerGetNode(ind, j).vrho/prevLayerGetNode(ind, j).rho;
    double E1 = prevLayerGetNode(ind, j).rhoE/prevLayerGetNode(ind, j).rho;
    double e1 = calcEps(E1,u1,v1);
    double p1 = calcP(prevLayerGetNode(ind, j).rho,e1);

    double u2 = prevLayerGetNode(ind+1, j).urho/prevLayerGetNode(ind+1, j).rho;
    double v2 = prevLayerGetNode(ind+1, j).vrho/prevLayerGetNode(ind+1, j).rho;
    double E2 = prevLayerGetNode(ind+1, j).rhoE/prevLayerGetNode(ind+1, j).rho;
    double e2 = calcEps(E2,u2,v2);
    double p2 = calcP(prevLayerGetNode(ind+1, j).rho,e2);

    return ((u1*u1*prevLayerGetNode(ind, j).rho + p1) + (u2*u2*prevLayerGetNode(ind, j).rho + p2))/2.0 - Dx(i,j)*(prevLayerGetNode(ind, j).urho - prevLayerGetNode(ind, j).urho)/2.0;
}

double SoDE::Q2_rho_u_v(int i, double j)
{
    int ind = (int) j;
    return (prevLayerGetNode(i, ind).urho*prevLayerGetNode(i, ind).vrho/prevLayerGetNode(i, ind).rho + prevLayerGetNode(i, ind+1).urho*prevLayerGetNode(i, ind+1).vrho/prevLayerGetNode(i, ind+1).rho)/2.0 - Dy(i,j)*(prevLayerGetNode(i, ind+1).urho - prevLayerGetNode(i, ind).urho)/2.0;
}

double SoDE::Q3_rho_v2_p(int i, double j)
{
    int ind = (int) j;
    double u1 = prevLayerGetNode(i, ind).urho/prevLayerGetNode(i, ind).rho;
    double v1 = prevLayerGetNode(i, ind).vrho/prevLayerGetNode(i, ind).rho;
    double E1 = prevLayerGetNode(i, ind).rhoE/prevLayerGetNode(i, ind).rho;
    double e1 = calcEps(E1,u1,v1);
    double p1 = calcP(prevLayerGetNode(i, ind+1).rho,e1);

    double u2 = prevLayerGetNode(i, ind+1).urho/prevLayerGetNode(i, ind+1).rho;
    double v2 = prevLayerGetNode(i, ind+1).vrho/prevLayerGetNode(i, ind+1).rho;
    double E2 = prevLayerGetNode(i, ind+1).rhoE/prevLayerGetNode(i, ind+1).rho;
    double e2 = calcEps(E2,u2,v2);
    double p2 = calcP(prevLayerGetNode(i, ind+1).rho,e2);

    return ((v1*v1*prevLayerGetNode(i, ind).rho + p1) + (v2*v2*prevLayerGetNode(i, ind+1).rho + p2))/2.0 - Dy(i,j)*(prevLayerGetNode(i, ind+1).vrho - prevLayerGetNode(i, ind).vrho)/2.0;
}

double SoDE::Q3_rho_u_v(double i, int j)
{
    int ind = (int) i;
    return (prevLayerGetNode(ind, j).urho*prevLayerGetNode(ind, j).vrho/prevLayerGetNode(ind, j).rho + prevLayerGetNode(ind+1, j).urho*prevLayerGetNode(ind+1, j).vrho/prevLayerGetNode(ind+1, j).rho)/2.0 - Dx(i,j)*(prevLayerGetNode(ind+1, j).vrho - prevLayerGetNode(ind, j).vrho)/2.0;
}

double SoDE::Q4_rho_u_E_p_u(double i, int j)
{
    int ind = (int) i;
    double u1 = prevLayerGetNode(ind, j).urho/prevLayerGetNode(ind, j).rho;
    double v1 = prevLayerGetNode(ind, j).vrho/prevLayerGetNode(ind, j).rho;
    double E1 = prevLayerGetNode(ind, j).rhoE/prevLayerGetNode(ind, j).rho;
    double e1 = calcEps(E1,u1,v1);
    double p1 = calcP(prevLayerGetNode(ind, j).rho,e1);

    double u2 = prevLayerGetNode(ind+1, j).urho/prevLayerGetNode(ind+1, j).rho;
    double v2 = prevLayerGetNode(ind+1, j).vrho/prevLayerGetNode(ind+1, j).rho;
    double E2 = prevLayerGetNode(ind+1, j).rhoE/prevLayerGetNode(ind+1, j).rho;
    double e2 = calcEps(E2,u2,v2);
    double p2 = calcP(prevLayerGetNode(ind+1, j).rho,e2);

    return ((prevLayerGetNode(ind, j).rhoE*u1+p1*u1) + (prevLayerGetNode(ind+1, j).rhoE*u2+p2*u2))/2.0 - Dx(i,j)*(prevLayerGetNode(ind+1, j).rhoE - prevLayerGetNode(ind, j).rhoE)/2.0;
}

double SoDE::Q4_rho_v_E_p_v(int i, double j)
{
    int ind = (int) j;
    double u1 = prevLayerGetNode(i, ind).urho/prevLayerGetNode(i, ind).rho;
    double v1 = prevLayerGetNode(i, ind).vrho/prevLayerGetNode(i, ind).rho;
    double E1 = prevLayerGetNode(i, ind).rhoE/prevLayerGetNode(i, ind).rho;
    double e1 = calcEps(E1,u1,v1);
    double p1 = calcP(prevLayerGetNode(i, ind+1).rho,e1);

    double u2 = prevLayerGetNode(i, ind+1).urho/prevLayerGetNode(i, ind+1).rho;
    double v2 = prevLayerGetNode(i, ind+1).vrho/prevLayerGetNode(i, ind+1).rho;
    double E2 = prevLayerGetNode(i, ind+1).rhoE/prevLayerGetNode(i, ind+1).rho;
    double e2 = calcEps(E2,u2,v2);
    double p2 = calcP(prevLayerGetNode(i, ind+1).rho,e2);

    return ((prevLayerGetNode(i, ind).rhoE*v1+p1*v1) + (prevLayerGetNode(i, ind+1).rhoE*v2+p2*v2))/2.0 - Dy(i,j)*(prevLayerGetNode(i, ind+1).rhoE - prevLayerGetNode(i, ind).rhoE)/2.0;
}

//-------------------------------------------------------------------------
//========================== END METHOD ===================================
//-------------------------------------------------------------------------

void SoDE::setBorderCond()
{
	if(!hasLeft())
	{

		double gammaTMP = gamma/(gamma-1.0);

		double M=2;

		double p0 = 1.0;
		double rho0 = 1.0;
		double D = M* sqrt(gamma*p0/rho0);
		double u0 = 0.0 - D;

		double A = (  0.5-gammaTMP   );
		double B = (  gammaTMP*(u0+p0/(rho0*u0))  );
		double C = -(  p0/(rho0*(gamma-1)) + p0/rho0 + u0*u0/2.0   );
		double Diskr = B*B - 4*A*C;

		double u1 = (-B-sqrt(Diskr))/(2.0*A); //+sqrt(Diskr)
		double p1 = rho0*u0*(u0-u1)+p0;
		double rho1 = rho0*u0/u1;
		u1 += D;
		u0 += D;

		for(int j=0; j<sizeLocY;j++)
			currLayer[0][j] = createNode(rho1, u1, 0, p1);

		/*for(int j=0; j<sizeLocY;j++)
			currLayer[0][j] = currLayer[1][j];*/
	}

	if(!hasRight())
	{
		for(int j=0; j<sizeLocY;j++)
			currLayer[sizeLocX-1][j] = currLayer[sizeLocX-2][j];
	}
	
	if(!hasBefore())
	{
		for(int i=0; i<sizeLocX;i++)
			currLayer[i][0] = currLayer[i][1];
	}

	if(!hasAfter())
	{
		for(int i=0; i<sizeLocX;i++)
			currLayer[i][sizeLocY-1] = currLayer[i][sizeLocY-2];
	}
}

void SoDE::setStartCond() // задаем начальные значения shock wave
{
	double gammaTMP = gamma/(gamma-1.0);

	double M=2;

	double p0 = 1.0;
	double rho0 = 1.0;
	double D = M* sqrt(gamma*p0/rho0);
	double u0 = 0.0 - D;

	double A = (  0.5-gammaTMP   );
	double B = (  gammaTMP*(u0+p0/(rho0*u0))  );
	double C = -(  p0/(rho0*(gamma-1)) + p0/rho0 + u0*u0/2.0   );
	double Diskr = B*B - 4*A*C;
	
	/*printf("A=%lf ", A);
	printf("B=%lf ", B);
	printf("C=%lf ", C);
	printf("Diskr=%lf\n", Diskr);*/

	double u1 = (-B-sqrt(Diskr))/(2.0*A); //+sqrt(Diskr)
	double p1 = rho0*u0*(u0-u1)+p0;
	double rho1 = rho0*u0/u1;
	u1 += D;
	u0 += D;

	//printf("rho1=%lf u1=%lf p1=%lf\n", rho1, u1, p1);
	
	/*double eps0 = p0/((gamma-1)*rho0);
	double eps1 = p1/((gamma-1)*rho1);
	if( rho0*(D-u0)!=rho1*(D-u1) || rho0*(D-u0)*(D-u0)+p0!=rho1*(D-u1)*(D-u1)+p1 || eps0+p0/rho0+(D-u0)*(D-u0)/2.0 != eps1+p1/rho1+(D-u1)*(D-u1)/2.0)
	{
		printf("Wrong Gugonio! %lf %lf  |    %lf %lf  |   %lf %lf  \n", rho0*(D-u0), rho1*(D-u1), rho0*(D-u0)*(D-u0)+p0, rho1*(D-u1)*(D-u1)+p1 , eps0+p0/rho0+(D-u0)*(D-u0)/2.0 , eps1+p1/rho1+(D-u1)*(D-u1)/2.0);
	}*/


	for(int i=0; i<sizeLocX;i++)
	{
		for(int j=0;j<sizeLocY;j++)
		{
			//createNode(double rho, double u, double v, double p)

			if(locX2globX(i)<0.1*sizeGlobX)
				currLayer[i][j] = createNode(rho1, u1, 0, p1);
			else
				currLayer[i][j] = createNode(rho0, u0, 0.0, p0);

			/*if( locX2globX(i)<sizeGlobX/10.0 )	// 1
			{
				currLayer[i][j] = createNode(2.0, 0.0, 2.0, 3.0);
			}
			else if( locX2globX(i)<0.9*sizeGlobX )	// 2
			{
				currLayer[i][j] = createNode(2.0-0.01*(locX2globX(i)-sizeGlobX/10.0), 0.0, 2.0, 3.0);//3.0-0.01*(locX2globX(i)-sizeGlobX/10.0));
			}
			else // 3
			{
				currLayer[i][j] = createNode(2.0-(0.8*sizeGlobX)*0.01, 0.0, 1.0, 3.0-0.01*(0.8*sizeGlobX));
			}*/

			/*if(locX2globX(i)<0.5*sizeGlobX)
				currLayer[i][j] = createNode(2.0, 0.0, 2.0, 3.0);
			else
				currLayer[i][j] = createNode(1.0, 0.0, 1.0, 1.0);*/
		}
	}
}

void SoDE::print2disck(Node** res, int sizeX, int sizeY)
{
	FILE *fsP = fopen("res_p.txt","w");
	FILE *fsU = fopen("res_u.txt","w");
	FILE *fsV = fopen("res_v.txt","w");
	FILE *fsRHO = fopen("res_rho.txt","w");

	for(int j=0; j < sizeY; j++)
	{
		for(int i=0; i < sizeX; i++)
		{
			fprintf(fsP,"%.2lf ", getP(res[i][j]) );
			fprintf(fsU,"%.2lf ", getU(res[i][j]) );
			fprintf(fsV,"%.2lf ", getV(res[i][j]) );
			fprintf(fsRHO,"%.2lf ", getRHO(res[i][j]) );
		}
		fprintf(fsP,"\n");
		fprintf(fsU,"\n");
		fprintf(fsV,"\n");
		fprintf(fsRHO,"\n");
	}
	fclose(fsP);
	fclose(fsU);
	fclose(fsV);
	fclose(fsRHO);

	//--------------------------------------------------------------------------------

	FILE *PLOTfsP = fopen("PLOTres_p.txt","w");
	FILE *PLOTfsU = fopen("PLOTres_u.txt","w");
	FILE *PLOTfsV = fopen("PLOTres_v.txt","w");
	FILE *PLOTfsRHO = fopen("PLOTres_rho.txt","w");

	for(int j=0; j < sizeY; j++)
	{
		for(int i=0; i < sizeX; i++)
		{
			fprintf(PLOTfsP,"%.2lf %.2lf %.2lf\n", edgeXLeft+(double)i*dx, edgeYBefore+(double)j*dy, getP(res[i][j]) );
			fprintf(PLOTfsU,"%.2lf %.2lf %.2lf\n", edgeXLeft+(double)i*dx, edgeYBefore+(double)j*dy, getU(res[i][j]) );
			fprintf(PLOTfsV,"%.2lf %.2lf %.2lf\n", edgeXLeft+(double)i*dx, edgeYBefore+(double)j*dy, getV(res[i][j]) );
			fprintf(PLOTfsRHO,"%.2lf %.2lf %.2lf\n", edgeXLeft+(double)i*dx, edgeYBefore+(double)j*dy, getRHO(res[i][j]) );
		}
		fprintf(PLOTfsP,"\n");
		fprintf(PLOTfsU,"\n");
		fprintf(PLOTfsV,"\n");
		fprintf(PLOTfsRHO,"\n");
	}
	fclose(PLOTfsP);
	fclose(PLOTfsU);
	fclose(PLOTfsV);
	fclose(PLOTfsRHO);


	system("gnuplot toplot.txt");

	//----------------------------------

	/*FILE *fsRHO = fopen("res_cons_rho.txt","w");
	FILE *fsRHOU = fopen("res_cons_rhou.txt","w");
	FILE *fsRHOV = fopen("res_cons_rhov.txt","w");
	FILE *fsRHOE = fopen("res_cons_rhoE.txt","w");

	for(int j=0; j < sizeY; j++)
	{
		for(int i=0; i < sizeX; i++)
		{
			fprintf(fsRHO,"%.2lf ", res[i][j].rho );
			fprintf(fsRHOU,"%.2lf ", res[i][j].urho );
			fprintf(fsRHOV,"%.2lf ", res[i][j].vrho );
			fprintf(fsRHOE,"%.2lf ", res[i][j].rhoE );
		}
		fprintf(fsRHO,"\n");
		fprintf(fsRHOU,"\n");
		fprintf(fsRHOV,"\n");
		fprintf(fsRHOE,"\n");
	}
	fclose(fsRHO);
	fclose(fsRHOU);
	fclose(fsRHOV);
	fclose(fsRHOE);*/
}

void SoDE::collectOnMain()
{
	MPI_Datatype MY_NODE;
	build_derived_type(prevLayer[0], &MY_NODE);
	
	if(idCurrProc!=idMainProc)
	{
		for(int i=0; i<sizeLocX;i++)
		{
			for(int j=0;j<sizeLocY;j++)
			{
				MPI_Send(&prevLayer[i][j], 1, MY_NODE, idMainProc, i*1000+j, MPI_COMM_WORLD);
			}
		}
	}

	if(idCurrProc==idMainProc)
	{
		Node** res;
		res = new Node*[sizeGlobX];
		for(int i=0;i<sizeGlobX;i++)
			res[i] = new Node[sizeGlobY];

		for(int i=0; i<sizeGlobX; i++)
		{
			for(int j=0;j<sizeGlobY;j++)
			{
				if(getProcID(i, j)==idMainProc)
					continue;
				MPI_Recv(&res[i][j], 1, MY_NODE, getProcID(i, j), globX2locX(i)*1000+ globY2locY(j), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		for(int i=0;i<sizeLocX;i++)
			for(int j=0; j<sizeLocY;j++)
				res[locX2globX(i)][locY2globY(j)] = prevLayer[i][j];

		// Collected, printing

		this->print2disck(res, sizeGlobX, sizeGlobY);

		/*fprintf(stdout, "Collected result on process %d:", idCurrProc);
		for(int j=0;j<sizeLocY;j++)
		{
			for(int i=0;i<sizeLocX;i++)
				printf("%.2lf ", prevLayer[i][j].rho);
			printf("\n");
		}*/
	}
}

SoDE::SoDE(){
	idMainProc = 0; //TODO get it by local file
	/* Получение номера процесса, вызвавшего функцию */
    MPI_Comm_rank(MPI_COMM_WORLD, &idCurrProc);
	/* Получение количества процессов, учавствующих в выполнении задания */
    MPI_Comm_size(MPI_COMM_WORLD, &maxProcCount);

	K = 200;
	sizeGlobX = 120;	//TODO get it normally
	sizeGlobY = 8;

	edgeXLeft = 0.0;
	edgeXRight = 1.0;
	edgeYBefore = 0.0;
	edgeYAfter = 0.1;


	//----------------------- Nx Ny -----------------------------------------
	//Now count cluster grid sizes knowing global sizes and process count
	Nx = 1;
	Ny = 1;
	
	while(true)
	{
		if( (double)Nx/(double)Ny-(double)sizeGlobX/(double)sizeGlobY > 0)
		{
			if( Nx*(Ny+1) > maxProcCount )
				break;
			if(Ny==sizeGlobY)
				break;
			Ny++;
		}
		
		if( (Nx+1)*Ny > maxProcCount )
			break;
		if(Nx==sizeGlobX)
			break;
		Nx++;
	}
	while( Nx*(Ny+1) <= maxProcCount)
	{
		if(Ny==sizeGlobY)
			break;
		Ny++;
	}
	while( (Nx+1)*Ny <= maxProcCount)
	{
		if(Nx==sizeGlobX)
			break;
		Nx++;
	}

	usedProcCount = Nx*Ny;

	
	//----------------------- END Nx Ny -----------------------------------------
	//Count local grid sizes
	sizeLocX = sizeGlobX/Nx;
	sizeLocY = sizeGlobY/Ny;
	if((idCurrProc+1)%Nx==0)
		sizeLocX = sizeGlobX - (Nx-1)*sizeLocX;
	if(idCurrProc >= usedProcCount-Nx)
		sizeLocY = sizeGlobY - (Ny-1)*sizeLocY;
	//----------------------------------------------

	//шаги сетки
	dx = (edgeXRight-edgeXLeft)/sizeGlobX;
	dy = (edgeYAfter-edgeYBefore)/sizeGlobY;
	dt = dx*dy*dx*dy;
	
	//----------------------------------------------


	// Allocate memory
	currLayer = new Node*[sizeLocX];
	prevLayer = new Node*[sizeLocX];
	for(int i=0;i<sizeLocX;i++)
	{
		currLayer[i] = new Node[sizeLocY];
		prevLayer[i] = new Node[sizeLocY];
	}

	borderLeft = new Node[sizeLocY];
	borderRight = new Node[sizeLocY];
	borderBefore = new Node[sizeLocX];
	borderAfter = new Node[sizeLocX];
	
	for(int i=0;i<sizeLocX;i++)
	{
		borderBefore[i].rho = -1;
		borderBefore[i].rho = -1;
	}
	for(int i=0;i<sizeLocY;i++)
	{
		borderLeft[i].rho = -1;
		borderRight[i].rho = -1;
	}
}

void SoDE::solve()
{
	setStartCond();
	//setStartCond_test();
	currL2prevL();
	synch();

	for(int k = 1; k < 1000; k++)
    {
		int i0 = (hasLeft()?0:1);
		int i1 = sizeLocX - (hasRight()?0:1);
		int j0 = (hasBefore()?0:1);
		int j1 = sizeLocY - (hasAfter()?0:1);

        for(int i = i0; i < i1; i++)
		{
            for(int j = j0; j < j1; j++)
            {
				//currLayer[i][j].rho += 100000;
                currLayer[i][j].rho = prevLayerGetNode(i, j).rho - dt*((Q1_rho_u(i+0.5,j)-Q1_rho_u(i-0.5,j))/dx - (Q1_rho_v(i,j+0.5)-Q1_rho_v(i,j-0.5))/dy);
                currLayer[i][j].urho = prevLayerGetNode(i, j).urho - dt*((Q2_rho_u2_p(i+0.5,j)-Q2_rho_u2_p(i-0.5,j))/dx - (Q2_rho_u_v(i,j+0.5)-Q2_rho_u_v(i,j-0.5))/dy);
                currLayer[i][j].vrho = prevLayerGetNode(i, j).vrho - dt*((Q3_rho_u_v(i+0.5,j)-Q3_rho_u_v(i-0.5,j))/dx - (Q3_rho_v2_p(i,j+0.5)-Q3_rho_v2_p(i,j-0.5))/dy);
                currLayer[i][j].rhoE = prevLayerGetNode(i, j).rhoE - dt*((Q4_rho_u_E_p_u(i+0.5,j)-Q4_rho_u_E_p_u(i-0.5,j))/dx - (Q4_rho_v_E_p_v(i,j+0.5)-Q4_rho_v_E_p_v(i,j-0.5))/dy);
            
			}
		}
		setBorderCond();
		currL2prevL();
		synch();

    }
	collectOnMain();
}


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	//double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffticks;
}


int main(int argc, char *argv[])
{

	clock_t seconds = clock();

	MPI_Init(&argc, &argv);

	SoDE krem_s_sodoy;
	krem_s_sodoy.solve();

	//krem_s_sodoy.print();
	//system("pause");
	MPI_Finalize();
	std::cout<<diffclock(clock(), seconds)<<std::endl;
	//printf ("%lf seconds\n", (time(NULL)-seconds));
	printf("END\n");
	//system("pause");
    return 0;
}
