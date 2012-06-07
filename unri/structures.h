#ifndef TCELL
#define TCELL 

struct Cell
{
	int i;
	int j;
	int bsX;
	int bsY;
	double start_x;
	double start_y;
	char border;

};


struct TCell
{
	double ro;
	double pressure;
	double speedU;
	double speedV;
	double eps;     // pressure[ij]/( (gamma - 1) * ro )
	double value1;  // ro
	double value2;  // ro * speedU
	double value3;  // ro * speedV
	double value4;  // E
};

struct TMessage
{
	TCell gcell;
	double start_x;
	double start_y;
};

#endif
