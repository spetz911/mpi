
#include "field.h"

Field::Field(int bx, int by)
{
	this->bx = bx;
	this->by = by;
	
	data.resize(bx + 2);
	for (auto &v : data)
		v.resize(by+2);
}

void
Field::init(double hx, double hy, double ht)
{
	this->hx = hx;
	this->hy = hy;
	this->ht = ht;

	for (int i = 0; i < bx + 2; ++i)
		for (int j = 0; j < by + 2; ++j) {
			data[i][j] = State();
		}
}


void
Field::show_result(const char *fname)
{
	FILE *F = fopen(fname, "w");
	
	if (!F) {
		printf("file not found\n");
		exit(1);
	}
	
	fprintf(F, "%s variable %d %d\n-------------------------------------\n", "rho", bx, by);


	// TODO full result
	
	for (int i = 1; i <= bx; ++i) {
		for (int j = 1; j <= by; ++j)
			fprintf(F, "%f ", data[i][j].p);

		fprintf(F, "\n");
		//fprintf("\n-----   line %d    --------\n", i);
	}



	fclose(F);
}

