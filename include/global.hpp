#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include <cmath>
#include <cstring>

typedef struct{
	double Mjd_UTC,Mjd_TT;
	int n,m,sun,moon,planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Matrix obs;

void eop19620101(int c);
void GGM03S(int n);
void DE430Coeff(int row,int column);
void AuxParamLoad();
void GEOS3(int nobs);

#endif