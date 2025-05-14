#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\R_z.hpp"
#include "..\include\R_y.hpp"
#include <cmath>
#include <tuple>

Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif