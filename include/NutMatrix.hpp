#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include <cmath>
#include <tuple>

Matrix NutMatrix(double Mjd_TT);

#endif