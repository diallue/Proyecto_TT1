#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include <cmath>
#include <tuple>

tuple<Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix> JPL_Eph_DE430(double Mjd_TDB);

#endif