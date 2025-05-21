#ifndef _DEINTEG_
#define _DEINTEG_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"
#include "..\include\global.hpp"
#include <tuple>
#include <iostream>
#include <cmath>
#include <limits>

Matrix DEInteg(Matrix f(double t, Matrix z), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif