#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Legendre.hpp"
#include <cmath>
#include <tuple>

Matrix AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max);

#endif