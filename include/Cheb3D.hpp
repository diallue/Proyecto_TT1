#ifndef _CHEB3D_
#define _CHEB3D_

#include "..\include\matrix.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);

#endif