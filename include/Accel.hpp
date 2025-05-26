#ifndef _ACCEL_
#define _ACCEL_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\global.hpp"
#include <vector>
#include <tuple>

/**
 * @details Calcula la aceleración total que actúa sobre un satélite considerando múltiples perturbaciones
 * @param x Tiempo desde la época de referencia [s]
 * @param Y Vector de estado (6x1) [posición; velocidad] [m, m/s]
 * @return Derivada del vector de estado (6x1) [velocidad; aceleración] [m/s, m/s²]
 */
Matrix Accel(double x, Matrix Y);

#endif