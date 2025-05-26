#ifndef _GMST_
#define _GMST_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Frac.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula el Tiempo Sideral Medio de Greenwich (GMST) en radianes para una fecha dada.
 * @param Mjd_UT1 Fecha en Tiempo Universal UT1 (Modified Julian Date)
 * @return Tiempo Sideral Medio de Greenwich en radianes (rango [0, 2π])
 */
double gmst(double Mjd_UT1);

#endif