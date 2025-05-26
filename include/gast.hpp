#ifndef _GAST_
#define _GAST_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\gmst.hpp"
#include <cmath>

/**
 * @details Calcula el Tiempo Sideral Aparente de Greenwich (GAST) en radianes para una fecha dada.
 * @param Mjd_UT1 Fecha en Tiempo Universal UT1 (Modified Julian Date)
 * @return Tiempo Sideral Aparente de Greenwich en radianes (rango [0, 2π])
 */
double gast(double Mjd_UT1);

#endif