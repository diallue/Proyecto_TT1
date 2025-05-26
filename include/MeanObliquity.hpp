#ifndef _MEANOBLIQUITY_
#define _MEANOBLIQUITY_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

/**
 * @details Calcula la oblicuidad media de la eclíptica para una fecha dada.
 * @param Mjd_TT Fecha en días julianos modificados (Tiempo Terrestre).
 * @return Oblicuidad media en radianes.
 */
double MeanObliquity(double Mjd_TT);

#endif