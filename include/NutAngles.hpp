#ifndef _NUTANGLES_
#define _NUTANGLES_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <vector>
#include <tuple>

/**
 * @details Calcula los ángulos de nutación (dψ, dε) según el modelo IAU 1980
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre (TT)
 * @return Tupla con:
 *         [0] dψ: Nutación en longitud [rad]
 *         [1] dε: Nutación en oblicuidad [rad]
 */
tuple<double, double> NutAngles(double Mjd_TT);

#endif