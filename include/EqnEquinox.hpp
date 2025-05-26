#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la ecuación de los equinoccios (diferencia entre tiempo sidéreo aparente y medio)
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre (TT)
 * @return Ecuación de los equinoccios en radianes [rad]
 */
double EqnEquinox(double Mjd_TT);

#endif