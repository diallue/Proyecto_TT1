#ifndef _ECCANOM_
#define _ECCANOM_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <iostream>

using namespace std;

/**
 * @details Calcula la anomalía excéntrica (E) a partir de la anomalía media (M) y la excentricidad (e).
 * @details Utiliza el método de Newton-Raphson para resolver la ecuación de Kepler.
 * @param M Anomalía media (en radianes).
 * @param e Excentricidad de la órbita (0 <= e < 1).
 * @return Anomalía excéntrica (E, en radianes).
 * @throw Error si no converge tras 15 iteraciones.
 */
double EccAnom(double M, double e);

#endif