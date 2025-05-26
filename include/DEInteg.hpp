#ifndef _DEINTEG_
#define _DEINTEG_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"
#include <tuple>
#include <iostream>
#include <cmath>
#include <limits>

/**
 * @details Integra ecuaciones diferenciales ordinarias (ODE) usando un método de diferencias hacia atrás de paso variable (método de Gear).
 * @param f Función que define la ODE (dy/dt = f(t,y)).
 * @param t Tiempo inicial.
 * @param tout Tiempo final deseado.
 * @param relerr Tolerancia relativa de error.
 * @param abserr Tolerancia absoluta de error.
 * @param n_eqn Número de ecuaciones.
 * @param y Vector de estado inicial (n_eqn x 1).
 * @return Vector de estado en el tiempo tout.
 * @throw Error si los parámetros son inválidos o se accede a índices fuera de rango.
 */
Matrix DEInteg(Matrix f(double t, Matrix z), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif