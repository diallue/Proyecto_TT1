#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "..\include\matrix.hpp"
#include <cmath>
#include <vector>

/**
 * @details Calcula los polinomios de Legendre asociados normalizados y sus derivadas.
 * @param n Grado máximo del polinomio
 * @param m Orden máximo del polinomio
 * @param fi Latitud geocéntrica en radianes
 * @param[out] pnm Matriz (n+1)x(m+1) de polinomios de Legendre normalizados
 * @param[out] dpnm Matriz (n+1)x(m+1) de derivadas de los polinomios
 */
tuple<Matrix, Matrix> Legendre(int n, int m, double fi);

#endif