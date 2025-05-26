#ifndef _AZELPA_
#define _AZELPA_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>
#include <vector>

/**
 * @details Calcula azimut, elevación y sus matrices de derivadas parciales para un vector de posición dado.
 * @param s Vector de posición en coordenadas topocéntricas [x, y, z] (metros)
 * @return Tupla con:
 *         [0] Azimut (radianes) [0, 2pi)
 *         [1] Elevación (radianes)
 *         [2] Matriz de derivadas parciales del azimut respecto a s (1x3)
 *         [3] Matriz de derivadas parciales de la elevación respecto a s (1x3)
 */
tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s);

#endif