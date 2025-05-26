#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Realiza la actualización de medida de un filtro de Kalman.
 * @param x Estado estimado antes de la medida (vector columna de dimensión n)
 * @param z Vector de medida observado (vector columna de dimensión m)
 * @param g Vector de medida estimado (predicción de la medida)
 * @param s Vector con las desviaciones típicas de la medida (dimensión m×1)
 * @param G Matriz de observación (dimensión m×n)
 * @param P Matriz de covarianza del estado estimado (dimensión n×n)
 * @param n Dimensión del vector de estado
 * @return Tupla con:
 *         - K: matriz de ganancia de Kalman (n×m)
 *         - x_new: nuevo vector de estado estimado (n×1)
 *         - P_new: nueva matriz de covarianza del estado (n×n)
 */
tuple<Matrix, Matrix, Matrix> MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n);

#endif


