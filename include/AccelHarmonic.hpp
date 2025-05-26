#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\global.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la aceleración gravitacional armónica en un punto dado, considerando los coeficientes de armónicos esféricos del geoide.
 * @param r Vector de posición en el sistema de coordenadas inercial [m] (3x1)
 * @param E Matriz de rotación del sistema inercial al sistema fijo a la Tierra (3x3)
 * @param n_max Grado máximo para la expansión en armónicos esféricos
 * @param m_max Orden máximo para la expansión en armónicos esféricos
 * @return Vector de aceleración en el sistema inercial [m/s²] (3x1)
 */
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif