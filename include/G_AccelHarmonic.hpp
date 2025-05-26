#ifndef _GACCELHARMONIC_
#define _GACCELHARMONIC_

#include "..\include\matrix.hpp"
#include "..\include\AccelHarmonic.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula numéricamente el gradiente de la aceleración debida al campo gravitatorio armónico.
 * @param r Vector de posición en coordenadas geocéntricas (ECEF), dimensión 3x1
 * @param U Matriz de coeficientes armónicos normalizados del campo gravitatorio
 * @param n_max Grado máximo del desarrollo armónico
 * @param m_max Orden máximo del desarrollo armónico
 * @return Matriz 3x3 del gradiente de la aceleración (jacobiano)
 */
Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max);

#endif