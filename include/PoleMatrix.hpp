#ifndef _POLEMATRIX_
#define _POLEMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la matriz de transformación debida al movimiento del polo terrestre.
 * @param xp Coordenada del polo en dirección X (en radianes)
 * @param yp Coordenada del polo en dirección Y (en radianes)
 * @return Matriz de transformación 3x3 para corregir el movimiento del polo
 */
Matrix PoleMatrix(double xp, double yp);

#endif