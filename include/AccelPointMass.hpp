#ifndef _ACCEL_POINT_MASS_
#define _ACCEL_POINT_MASS_

#include "..\include\matrix.hpp"
#include <cmath>

/**
 * @details Calcula la aceleración gravitacional en un punto debido a una masa puntual.
 * @param r Vector de posición del punto (3x1, metros).
 * @param s Vector de posición de la masa (3x1, metros).
 * @param GM Parámetro gravitacional (m³/s²).
 * @return Vector de aceleración (3x1, m/s²).
 * @throw Error si r o s no son 3x1 o si normas son cero.
 */
Matrix AccelPointMass(Matrix& r, Matrix& s, double GM);

#endif