#ifndef _LTC_
#define _LTC_

#include "..\include\matrix.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la matriz de transformación Local Tangent Coordinate (LTC) en función de la latitud y longitud geodésica del observador sobre la Tierra.
 * @param lon Longitud geodésica en radianes
 * @param lat Latitud geodésica en radianes
 * @return Matriz de transformación 3x3 del sistema ECEF al sistema local ENZ
 */
Matrix LTC(double lon, double lat);

#endif