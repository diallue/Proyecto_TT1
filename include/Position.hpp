#ifndef _POSITION_
#define _POSITION_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

/**
 * @details Calcula el vector de posición en coordenadas cartesianas a partir de coordenadas geodésicas.
 * @param lon Longitud geodésica (en radianes).
 * @param lat Latitud geodésica (en radianes).
 * @param h Altura sobre el elipsoide (en metros).
 * @return Vector de posición cartesiana (3x1, en metros).
 */
Matrix Position(double lon, double lat, double h);

#endif