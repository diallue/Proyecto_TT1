#ifndef _GHAMATRIX_
#define _GHAMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\gast.hpp"
#include "..\include\R_z.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la matriz de rotación del ángulo horario de Greenwich (GHA), que transforma coordenadas del sistema inercial (ECI) al sistema fijo a la Tierra (ECEF).
 * @param Mjd_UT1 Fecha juliana modificada (UT1)
 * @return Matriz de rotación 3x3 que representa la transformación ECI -> ECEF
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif