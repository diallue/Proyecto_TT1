#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\R_z.hpp"
#include "..\include\R_y.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la matriz de precesión que transforma coordenadas del sistema de referencia ecuatorial medio de una época a otra.
 * @param Mjd_1 Fecha inicial en Tiempo Terrestre (Modified Julian Date)
 * @param Mjd_2 Fecha final en Tiempo Terrestre (Modified Julian Date)
 * @return Matriz de transformación 3x3 de precesión entre las dos fechas dadas
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif