#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula la matriz de nutación que transforma coordenadas del sistema medio al sistema verdadero.
 * @param Mjd_TT Fecha en Tiempo Terrestre (Modified Julian Date)
 * @return Matriz de transformación 3x3 que representa la nutación
 */
Matrix NutMatrix(double Mjd_TT);

#endif