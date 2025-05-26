#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>
#include <vector>

/**
 * @details Obtiene parámetros de orientación terrestre (EOP) del IERS para una fecha dada.
 * @param eop Matriz de datos EOP (Earth Orientation Parameters)
 * @param Mjd_UTC Fecha Juliana Modificada en UTC
 * @param interp Método de interpolación: 'l'=lineal, otro=sin interpolación
 * @return Tupla con 9 parámetros de orientación terrestre:
 *         [0] x_pole: Coordenada del polo [rad]
 *         [1] y_pole: Coordenada del polo [rad] 
 *         [2] UT1_UTC: Diferencia UT1-UTC [s]
 *         [3] LOD: Length of day [s]
 *         [4] dpsi: Corrección nutación en longitud [rad]
 *         [5] deps: Corrección nutación en oblicuidad [rad]
 *         [6] dx_pole: Tasa de cambio de x_pole [rad]
 *         [7] dy_pole: Tasa de cambio de y_pole [rad]
 *         [8] TAI_UTC: Diferencia TAI-UTC [s]
 * @throws Exit failure si no se encuentra la fecha en los datos EOP
 */
tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp = 'n');

#endif