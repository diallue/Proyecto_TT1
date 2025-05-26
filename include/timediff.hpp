#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include "..\include\matrix.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula diferencias entre escalas de tiempo astronómicas y de navegación.
 * @param UT1_UTC Diferencia UT1-UTC [s] (DUT1)
 * @param TAI_UTC Diferencia TAI-UTC [s] (número de segundos intercalares + 32.184)
 * @return Tupla con diferencias de tiempo:
 *         [0] UT1-TAI [s]
 *         [1] UTC-GPS [s]
 *         [2] UT1-GPS [s]
 *         [3] TT-UTC [s]
 *         [4] GPS-UTC [s]
 */
tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC);

#endif