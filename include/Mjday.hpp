#ifndef _MJDAY_
#define _MJDAY_

#include "..\include\matrix.hpp"
#include <cmath>

/**
 * @details Calcula el día juliano modificado (MJD) a partir de una fecha y hora.
 * @param yr Año (entero).
 * @param mon Mes (1 a 12).
 * @param day Día del mes (1 a 31).
 * @param hr Hora (0 a 23, por defecto 0).
 * @param min Minutos (0 a 59, por defecto 0).
 * @param sec Segundos (0 a 59, por defecto 0).
 * @return Día juliano modificado (MJD).
 */
double Mjday(double yr, double mon, double day, double hr, double min, double sec);

#endif