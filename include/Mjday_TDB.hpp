#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

/**
 * @details Convierte el día juliano modificado en Tiempo Terrestre (Mjd_TT) a Tiempo Dinámico Baricéntrico (Mjd_TDB).
 * @param Mjd_TT Día juliano modificado en Tiempo Terrestre (MJD TT).
 * @return Día juliano modificado en Tiempo Dinámico Baricéntrico (MJD TDB).
 */
double Mjday_TDB(double Mjd_TT);

#endif