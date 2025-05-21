#include "..\include\Mjday.hpp"

/**
 * Calcula el día juliano modificado (MJD) a partir de una fecha y hora.
 * @param yr Año (entero).
 * @param mon Mes (1 a 12).
 * @param day Día del mes (1 a 31).
 * @param hr Hora (0 a 23, por defecto 0).
 * @param min Minutos (0 a 59, por defecto 0).
 * @param sec Segundos (0 a 59, por defecto 0).
 * @return Día juliano modificado (MJD).
 */
double Mjday(double yr, double mon, double day, double hr = 0, double min = 0, double sec = 0) {
    double jd = 367.0 * yr
        - floor( (7 * (yr + floor( (mon + 9) / 12.0 ) ) ) * 0.25 )
        + floor( 275 * mon / 9.0 )
        + day + 1721013.5
        + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

    return jd - 2400000.5;
}