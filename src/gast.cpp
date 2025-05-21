#include "..\include\gast.hpp"

/**
 * Calcula el Tiempo Sideral Aparente de Greenwich (GAST) en radianes para una fecha dada.
 * 
 * @param Mjd_UT1 Fecha en Tiempo Universal UT1 (Modified Julian Date)
 * @return Tiempo Sideral Aparente de Greenwich en radianes (rango [0, 2π])
 */
double gast(double Mjd_UT1) {
    return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), PI2);
}