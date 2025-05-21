#include "..\include\EqnEquinox.hpp"

/**
 * Calcula la ecuación de los equinoccios (diferencia entre tiempo sidéreo aparente y medio)
 * 
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre (TT)
 * @return Ecuación de los equinoccios en radianes [rad]
 */
double EqnEquinox(double Mjd_TT) {
    auto [dpsi, deps] = NutAngles(Mjd_TT);

    double EqE = dpsi * cos(MeanObliquity(Mjd_TT));

    return EqE;
}