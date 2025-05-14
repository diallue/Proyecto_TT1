#include "..\include\PrecMatrix.hpp"

Matrix PrecMatrix(double Mjd_1, double Mjd_2) {
    double T = (Mjd_1 - MJD_J2000) / 36525.0;
    double dT = (Mjd_2 - Mjd_1) / 36525.0;

    double zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) + ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / ARCS;
    double z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / ARCS;
    double theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) - ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / ARCS;

    Matrix PrecMat = R_z(-z) * R_y(theta) * R_z(-zeta);

    return PrecMat;
}