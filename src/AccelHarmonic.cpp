#include "..\include\AccelHarmonic.hpp"

static Matrix Cnm;
static Matrix Snm;

Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {
    const double r_ref = 6378.1363e3;     // Radio de referencia (m)
    const double gm = 398600.4415e9;      // Parámetro gravitacional (m^3/s^2)

    // Vector posición en el sistema del cuerpo (body-fixed)
    Matrix r_bf = E * r;

    double x = r_bf(1,1);
    double y = r_bf(2,1);
    double z = r_bf(3,1);

    double d = sqrt(x*x + y*y + z*z);        // Norma del vector
    double latgc = asin(z / d);              // Latitud geocéntrica
    double lon = atan2(y, x);                // Longitud

    Matrix pnm(n_max+2, m_max+2);
    Matrix dpnm(n_max+2, m_max+2);
    Legendre(n_max, m_max, latgc, pnm, dpnm);

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;

    for (int n = 0; n <= n_max; ++n) {
        double b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * pow(r_ref / d, n);
        double b3 = (gm / d) * pow(r_ref / d, n);

        double q1 = 0.0;
        double q2 = 0.0;
        double q3 = 0.0;

        for (int m = 0; m <= m_max; ++m) {
            double cosml = cos(m * lon);
            double sinml = sin(m * lon);

            double C = Cnm(n+1, m+1);
            double S = Snm(n+1, m+1);

            q1 += pnm(n+1, m+1) * (C * cosml + S * sinml);
            q2 += dpnm(n+1, m+1) * (C * cosml + S * sinml);
            q3 += m * pnm(n+1, m+1) * (S * cosml - C * sinml);
        }

        dUdr += b1 * q1;
        dUdlatgc += b2 * q2;
        dUdlon += b3 * q3;
    }

    double r2xy = x*x + y*y;
    double rho = sqrt(r2xy);

    double ax = (1.0 / d * dUdr - z / (d * d * rho) * dUdlatgc) * x - (1.0 / r2xy) * dUdlon * y;
    double ay = (1.0 / d * dUdr - z / (d * d * rho) * dUdlatgc) * y + (1.0 / r2xy) * dUdlon * x;
    double az = (1.0 / d * dUdr) * z + rho / (d * d) * dUdlatgc;

    Matrix a_bf(3, 1);
    a_bf(1,1) = ax;
    a_bf(2,1) = ay;
    a_bf(3,1) = az;

    // Volvemos al sistema inercial aplicando E^T porque r_bf = E * r
    Matrix a = transpose(E) * a_bf;

    return a;
}