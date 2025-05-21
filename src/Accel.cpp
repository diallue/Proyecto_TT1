#include "..\include\Accel.hpp"

/**
 * Calcula la aceleración total que actúa sobre un satélite considerando múltiples perturbaciones
 * 
 * @param x Tiempo desde la época de referencia [s]
 * @param Y Vector de estado (6x1) [posición; velocidad] [m, m/s]
 * @return Derivada del vector de estado (6x1) [velocidad; aceleración] [m/s, m/s²]
 */
Matrix Accel(double x, Matrix Y) {
    if (Y.n_row < Y.n_column) {
        Y = transpose(Y);
    }
    if (Y.n_row != 6 || Y.n_column != 1) {
        std::cout << "Accel: Y debe ser 6x1, dims = " << Y.n_row << "x" << Y.n_column << std::endl;
        exit(EXIT_FAILURE);
    }

    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = 
        IERS(eopdata, AuxParam.Mjd_UTC + x / 86400.0, 'l');

    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = 
        timediff(UT1_UTC, TAI_UTC);

    auto Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400.0 + UT1_UTC / 86400.0;
    auto Mjd_TT = AuxParam.Mjd_UTC + x / 86400.0 + TT_UTC / 86400.0;

    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    auto MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, 
           r_Neptune, r_Pluto, r_Moon, r_Sun;
    std::tie(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, 
             r_Neptune, r_Pluto, r_Moon, r_Sun) = JPL_Eph_DE430(MJD_TDB);

    Matrix r(3, 1);
    for (int i = 1; i <= 3; ++i) {
        r(i, 1) = Y(i, 1);
    }

    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    if (AuxParam.sun) {
        a = a + AccelPointMass(r, r_Sun, GM_SUN);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r, r_Moon, GM_MOON);
    }

    if (AuxParam.planets) {
        a = a + AccelPointMass(r, r_Mercury, GM_MERCURY);
        a = a + AccelPointMass(r, r_Venus, GM_VENUS);
        a = a + AccelPointMass(r, r_Mars, GM_MARS);
        a = a + AccelPointMass(r, r_Jupiter, GM_JUPITER);
        a = a + AccelPointMass(r, r_Saturn, GM_SATURN);
        a = a + AccelPointMass(r, r_Uranus, GM_URANUS);
        a = a + AccelPointMass(r, r_Neptune, GM_NEPTUNE);
        a = a + AccelPointMass(r, r_Pluto, GM_PLUTO);
    }

    Matrix dY = zeros(6, 1);
    for (int i = 1; i <= 3; ++i) {
        dY(i, 1) = Y(i + 3, 1);
        dY(i + 3, 1) = a(i, 1);
    }

    return dY;
}