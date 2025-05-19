#include "..\include\IERS.hpp"
#include <iostream>

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp) {
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    if (interp == 'l') {
        double mjd = floor(Mjd_UTC);
        int i = 0;
        for (int j = 1; j <= eop.n_column; ++j) {
            if (eop(4, j) == mjd) {
                i = j;
                break;
            }
        }
        if (i == 0 || i > eop.n_column) {
            std::cout << "IERS: No se encontró MJD " << mjd << " en eopdata\n";
            exit(EXIT_FAILURE);
        }
        if (i + 1 > eop.n_column) {
            std::cout << "IERS: No hay columna siguiente para interpolación (i=" << i << ", n_column=" << eop.n_column << ")\n";
            exit(EXIT_FAILURE);
        }
        Matrix preeop = eop.extract_column(i);
        Matrix nexteop = eop.extract_column(i + 1);
        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;
        x_pole  = preeop(5, 1) + (nexteop(5, 1) - preeop(5, 1)) * fixf;
        y_pole  = preeop(6, 1) + (nexteop(6, 1) - preeop(6, 1)) * fixf;
        UT1_UTC = preeop(7, 1) + (nexteop(7, 1) - preeop(7, 1)) * fixf;
        LOD     = preeop(8, 1) + (nexteop(8, 1) - preeop(8, 1)) * fixf;
        dpsi    = preeop(9, 1) + (nexteop(9, 1) - preeop(9, 1)) * fixf;
        deps    = preeop(10, 1) + (nexteop(10, 1) - preeop(10, 1)) * fixf;
        dx_pole = preeop(11, 1) + (nexteop(11, 1) - preeop(11, 1)) * fixf;
        dy_pole = preeop(12, 1) + (nexteop(12, 1) - preeop(12, 1)) * fixf;
        TAI_UTC = preeop(13, 1);

        x_pole  = x_pole / ARCS;
        y_pole  = y_pole / ARCS;
        dpsi    = dpsi / ARCS;
        deps    = deps / ARCS;
        dx_pole = dx_pole / ARCS;
        dy_pole = dy_pole / ARCS;
    } else {
        double mjd = floor(Mjd_UTC);
        int i = 0;
        for (int j = 1; j <= eop.n_column; ++j) {
            if (eop(4, j) == mjd) {
                i = j;
                break;
            }
        }
        if (i == 0 || i > eop.n_column) {
            std::cout << "IERS: No se encontró MJD " << mjd << " en eopdata\n";
            exit(EXIT_FAILURE);
        }
        Matrix eop_col = eop.extract_column(i);
        x_pole  = eop_col(5, 1) / ARCS;
        y_pole  = eop_col(6, 1) / ARCS;
        UT1_UTC = eop_col(7, 1);                
        LOD     = eop_col(8, 1);                 
        dpsi    = eop_col(9, 1) / ARCS;
        deps    = eop_col(10, 1) / ARCS;
        dx_pole = eop_col(11, 1) / ARCS; 
        dy_pole = eop_col(12, 1) / ARCS; 
        TAI_UTC = eop_col(13, 1);                
    }

    return make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}