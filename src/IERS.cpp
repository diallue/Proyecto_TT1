#include "..\include\IERS.hpp"
#include <iostream>

/**
 * Obtiene parámetros de orientación terrestre (EOP) del IERS para una fecha dada.
 * @param eop Matriz de datos EOP (Earth Orientation Parameters)
 * @param Mjd_UTC Fecha Juliana Modificada en UTC
 * @param interp Método de interpolación: 'l'=lineal, otro=sin interpolación
 * @return Tupla con 9 parámetros de orientación terrestre:
 *         [0] x_pole: Coordenada del polo [rad]
 *         [1] y_pole: Coordenada del polo [rad] 
 *         [2] UT1_UTC: Diferencia UT1-UTC [s]
 *         [3] LOD: Length of day [s]
 *         [4] dpsi: Corrección nutación en longitud [rad]
 *         [5] deps: Corrección nutación en oblicuidad [rad]
 *         [6] dx_pole: Tasa de cambio de x_pole [rad]
 *         [7] dy_pole: Tasa de cambio de y_pole [rad]
 *         [8] TAI_UTC: Diferencia TAI-UTC [s]
 * @throws Exit failure si no se encuentra la fecha en los datos EOP
 */
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