#include "..\include\VarEqn.hpp"
#include <iostream>

/**
 * Calcula la ecuación variacional para la propagación de órbitas y matrices de transición de estado
 * 
 * @param x Tiempo desde la época de referencia [s]
 * @param yPhi Vector de estado extendido (42x1) que contiene:
 *             - Posición (3 elementos)
 *             - Velocidad (3 elementos)
 *             - Matriz de transición de estado (36 elementos)
 * @return Derivada del vector de estado extendido (42x1)
 */
Matrix VarEqn(double x, Matrix& yPhi) {
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = 
        IERS(eopdata, AuxParam.Mjd_UTC, 'l');
    
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
    
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC)/86400.0;

    Matrix P = PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x/86400.0);

    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);

    Matrix T = P * N;

    Matrix E = PoleMatrix(x_pole, y_pole);

    E = E * GHAMatrix(Mjd_UT1) * T;

    Matrix r(3, 1);
    r(1,1) = yPhi(1,1);
    r(2,1) = yPhi(2,1); 
    r(3,1) = yPhi(3,1); 
    
    Matrix v(3, 1);
    v(1,1) = yPhi(4,1); 
    v(2,1) = yPhi(5,1); 
    v(3,1) = yPhi(6,1);
    
    Matrix Phi(6, 6);
    for (int j = 1; j <= 6; ++j) {
        Matrix col(6, 1);
        for (int i = 1; i <= 6; ++i) {
            int index = 6*(j-1) + i + 6;
            col(i,1) = yPhi(index, 1);
        }
        Phi.assign_column(j, col);
    }

    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    Matrix G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    Matrix yPhip(42, 1);
    Matrix dfdy(6, 6);

    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            dfdy(i, j) = 0.0;                 
            dfdy(i+3, j) = G(i, j);           
            
            if (i == j) {
                dfdy(i, j+3) = 1.0;           
            } else {
                dfdy(i, j+3) = 0.0;          
            }
            
            dfdy(i+3, j+3) = 0.0;             
        }
    }

    Matrix Phip = dfdy * Phi;

    for (int i = 1; i <= 3; ++i) {
        yPhip(i, 1) = v(i, 1);               
        yPhip(i+3, 1) = a(i, 1);              
    }

    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 6; ++j) {
            yPhip(6*j + i, 1) = Phip(i, j);   
        }
    }

    return yPhip;
}