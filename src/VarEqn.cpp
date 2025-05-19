#include "..\include\VarEqn.hpp"

Matrix VarEqn(double x, Matrix& yPhi) {
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = 
        IERS(eopdata, AuxParam.Mjd_UTC, 'l');
    
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
    
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC)/86400.0;

    Matrix P = PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x/86400.0);
    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix r = yPhi.extract_vector(1, 3, false);
    Matrix v = yPhi.extract_vector(4, 6, false);
    
    Matrix Phi(6, 6);
    for (int j = 1; j <= 6; ++j) {
        Phi.assign_column(j, yPhi.extract_vector(6*j+1, 6*j+6, false));
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