#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Position.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\LTC.hpp"
#include "..\include\gmst.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"
#include <tuple>
#include <cmath>
#include <iomanip>
#include <fstream>

int main() {
	eop19620101(24324);
    GGM03S(181);
    DE430Coeff(2285, 1020);
    AuxParamLoad();
	GEOS3(46);

    double sigma_range = 92.5;         
    double sigma_az = 0.0224 * RAD;    
    double sigma_el = 0.0139 * RAD;    

    double lat = RAD * 21.5748;        
    double lon = RAD * (-158.2706);    
    double alt = 300.20;               

    Matrix pos = Position(lon, lat, alt);
    if (pos.n_row != 3 || pos.n_column != 1) {
        std::cerr << "Error: pos tiene dimensiones incorrectas (" << pos.n_row << "x" << pos.n_column << ")\n";
        return 1;
    }
    Matrix Rs = pos;
    if (Rs.n_row != 3 || Rs.n_column != 1) {
        std::cerr << "Error: Rs tiene dimensiones incorrectas (" << Rs.n_row << "x" << Rs.n_column << ")\n";
        return 1;
    }

    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

    Matrix r2(3, 1), v2(3, 1);
    r2(1, 1) = 6221397.62857869;
    r2(2, 1) = 2867713.77965738;
    r2(3, 1) = 3006155.98509949;
    v2(1, 1) = 4645.04725161806;
    v2(2, 1) = -2752.21591588204;
    v2(3, 1) = -7507.99940987031;

    Matrix Y0_apr(6, 1);
    for (int i = 1; i <= 3; ++i) {
        Y0_apr(i, 1) = r2(i, 1);
        Y0_apr(i + 3, 1) = v2(i, 1);
    }

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = obs(9, 1);
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = true;
    AuxParam.moon = true;
    AuxParam.planets = true;
    int n_eqn = 6;

    Matrix Y = DEInteg(Accel, 0, -(Mjd_UTC - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
    if (Y.n_row != 6 || Y.n_column != 1) {
        std::cerr << "Error: Y tiene dimensiones incorrectas después de la primera integración (" << Y.n_row << "x" << Y.n_column << ")\n";
        return 1;
    }

    Matrix P = zeros(6, 6);
    for (int i = 1; i <= 3; ++i) P(i, i) = 1e8;
    for (int i = 4; i <= 6; ++i) P(i, i) = 1e3; 

    Matrix LT = LTC(lon, lat);
    if (LT.n_row != 3 || LT.n_column != 3) {
        std::cerr << "Error: LT tiene dimensiones incorrectas (" << LT.n_row << "x" << LT.n_column << ")\n";
        return 1;
    }

    Matrix yPhi = zeros(42, 1);
    Matrix Phi = zeros(6, 6);

    double t = 0;

    for (int i = 1; i <= 46; ++i) {
        double t_old = t;
        Matrix Y_old = Y;

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        std::tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) =
            IERS(eopdata, Mjd_UTC, 'l');
        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        std::tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) =
            timediff(UT1_UTC, TAI_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; ++j) {
                yPhi(6 * j + ii, 1) = (ii == j) ? 1.0 : 0.0;
            }
        }

        yPhi = DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
        if (yPhi.n_row != 42 || yPhi.n_column != 1) {
            std::cerr << "Error: yPhi tiene dimensiones incorrectas en la iteración " << i << " (" << yPhi.n_row << "x" << yPhi.n_column << ")\n";
            return 1;
        }

        for (int j = 1; j <= 6; ++j) {
            Matrix col(6, 1);
            for (int k = 1; k <= 6; ++k) {
                col(k, 1) = yPhi(6 * j + k, 1);
            }
            Phi.assign_column(j, col);
        }

        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
        if (Y.n_row != 6 || Y.n_column != 1) {
            std::cerr << "Error: Y tiene dimensiones incorrectas en la iteración " << i << " (" << Y.n_row << "x" << Y.n_column << ")\n";
            return 1;
        }

        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        if (U.n_row != 3 || U.n_column != 3) {
            std::cerr << "Error: U tiene dimensiones incorrectas (" << U.n_row << "x" << U.n_column << ")\n";
            return 1;
        }
        Matrix r(3, 1);
        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Matrix Ur = U * r;
        if (Ur.n_row != 3 || Ur.n_column != 1) {
            std::cerr << "Error: U*r tiene dimensiones incorrectas (" << Ur.n_row << "x" << Ur.n_column << ")\n";
            return 1;
        }
        if (Ur.n_row != Rs.n_row || Ur.n_column != Rs.n_column) {
            std::cerr << "Error: Dimensiones incompatibles para U*r - Rs en la iteración " << i
                      << " (U*r: " << Ur.n_row << "x" << Ur.n_column << ", Rs: " << Rs.n_row << "x" << Rs.n_column << ")\n";
            return 1;
        }
        Matrix s = LT * (Ur - Rs);
        if (s.n_row != 3 || s.n_column != 1) {
            std::cerr << "Error: s tiene dimensiones incorrectas (" << s.n_row << "x" << s.n_column << ")\n";
            return 1;
        }

        P = TimeUpdate(P, Phi);

        double Azim, Elev;
        Matrix dAds(1, 3), dEds(1, 3);
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);
        Matrix zeros_1x3 = zeros(1, 3);
        Matrix dAdY = dAds * LT * U;
        dAdY = dAdY.union_vector(zeros_1x3, true);
        Matrix K, Y_new, P_new;
        Matrix z(1, 1), g(1, 1), s_sigma(1, 1);
        z(1, 1) = obs(i, 2);
        g(1, 1) = Azim;
        s_sigma(1, 1) = sigma_az;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dAdY, P, 6);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);
        Matrix dEdY = dEds * LT * U;
        dEdY = dEdY.union_vector(zeros_1x3, true); 
        z(1, 1) = obs(i, 3);
        g(1, 1) = Elev;
        s_sigma(1, 1) = sigma_el;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dEdY, P, 6);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
        double Dist = s.norm();
        Matrix dDds = transpose(s / Dist);
        Matrix dDdY = dDds * LT * U;
        dDdY = dDdY.union_vector(zeros_1x3, true); 
        z(1, 1) = obs(i, 4);
        g(1, 1) = Dist;
        s_sigma(1, 1) = sigma_range;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dDdY, P, 6);
        Y = Y_new;
        P = P_new;
    }

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    std::tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) =
        IERS(eopdata, obs(46, 1), 'l');
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    std::tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) =
        timediff(UT1_UTC, TAI_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg(Accel, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);

    double aux[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true(6, 1);
    for (int i = 0; i < 6; ++i) {
        Y_true(i + 1, 1) = aux[i];
    }

    std::cout << "\nError de Estimación de Posición\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "dX " << std::setw(10) << Y0(1, 1) - Y_true(1, 1) << " [m]\n";
    std::cout << "dY " << std::setw(10) << Y0(2, 1) - Y_true(2, 1) << " [m]\n";
    std::cout << "dZ " << std::setw(10) << Y0(3, 1) - Y_true(3, 1) << " [m]\n";
    std::cout << "\nError de Estimación de Velocidad\n";
    std::cout << "dVx " << std::setw(8) << Y0(4, 1) - Y_true(4, 1) << " [m/s]\n";
    std::cout << "dVy " << std::setw(8) << Y0(5, 1) - Y_true(5, 1) << " [m/s]\n";
    std::cout << "dVz " << std::setw(8) << Y0(6, 1) - Y_true(6, 1) << " [m/s]\n";
	
    return 0;
}