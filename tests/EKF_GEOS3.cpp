#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"
#include <tuple>
#include <cmath>
#include <iomanip>

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
    Matrix Rs = pos;

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
    for (int i = 1; i <= 3; i++) {
        Y0_apr(i, 1) = r2(i, 1);
        Y0_apr(i + 3, 1) = v2(i, 1);
    }

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = obs(9, 1);
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    double n_eqn = 6;

    Matrix Y = DEInteg(Accel, 0, -(Mjd_UTC - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);

    Matrix P = zeros(6, 6);
    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT = LTC(lon, lat);

    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    double t = 0;

    for (int i = 1; i <= 46; i++) {
        double t_old = t;
        Matrix Y_old = Y;

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        std::tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(eopdata, Mjd_UTC, 'l');

        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        std::tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; j++) {
                yPhi(6 * j + ii, 1) = (ii == j) ? 1.0 : 0.0;
            }
        }

        yPhi = DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);

        for (int j = 1; j <= 6; j++) {
            for (int k = 1; k <= 6; k++) {
                Phi(k, j) = yPhi(6 * j + k, 1);
            }
        }

        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);

        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        Matrix r(3, 1);
        for (int j = 1; j <= 3; j++) {
            r(j, 1) = Y(j, 1);
        }

        Matrix Ur = U * r;

        Matrix s = LT * (Ur - Rs);

        Matrix Qdt = zeros(6, 6);
        TimeUpdate(P, Phi, Qdt);

        double Azim, Elev;
        Matrix dAds(1, 3), dEds(1, 3);
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);
        Matrix dAdY(1, 6);
        Matrix LT_U = LT * U;
        for (int j = 1; j <= 3; j++) {
            dAdY(1, j) = dAds(1, j) * LT_U(1, j);
            dAdY(1, j + 3) = 0.0;
        }

        Matrix K, Y_new, P_new;
        Matrix z(1, 1), g(1, 1), s_sigma(1, 1);
        z(1, 1) = obs(i, 2);
        g(1, 1) = Azim;
        s_sigma(1, 1) = sigma_az;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dAdY, P, 6);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; j++) {
            r(j, 1) = Y(j, 1);
        }
        s = LT * (U * r - Rs);
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);
        Matrix dEdY(1, 6);
        for (int j = 1; j <= 3; j++) {
            dEdY(1, j) = dEds(1, j) * LT_U(1, j);
            dEdY(1, j + 3) = 0.0;
        }

        z(1, 1) = obs(i, 3);
        g(1, 1) = Elev;
        s_sigma(1, 1) = sigma_el;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dEdY, P, 6);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; j++) {
            r(j, 1) = Y(j, 1);
        }
        s = LT * (U * r - Rs);
        double Dist = s.norm();
        Matrix dDds = transpose(s / Dist);
        Matrix dDdY(1, 6);
        for (int j = 1; j <= 3; j++) {
            dDdY(1, j) = dDds(1, j) * LT_U(1, j);
            dDdY(1, j + 3) = 0.0;
        }

        z(1, 1) = obs(i, 4);
        g(1, 1) = Dist;
        s_sigma(1, 1) = sigma_range;
        std::tie(K, Y_new, P_new) = MeasUpdate(Y, z, g, s_sigma, dDdY, P, 6);
        Y = Y_new;
        P = P_new;
    }

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    std::tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(eopdata, obs(46, 1), 'l');

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    std::tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg(Accel, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);

    Matrix Y_true(6, 1);
    Y_true(1, 1) = 5753.173e3;
    Y_true(2, 1) = 2673.361e3;
    Y_true(3, 1) = 3440.304e3;
    Y_true(4, 1) = 4.324207e3;
    Y_true(5, 1) = -1.924299e3;
    Y_true(6, 1) = -5.728216e3;

    std::cout << "\nError de Estimación de Posición\n";
    std::cout << "dX " << std::setw(10) << std::fixed << std::setprecision(1) << Y0(1, 1) - Y_true(1, 1) << " [m]\n";
    std::cout << "dY " << std::setw(10) << std::fixed << std::setprecision(1) << Y0(2, 1) - Y_true(2, 1) << " [m]\n";
    std::cout << "dZ " << std::setw(10) << std::fixed << std::setprecision(1) << Y0(3, 1) - Y_true(3, 1) << " [m]\n";
    std::cout << "\nError de Estimación de Velocidad\n";
    std::cout << "dVx " << std::setw(8) << std::fixed << std::setprecision(1) << Y0(4, 1) - Y_true(4, 1) << " [m/s]\n";
    std::cout << "dVy " << std::setw(8) << std::fixed << std::setprecision(1) << Y0(5, 1) - Y_true(5, 1) << " [m/s]\n";
    std::cout << "dVz " << std::setw(8) << std::fixed << std::setprecision(1) << Y0(6, 1) - Y_true(6, 1) << " [m/s]\n";

    return 0;
}