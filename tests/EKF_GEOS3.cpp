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
    GGM03S(16471);
    DE430Coeff(2285, 1020);
    AuxParamLoad();
    GEOS3(46);
	
	std::cout << "Validating dimensions...\n";
    if (eopdata.n_row != 13 || eopdata.n_column != 24324) {
        std::cerr << "Invalid eopdata dimensions: " << eopdata.n_row << " x " << eopdata.n_column << "\n";
        return 1;
    }
    if (Cnm.n_row != 16471 || Cnm.n_column != 16471 || Snm.n_row != 16471 || Snm.n_column != 16471) {
        std::cerr << "Invalid Cnm/Snm dimensions\n";
        return 1;
    }
    if (PC.n_row != 2285 || PC.n_column != 1020) {
        std::cerr << "Invalid PC dimensions: " << PC.n_row << " x " << PC.n_column << "\n";
        return 1;
    }
    if (obs.n_row != 46 || obs.n_column != 4) {
        std::cerr << "Invalid obs dimensions: " << obs.n_row << " x " << obs.n_column << "\n";
        return 1;
    }

    double sigma_range = 92.5;        
    double sigma_az = 0.0224 * RAD;    
    double sigma_el = 0.0139 * RAD;   

    double lat = RAD * 21.5748;       
    double lon = RAD * (-158.2706);    
    double alt = 300.20;              

	std::cout << "Computing station position...\n";
    Matrix Rs = Position(lon, lat, alt);
	std::cout << "Rs dimensions: " << Rs.n_row << " x " << Rs.n_column << "\n";

    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

    Matrix r2(3, 1), v2(3, 1);
    r2(1, 1) = 6221397.62857869;
    r2(2, 1) = 2867713.77965738;
    r2(3, 1) = 3006155.98509949;
    v2(1, 1) = 3837.62169313926;
    v2(2, 1) = -1382.03320148307;
    v2(3, 1) = -5978.16810402366;
	std::cout << "r2 dimensions: " << r2.n_row << " x " << r2.n_column << "\n";
    std::cout << "v2 dimensions: " << v2.n_row << " x " << v2.n_column << "\n";
    //anglesg(obs(1, 2), obs(9, 2), obs(18, 2), obs(1, 3), obs(9, 3), obs(18, 3), Mjd1, Mjd2, Mjd3, Rs, Rs, Rs, r2, v2);

    Matrix Y0_apr(6, 1);
    for (int i = 1; i <= 3; ++i) {
        Y0_apr(i, 1) = r2(i, 1);
        Y0_apr(i + 3, 1) = v2(i, 1);
    }
	std::cout << "Y0_apr dimensions: " << Y0_apr.n_row << " x " << Y0_apr.n_column << "\n";

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);
    double Mjd_UTC = obs(9, 1);
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

	std::cout << "Integrating initial orbit...\n";
    Matrix Y = DEInteg(Accel, 0, -(obs(9, 1) - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
	std::cout << "Y dimensions: " << Y.n_row << " x " << Y.n_column << "\n";

    Matrix P(6, 6);
    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 6; ++j) {
            P(i, j) = 0.0;
        }
    }
    for (int i = 1; i <= 3; ++i) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; ++i) {
        P(i, i) = 1e3;
    }
	std::cout << "P dimensions: " << P.n_row << " x " << P.n_column << "\n";

    Matrix LT = LTC(lon, lat);
    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);
    Matrix Qdt(6, 6);
    double t = 0;
	std::cout << "LT dimensions: " << LT.n_row << " x " << LT.n_column << "\n";
    std::cout << "yPhi dimensions: " << yPhi.n_row << " x " << yPhi.n_column << "\n";
    std::cout << "Phi dimensions: " << Phi.n_row << " x " << Phi.n_column << "\n";
    std::cout << "Qdt dimensions: " << Qdt.n_row << " x " << Qdt.n_column << "\n";

    for (int i = 1; i <= 46; ++i) {
		std::cout << "\nProcessing measurement " << i << "...\n";
        double t_old = t;
        Matrix Y_old = Y;

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;
		std::cout << "Mjd_UTC: " << Mjd_UTC << ", t: " << t << "\n";

		std::cout << "Calling IERS...\n";
        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        auto iers_result = IERS(eopdata, Mjd_UTC, 'l');
        x_pole = std::get<0>(iers_result);
        y_pole = std::get<1>(iers_result);
        UT1_UTC = std::get<2>(iers_result);
        LOD = std::get<3>(iers_result);
        dpsi = std::get<4>(iers_result);
        deps = std::get<5>(iers_result);
        dx_pole = std::get<6>(iers_result);
        dy_pole = std::get<7>(iers_result);
        TAI_UTC = std::get<8>(iers_result);

		std::cout << "Calling timediff...\n";
        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        auto time_result = timediff(UT1_UTC, TAI_UTC);
        UT1_TAI = std::get<0>(time_result);
        UTC_GPS = std::get<1>(time_result);
        UT1_GPS = std::get<2>(time_result);
        TT_UTC = std::get<3>(time_result);
        GPS_UTC = std::get<4>(time_result);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

		std::cout << "Initializing yPhi...\n";
        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; ++j) {
                yPhi(6 * (j - 1) + ii, 1) = (ii == j) ? 1.0 : 0.0;
            }
        }

		std::cout << "Integrating variational equations...\n";
        yPhi = DEInteg([](double t, Matrix z) { return VarEqn(t, z); }, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
		std::cout << "yPhi dimensions: " << yPhi.n_row << " x " << yPhi.n_column << "\n";

		std::cout << "Extracting Phi...\n";
        for (int j = 1; j <= 6; ++j) {
            for (int k = 1; k <= 6; ++k) {
                Phi(k, j) = yPhi(6 * (j - 1) + k, 1);
            }
        }

		std::cout << "Integrating state...\n";
        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
		std::cout << "Y dimensions: " << Y.n_row << " x " << Y.n_column << "\n";

		std::cout << "Computing topocentric coordinates...\n";
        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        Matrix r(3, 1);
        for (int j = 1; j <= 3; ++j) {
            r(j, 1) = Y(j, 1);
        }
        Matrix s = LT * (U * r - Rs);
		std::cout << "U dimensions: " << U.n_row << " x " << U.n_column << "\n";
        std::cout << "r dimensions: " << r.n_row << " x " << r.n_column << "\n";
        std::cout << "s dimensions: " << s.n_row << " x " << s.n_column << "\n";

        std::cout << "Performing time update...\n";
        TimeUpdate(P, Phi, Qdt);

		std::cout << "Computing azimuth and elevation...\n";
        double Azim, Elev;
        Matrix dAds(1, 3), dEds(1, 3);
        auto azel_result = AzElPa(s);
        Azim = std::get<0>(azel_result);
        Elev = std::get<1>(azel_result);
        dAds = std::get<2>(azel_result);
        dEds = std::get<3>(azel_result);
		std::cout << "dAds dimensions: " << dAds.n_row << " x " << dAds.n_column << "\n";
        std::cout << "dEds dimensions: " << dEds.n_row << " x " << dEds.n_column << "\n";

        std::cout << "Computing azimuth partials...\n";
        Matrix dAdY(1, 6);
        Matrix temp = dAds * LT * U;
        for (int j = 1; j <= 3; ++j) {
            dAdY(1, j) = temp(1, j);
            dAdY(1, j + 3) = 0.0;
        }
		std::cout << "dAdY dimensions: " << dAdY.n_row << " x " << dAdY.n_column << "\n";

        std::cout << "Performing azimuth measurement update...\n";
        Matrix K;
		Matrix z(1, 1); z(1, 1) = obs(i, 2);
		Matrix g(1, 1); g(1, 1) = Azim;
		Matrix s_sigma(1, 1); s_sigma(1, 1) = sigma_az;
		std::cout << "z dimensions: " << z.n_row << " x " << z.n_column << "\n";
		std::cout << "g dimensions: " << g.n_row << " x " << g.n_column << "\n";
		std::cout << "s_sigma dimensions: " << s_sigma.n_row << " x " << s_sigma.n_column << "\n";
		auto meas_result = MeasUpdate(Y, z, g, s_sigma, dAdY, P, 6);
        K = std::get<0>(meas_result);
        Y = std::get<1>(meas_result);
        P = std::get<2>(meas_result);

		std::cout << "Updating topocentric coordinates for elevation...\n";
        for (int j = 1; j <= 3; ++j) {
            r(j, 1) = Y(j, 1);
        }
        s = LT * (U * r - Rs);
        azel_result = AzElPa(s);
        Azim = std::get<0>(azel_result);
        Elev = std::get<1>(azel_result);
        dAds = std::get<2>(azel_result);
        dEds = std::get<3>(azel_result);

		std::cout << "Computing elevation partials...\n";
        Matrix dEdY(1, 6);
        temp = dEds * LT * U;
        for (int j = 1; j <= 3; ++j) {
            dEdY(1, j) = temp(1, j);
            dEdY(1, j + 3) = 0.0;
        }
		std::cout << "dEdY dimensions: " << dEdY.n_row << " x " << dEdY.n_column << "\n";

        std::cout << "Performing elevation measurement update...\n";
        z(1, 1) = obs(i, 3);
		g(1, 1) = Elev;
		s_sigma(1, 1) = sigma_el;
		std::cout << "z dimensions: " << z.n_row << " x " << z.n_column << "\n";
		std::cout << "g dimensions: " << g.n_row << " x " << g.n_column << "\n";
		std::cout << "s_sigma dimensions: " << s_sigma.n_row << " x " << s_sigma.n_column << "\n";
		meas_result = MeasUpdate(Y, z, g, s_sigma, dEdY, P, 6);
        K = std::get<0>(meas_result);
        Y = std::get<1>(meas_result);
        P = std::get<2>(meas_result);

		std::cout << "Computing range...\n";
        for (int j = 1; j <= 3; ++j) {
            r(j, 1) = Y(j, 1);
        }
        s = LT * (U * r - Rs);
        double Dist = s.norm();
        Matrix dDds = transpose(s / Dist);
        Matrix dDdY(1, 6);
        temp = dDds * LT * U;
        for (int j = 1; j <= 3; ++j) {
            dDdY(1, j) = temp(1, j);
            dDdY(1, j + 3) = 0.0;
        }
		std::cout << "dDds dimensions: " << dDds.n_row << " x " << dDds.n_column << "\n";
        std::cout << "dDdY dimensions: " << dDdY.n_row << " x " << dDdY.n_column << "\n";

        std::cout << "Performing range measurement update...\n";
        z(1, 1) = obs(i, 4);
		g(1, 1) = Dist;
		s_sigma(1, 1) = sigma_range;
		std::cout << "z dimensions: " << z.n_row << " x " << z.n_column << "\n";
		std::cout << "g dimensions: " << g.n_row << " x " << g.n_column << "\n";
		std::cout << "s_sigma dimensions: " << s_sigma.n_row << " x " << s_sigma.n_column << "\n";
		meas_result = MeasUpdate(Y, z, g, s_sigma, dDdY, P, 6);
        K = std::get<0>(meas_result);
        Y = std::get<1>(meas_result);
        P = std::get<2>(meas_result);
    }

	std::cout << "Final IERS call...\n";
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    auto iers_result = IERS(eopdata, obs(46, 1), 'l');
    x_pole = std::get<0>(iers_result);
    y_pole = std::get<1>(iers_result);
    UT1_UTC = std::get<2>(iers_result);
    LOD = std::get<3>(iers_result);
    dpsi = std::get<4>(iers_result);
    deps = std::get<5>(iers_result);
    dx_pole = std::get<6>(iers_result);
    dy_pole = std::get<7>(iers_result);
    TAI_UTC = std::get<8>(iers_result);

	std::cout << "Final timediff call...\n";
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    auto time_result = timediff(UT1_UTC, TAI_UTC);
    UT1_TAI = std::get<0>(time_result);
    UTC_GPS = std::get<1>(time_result);
    UT1_GPS = std::get<2>(time_result);
    TT_UTC = std::get<3>(time_result);
    GPS_UTC = std::get<4>(time_result);

    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

	std::cout << "Final integration...\n";
    Matrix Y0 = DEInteg(Accel, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);
	std::cout << "Y0 dimensions: " << Y0.n_row << " x " << Y0.n_column << "\n";

    Matrix Y_true(6, 1);
    Y_true(1, 1) = 5753.173e3;
    Y_true(2, 1) = 2673.361e3;
    Y_true(3, 1) = 3440.304e3;
    Y_true(4, 1) = 4.324207e3;
    Y_true(5, 1) = -1.924299e3;
    Y_true(6, 1) = -5.728216e3;

    std::cout << "\nError of Position Estimation\n";
    std::cout << "dX" << std::setw(10) << std::fixed << std::setprecision(1) << Y0(1, 1) - Y_true(1, 1) << " [m]\n";
    std::cout << "dY" << std::setw(10) << std::fixed << std::setprecision(1) << Y0(2, 1) - Y_true(2, 1) << " [m]\n";
    std::cout << "dZ" << std::setw(10) << std::fixed << std::setprecision(1) << Y0(3, 1) - Y_true(3, 1) << " [m]\n";
    std::cout << "\nError of Velocity Estimation\n";
    std::cout << "dVx" << std::setw(8) << std::fixed << std::setprecision(1) << Y0(4, 1) - Y_true(4, 1) << " [m/s]\n";
    std::cout << "dVy" << std::setw(8) << std::fixed << std::setprecision(1) << Y0(5, 1) - Y_true(5, 1) << " [m/s]\n";
    std::cout << "dVz" << std::setw(8) << std::fixed << std::setprecision(1) << Y0(6, 1) - Y_true(6, 1) << " [m/s]\n";

    return 0;
}