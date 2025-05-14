#include "..\include\EqnEquinox.hpp"

double EqnEquinox(double Mjd_TT) {
    auto [dpsi, deps] = NutAngles(Mjd_TT);

    double EqE = dpsi * cos(MeanObliquity(Mjd_TT));

    return EqE;
}