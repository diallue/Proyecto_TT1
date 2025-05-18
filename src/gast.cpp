#include "..\include\gast.hpp"

double gast(double Mjd_UT1) {
    return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), PI2);
}