#include "..\include\GHAMatrix.hpp"

Matrix GHAMatrix(double Mjd_UT1) {
    double gast_value = gast(Mjd_UT1);
    Matrix GHAmat = R_z(gast_value);
    return GHAmat;
}