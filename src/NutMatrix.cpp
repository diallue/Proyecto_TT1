#include "..\include\NutMatrix.hpp"

Matrix NutMatrix(double Mjd_TT) {
	double eps = MeanObliquity(Mjd_TT);

    auto [dpsi, deps] = NutAngles(Mjd_TT);

    Matrix NutMat = R_x(-eps - deps) * R_z(-dpsi) * R_x(eps);

    return NutMat;
}