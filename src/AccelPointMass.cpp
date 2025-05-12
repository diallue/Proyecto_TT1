#include "..\include\AccelPointMass.hpp"

Matrix AccelPointMass(const Matrix& r, const Matrix& s, double GM) {
    Matrix d = r - s;

    double norm_d_cubed = pow(d.norm(), 3);
    double norm_s_cubed = pow(s.norm(), 3);
	
    Matrix a = -GM * ((d / norm_d_cubed) + (s / norm_s_cubed));

    return a;
}