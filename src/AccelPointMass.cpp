#include "..\include\AccelPointMass.hpp"

Matrix AccelPointMass(Matrix& r, Matrix& s, double GM) {
    Matrix d = r - s;

    double norm_d_cubed = pow(d.norm(), 3);
    double norm_s_cubed = pow(s.norm(), 3);
	
    Matrix a = ((d / norm_d_cubed) + (s / norm_s_cubed)) * (-GM);

    return a;
}