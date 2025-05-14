#include "..\include\PoleMatrix.hpp"

Matrix PoleMatrix(double xp, double yp) {
	Matrix PoleMat = R_y(-xp) * R_x(-yp);

    return PoleMat;
}