#include "..\include\R_y.hpp"

/**
 * Genera una matriz de rotación alrededor del eje Y.
 * @param angle Ángulo de rotación (en radianes).
 * @return Matriz de rotación 3x3.
 */
Matrix R_y(double angle) {
	double C = cos(angle);
	double S = sin(angle);
	Matrix rotmat = zeros(3,3);

	rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
	rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
	rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;
	
	return rotmat;
}