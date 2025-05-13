#include "..\include\AzElPa.hpp"

tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s) {
	double rho = sqrt(s(1,1)*s(1,1)+s(2,1)*s(2,1));
	
	double Az = atan2(s(1,1), s(2,1));
	
	if (Az < 0.0) {
		Az = Az + PI2;
	}
	
	double El = atan(s(3,1)/rho);
	
	Matrix dAds(1, 3);
    dAds(1,1) = s(2,1) / (rho * rho);
    dAds(1,2) = -s(1,1) / (rho * rho);
    dAds(1,3) = 0.0;

    Matrix dEds(1, 3);
    double s_dot_s = s(1,1) * s(1,1) + s(2,1) * s(2,1) + s(3,1) * s(3,1);
    dEds(1,1) = -s(1,1) * s(3,1) / rho;
    dEds(1,2) = -s(2,1) * s(3,1) / rho;
    dEds(1,3) = rho;
    dEds = dEds / s_dot_s;
	
	return make_tuple(Az, El, dAds, dEds);
}