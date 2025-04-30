#include "..\include\AzElPa.hpp"

tuple<double, double, std::vector<double>, vector<double>> AzElPa(double s) {
	double rho = sqrt(s[0]*s[0]+s[1]*s[1]);
	
	double Az = atan2(s[0], s[1]);
	
	if (Az < 0.0) {
		Az = Az + PI2;
	}
	
	double El = atan(s[2]/rho);
	
	vector<double> dAds = {s[1]/(rho * rho), -s[0]/(rho * rho), 0.0 };
    vector<double> dEds = {-s[0]*s[2]/rho, -s[1]*s[2]/rho, rho } / dot(s, s);
	
	return {Az, El, dAds, dEds};
}