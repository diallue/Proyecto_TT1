#include "..\include\AccelPointMass.hpp"

double AccelPointMass(double r, double s, double GM) {
	// Relative position vector of satellite w.r.t. point mass 
	double d = r - s;

	// Acceleration 
	double a = -GM * ( d/(norm(d)^3) + s/(norm(s)^3) );
	
	return a;
}