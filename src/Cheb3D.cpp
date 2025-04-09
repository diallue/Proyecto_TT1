#include "..\include\Cheb3D.hpp"

double Cheb3D(double t, double N, double Ta, double Tb, double Cx, double Cy, double Cz) {
	// Check validity
	if ( (t<Ta) || (Tb<t) ) {
		error('ERROR: Time out of range in Cheb3D::Value\n');
	}

	// Clenshaw algorithm
	double tau = (2*t-Ta-Tb)/(Tb-Ta);

	double f1 = zeros(1,3);
	double f2 = zeros(1,3);

	for (int i = N - 1; i >= -1; i--) {
		double old_f1 = f1;
		f1 = 2*tau*f1-f2+[Cx(i),Cy(i),Cz(i)];
		f2 = old_f1;
	}

	double ChebApp = tau*f1-f2+[Cx(1),Cy(1),Cz(1)];
	
	return ChebApp;
}