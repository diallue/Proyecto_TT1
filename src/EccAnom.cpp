#include "..\include\R_x.hpp"

Matrix& EccAnom(double M, double e) {
	double maxit = 15;
	double i = 1;

	// Starting value
	Matrix M = mod(M, 2.0*pi);
	Matrix E = zeros(3,3);

	if (e<0.8) {
		E = M; 
	} else{
		E = pi;
	}

	double f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );

	// Iteration
	while (abs(f) > 1e2*eps) {  
		f = E - e*sin(E) - M;
		E = E - f / ( 1.0 - e*cos(E) );
		i = i+1;
		if (i==maxit) {
			error(' convergence problems in EccAnom');
		}  
	}
	
	return E;
}