#include "..\include\Frac.hpp"

Matrix& Frac(double x) {
	Matrix& res = x-floor(x);
	
	return res;
}