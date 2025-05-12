#include "..\include\sign_.hpp"

double sign_(double a, double b) {
	double result;
	
	if (b >= 0.0) {
		result = fabs(a);
	} else {
		result = -fabs(a);
	}
	
	return result;
}