Matrix& sign_(double a, double b) {
	Matrix& result;
	
	if (b >= 0.0) {
		result = fabs(a);
	} else {
		result = -fabs(a);
	}
	
	return result;
}