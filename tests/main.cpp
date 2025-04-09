#include "..\include\matrix.h"
#include "..\include\R_x.h"
#include "..\include\R_y.h"
#include "..\include\R_z.h"
#include <iostream>

using namespace std;

int main() {
	Matrix auxX = R_x(3);
	cout << "auxX\n" << auxX << "\n";
	
	Matrix auxY = R_y(3);
	cout << "auxY\n" << auxY << "\n";
	
	Matrix auxZ = R_z(3);
	cout << "auxZ\n" << auxZ << "\n";
	
    Matrix M1(3, 2);
	M1(1,1) = 5;
	
    Matrix M2(3, 2);
	M2(1,1) = -3;
	
    Matrix M3 = M1 - M2;

    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	
	cout << M1(1,1) << "\n";

    return 0;
}