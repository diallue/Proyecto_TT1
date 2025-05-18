#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\gloabl.hpp"
#include <iostream>

using namespace std;

int main() {
	eop19620101(21413); // c = 21413
	cout << eop19620101;
	
	GGM03S(16471);
	
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