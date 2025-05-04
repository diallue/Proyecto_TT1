#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_constructor_01() {
    Matrix A(3, 3);
    _assert(A.n_row == 3 && A.n_column == 3);
    return 0;
}

int m_constructor_02() {
    Matrix A;
    _assert(A.n_row == 0 && A.n_column == 0);
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_scalar_sum_01() {
    Matrix A(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    C(1,1)=6; C(1,2)=7; C(2,1)=8; C(2,2)=9;
    
    Matrix R = A + 5;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_scalar_sub_01() {
    Matrix A(2, 2), C(2, 2);
    A(1,1)=5; A(1,2)=7; A(2,1)=9; A(2,2)=11;
    C(1,1)=3; C(1,2)=5; C(2,1)=7; C(2,2)=9;
    
    Matrix R = A - 2;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_mult_01() {
    Matrix A(2, 3), B(3, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    B(1,1)=7; B(1,2)=8;
    B(2,1)=9; B(2,2)=10;
    B(3,1)=11; B(3,2)=12;
    
    C(1,1)=58; C(1,2)=64;
    C(2,1)=139; C(2,2)=154;
    
    Matrix R = A * B;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_scalar_mult_01() {
    Matrix A(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    C(1,1)=3; C(1,2)=6; C(2,1)=9; C(2,2)=12;
    
    Matrix R = A * 3;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_div_01() {
    Matrix A(2, 2), B(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    B(1,1)=4; B(1,2)=3; B(2,1)=2; B(2,2)=1;
    
    Matrix R = A / B;
    Matrix expected = A * inv(B);
    _assert(m_equals(R, expected, 1e-10));
    return 0;
}

int m_scalar_div_01() {
    Matrix A(2, 2), C(2, 2);
    A(1,1)=2; A(1,2)=4; A(2,1)=6; A(2,2)=8;
    C(1,1)=1; C(1,2)=2; C(2,1)=3; C(2,2)=4;
    
    Matrix R = A / 2;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_transpose_01() {
    Matrix A(2, 3), C(3, 2);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    C(1,1)=1; C(1,2)=4;
    C(2,1)=2; C(2,2)=5;
    C(3,1)=3; C(3,2)=6;
    
    Matrix R = transpose(A);
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_inv_01() {
    Matrix A(2, 2), I(2, 2);
    A(1,1)=4; A(1,2)=7;
    A(2,1)=2; A(2,2)=6;
    
    Matrix A_inv = inv(A);
    Matrix product = A * A_inv;
    I(1,1)=1; I(1,2)=0;
    I(2,1)=0; I(2,2)=1;
    
    _assert(m_equals(product, I, 1e-10));
    return 0;
}

int m_eye_01() {
    Matrix I = Matrix::eye(3);
    _assert(I(1,1)==1 && I(1,2)==0 && I(1,3)==0 &&
            I(2,1)==0 && I(2,2)==1 && I(2,3)==0 &&
            I(3,1)==0 && I(3,2)==0 && I(3,3)==1);
    return 0;
}

int m_norm_01() {
    Matrix A(2, 2);
    A(1,1)=3; A(1,2)=4;
    A(2,1)=0; A(2,2)=0;
    
    _assert(fabs(A.norm() - 5) < 1e-10);
    return 0;
}

int m_dot_01() {
    Matrix A(1, 3), B(1, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    B(1,1)=4; B(1,2)=5; B(1,3)=6;
    
    _assert(fabs(A.dot(B) - 32) < 1e-10);
    return 0;
}

int m_cross_01() {
    Matrix A(3, 1), B(3, 1), C(3, 1);
    A(1,1)=1; A(2,1)=0; A(3,1)=0;
    B(1,1)=0; B(2,1)=1; B(3,1)=0;
    C(1,1)=0; C(2,1)=0; C(3,1)=1;
    
    Matrix R = A.cross(B);
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_extract_row_01() {
    Matrix A(2, 3), row(1, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    row(1,1)=4; row(1,2)=5; row(1,3)=6;
    
    Matrix R = A.extract_row(2);
    _assert(m_equals(row, R, 1e-10));
    return 0;
}

int m_extract_column_01() {
    Matrix A(2, 3), col(2, 1);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    col(1,1)=2; col(2,1)=5;
    
    Matrix R = A.extract_column(2);
    _assert(m_equals(col, R, 1e-10));
    return 0;
}

int m_assign_row_01() {
    Matrix A(2, 3), row(1, 3), expected(2, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    row(1,1)=7; row(1,2)=8; row(1,3)=9;
    
    expected(1,1)=1; expected(1,2)=2; expected(1,3)=3;
    expected(2,1)=7; expected(2,2)=8; expected(2,3)=9;
    
    A.assign_row(2, row);
    _assert(m_equals(expected, A, 1e-10));
    return 0;
}

int m_assign_column_01() {
    Matrix A(2, 3), col(2, 1), expected(2, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    
    col(1,1)=7; col(2,1)=8;
    
    expected(1,1)=1; expected(1,2)=7; expected(1,3)=3;
    expected(2,1)=4; expected(2,2)=8; expected(2,3)=6;
    
    A.assign_column(2, col);
    _assert(m_equals(expected, A, 1e-10));
    return 0;
}

int m_union_vector_01() {
    Matrix A(2, 2), B(2, 2), expected(2, 4);
    A(1,1)=1; A(1,2)=2;
    A(2,1)=3; A(2,2)=4;
    
    B(1,1)=5; B(1,2)=6;
    B(2,1)=7; B(2,2)=8;
    
    expected(1,1)=1; expected(1,2)=2; expected(1,3)=5; expected(1,4)=6;
    expected(2,1)=3; expected(2,2)=4; expected(2,3)=7; expected(2,4)=8;
    
    Matrix R = A.union_vector(B, true);
    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

int all_tests() {
    _verify(m_constructor_01);
    _verify(m_constructor_02);
    _verify(m_sum_01);
    _verify(m_scalar_sum_01);
    _verify(m_sub_01);
    _verify(m_scalar_sub_01);
    _verify(m_mult_01);
    _verify(m_scalar_mult_01);
    _verify(m_div_01);
    _verify(m_scalar_div_01);
    _verify(m_transpose_01);
    _verify(m_inv_01);
    _verify(m_eye_01);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
    _verify(m_extract_row_01);
    _verify(m_extract_column_01);
    _verify(m_assign_row_01);
    _verify(m_assign_column_01);
    _verify(m_union_vector_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
