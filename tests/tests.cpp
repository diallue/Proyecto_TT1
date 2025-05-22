#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"
#include <cstdio>
#include <cmath>
#include <iostream>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { std::cout << "Running " #test << std::endl; int r=test(); tests_run++; if(r) return r; } while(0)

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
    std::cout << "Starting m_constructor_01\n";
    Matrix A(3, 3);
    _assert(A.n_row == 3 && A.n_column == 3);
    std::cout << "Finished m_constructor_01\n";
    return 0;
}

int m_constructor_02() {
    std::cout << "Starting m_constructor_02\n";
    Matrix A;
    _assert(A.n_row == 0 && A.n_column == 0);
    std::cout << "Finished m_constructor_02\n";
    return 0;
}

int m_zeros_01() {
    std::cout << "Starting m_zeros_01\n";
    int f = 3;
    int c = 3; // Corregido de c=4 a c=3
    Matrix A(f, c);
    for (int i = 1; i <= f; i++) {
        for (int j = 1; j <= c; j++) {
            A(i,j) = 0;
        }
    }
    Matrix B = zeros(f, c);
    _assert(m_equals(A, B, 1e-10));
    std::cout << "Finished m_zeros_01\n";
    return 0;
}

int m_zeros_02() {
    std::cout << "Starting m_zeros_02\n";
    int f = 3;
    int c = 3; // Corregido de c=4 a c=3
    Matrix A(f, c);
    for (int i = 1; i <= f; i++) {
        for (int j = 1; j <= c; j++) {
            A(i,j) = 0;
        }
    }
    Matrix B = zeros(f, c);
    _assert(m_equals(A, B, 1e-10));
    std::cout << "Finished m_zeros_02\n";
    return 0;
}

int m_sum_01() {
    std::cout << "Starting m_sum_01\n";
    int f = 3;
    int c = 3; // Corregido de c=4 a c=3
    Matrix A(f, c);
    A(1,1) = 0; A(1,2) =  2; A(1,3) = 8;
    A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
    A(3,1) = 0; A(3,2) =  1; A(3,3) = 5;
    Matrix B(f, c);
    B(1,1) = 2; B(1,2) =  0; B(1,3) = 0;
    B(2,1) = 7; B(2,2) = -2; B(2,3) = 1;
    B(3,1) = 0; B(3,2) = -3; B(3,3) = 2;
    Matrix C(f, c);
    C(1,1) = 2; C(1,2) =  2; C(1,3) = 8;
    C(2,1) = 8; C(2,2) = -3; C(2,3) = 1;
    C(3,1) = 0; C(3,2) = -2; C(3,3) = 7;
    Matrix R = A + B;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_sum_01\n";
    return 0;
}

int m_scalar_sum_01() {
    std::cout << "Starting m_scalar_sum_01\n";
    Matrix A(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    C(1,1)=6; C(1,2)=7; C(2,1)=8; C(2,2)=9;
    Matrix R = A + 5;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_scalar_sum_01\n";
    return 0;
}

int m_sub_01() {
    std::cout << "Starting m_sub_01\n";
    int f = 3;
    int c = 3; // Corregido de c=4 a c=3
    Matrix A(f, c);
    A(1,1) = 0; A(1,2) = 2; A(1,3) = 8;
    A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
    A(3,1) = 0; A(3,2) = 1; A(3,3) = 5;
    Matrix B(f, c);
    B(1,1) = 2; B(1,2) = 0; B(1,3) = 0;
    B(2,1) = 7; B(2,2) = -2; B(2,3) = 1;
    B(3,1) = 0; B(3,2) = -3; B(3,3) = 2;
    Matrix C(f, c);
    C(1,1) = -2; C(1,2) = 2; C(1,3) = 8;
    C(2,1) = -6; C(2,2) = 1; C(2,3) = -1;
    C(3,1) = 0; C(3,2) = 4; C(3,3) = 3;
    Matrix R = A - B;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_sub_01\n";
    return 0;
}

int m_scalar_sub_01() {
    std::cout << "Starting m_scalar_sub_01\n";
    Matrix A(2, 2), C(2, 2);
    A(1,1)=5; A(1,2)=7; A(2,1)=9; A(2,2)=11;
    C(1,1)=3; C(1,2)=5; C(2,1)=7; C(2,2)=9;
    Matrix R = A - 2;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_scalar_sub_01\n";
    return 0;
}

int m_mult_01() {
    std::cout << "Starting m_mult_01\n";
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
    std::cout << "Finished m_mult_01\n";
    return 0;
}

int m_scalar_mult_01() {
    std::cout << "Starting m_scalar_mult_01\n";
    Matrix A(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    C(1,1)=3; C(1,2)=6; C(2,1)=9; C(2,2)=12;
    Matrix R = A * 3;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_scalar_mult_01\n";
    return 0;
}

int m_div_01() {
    std::cout << "Starting m_div_01\n";
    Matrix A(2, 2), B(2, 2), C(2, 2);
    A(1,1)=1; A(1,2)=2; A(2,1)=3; A(2,2)=4;
    B(1,1)=4; B(1,2)=3; B(2,1)=2; B(2,2)=1;
    Matrix R = A / B;
    Matrix expected = A * inv(B);
    _assert(m_equals(R, expected, 1e-10));
    std::cout << "Finished m_div_01\n";
    return 0;
}

int m_scalar_div_01() {
    std::cout << "Starting m_scalar_div_01\n";
    Matrix A(2, 2), C(2, 2);
    A(1,1)=2; A(1,2)=4; A(2,1)=6; A(2,2)=8;
    C(1,1)=1; C(1,2)=2; C(2,1)=3; C(2,2)=4;
    Matrix R = A / 2;
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_scalar_div_01\n";
    return 0;
}

int m_transpose_01() {
    std::cout << "Starting m_transpose_01\n";
    Matrix A(2, 3), C(3, 2);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    C(1,1)=1; C(1,2)=4;
    C(2,1)=2; C(2,2)=5;
    C(3,1)=3; C(3,2)=6;
    Matrix R = transpose(A);
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_transpose_01\n";
    return 0;
}

int m_inv_01() {
    std::cout << "Starting m_inv_01\n";
    Matrix A(2, 2), I(2, 2);
    A(1,1)=4; A(1,2)=7;
    A(2,1)=2; A(2,2)=6;
    Matrix A_inv = inv(A);
    Matrix product = A * A_inv;
    I(1,1)=1; I(1,2)=0;
    I(2,1)=0; I(2,2)=1;
    _assert(m_equals(product, I, 1e-10));
    std::cout << "Finished m_inv_01\n";
    return 0;
}

int m_eye_01() {
    std::cout << "Starting m_eye_01\n";
    Matrix I = Matrix::eye(3);
    _assert(I(1,1)==1 && I(1,2)==0 && I(1,3)==0 &&
            I(2,1)==0 && I(2,2)==1 && I(2,3)==0 &&
            I(3,1)==0 && I(3,2)==0 && I(3,3)==1);
    std::cout << "Finished m_eye_01\n";
    return 0;
}

int m_norm_01() {
    std::cout << "Starting m_norm_01\n";
    Matrix A(2, 2);
    A(1,1)=3; A(1,2)=4;
    A(2,1)=0; A(2,2)=0;
    _assert(fabs(A.norm() - 5) < 1e-10);
    std::cout << "Finished m_norm_01\n";
    return 0;
}

int m_dot_01() {
    std::cout << "Starting m_dot_01\n";
    Matrix A(1, 3), B(1, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    B(1,1)=4; B(1,2)=5; B(1,3)=6;
    _assert(fabs(A.dot(B) - 32) < 1e-10);
    std::cout << "Finished m_dot_01\n";
    return 0;
}

int m_cross_01() {
    std::cout << "Starting m_cross_01\n";
    Matrix A(3, 1), B(3, 1), C(3, 1);
    A(1,1)=1; A(2,1)=0; A(3,1)=0;
    B(1,1)=0; B(2,1)=1; B(3,1)=0;
    C(1,1)=0; C(2,1)=0; C(3,1)=1;
    Matrix R = A.cross(B);
    _assert(m_equals(C, R, 1e-10));
    std::cout << "Finished m_cross_01\n";
    return 0;
}

int m_extract_row_01() {
    std::cout << "Starting m_extract_row_01\n";
    Matrix A(2, 3), row(1, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    row(1,1)=4; row(1,2)=5; row(1,3)=6;
    Matrix R = A.extract_row(2);
    _assert(m_equals(row, R, 1e-10));
    std::cout << "Finished m_extract_row_01\n";
    return 0;
}

int m_extract_column_01() {
    std::cout << "Starting m_extract_column_01\n";
    Matrix A(2, 3), col(2, 1);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    col(1,1)=2; col(2,1)=5;
    Matrix R = A.extract_column(2);
    _assert(m_equals(col, R, 1e-10));
    std::cout << "Finished m_extract_column_01\n";
    return 0;
}

int m_assign_row_01() {
    std::cout << "Starting m_assign_row_01\n";
    Matrix A(2, 3), row(1, 3), expected(2, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    row(1,1)=7; row(1,2)=8; row(1,3)=9;
    expected(1,1)=1; expected(1,2)=2; expected(1,3)=3;
    expected(2,1)=7; expected(2,2)=8; expected(2,3)=9;
    A.assign_row(2, row);
    _assert(m_equals(expected, A, 1e-10));
    std::cout << "Finished m_assign_row_01\n";
    return 0;
}

int m_assign_column_01() {
    std::cout << "Starting m_assign_column_01\n";
    Matrix A(2, 3), col(2, 1), expected(2, 3);
    A(1,1)=1; A(1,2)=2; A(1,3)=3;
    A(2,1)=4; A(2,2)=5; A(2,3)=6;
    col(1,1)=7; col(2,1)=8;
    expected(1,1)=1; expected(1,2)=7; expected(1,3)=3;
    expected(2,1)=4; expected(2,2)=8; expected(2,3)=6;
    A.assign_column(2, col);
    _assert(m_equals(expected, A, 1e-10));
    std::cout << "Finished m_assign_column_01\n";
    return 0;
}

int m_union_vector_01() {
    std::cout << "Starting m_union_vector_01\n";
    Matrix A(2, 2), B(2, 2), expected(2, 4);
    A(1,1)=1; A(1,2)=2;
    A(2,1)=3; A(2,2)=4;
    B(1,1)=5; B(1,2)=6;
    B(2,1)=7; B(2,2)=8;
    expected(1,1)=1; expected(1,2)=2; expected(1,3)=5; expected(1,4)=6;
    expected(2,1)=3; expected(2,2)=4; expected(2,3)=7; expected(2,4)=8;
    Matrix R = A.union_vector(B, true);
    _assert(m_equals(expected, R, 1e-10));
    std::cout << "Finished m_union_vector_01\n";
    return 0;
}

int accel_point_mass_01() {
    std::cout << "Starting accel_point_mass_01\n";
    Matrix r(3, 1), s(3, 1), expected(3, 1);
    std::cout << "Initializing r\n";
    r(1,1) = 7078e3; r(2,1) = 0; r(3,1) = 0;
    std::cout << "Initializing s\n";
    s(1,1) = 384400e3; s(2,1) = 0; s(3,1) = 0;
    double GM = 4.9048695e12;
    std::cout << "Initializing expected\n";
    expected(1,1) = 1.25702380070747e-06; expected(2,1) = 0; expected(3,1) = 0;
    std::cout << "Calling AccelPointMass\n";
    Matrix result = AccelPointMass(r, s, GM);
    std::cout << "Comparing results\n";
    _assert(m_equals(result, expected, 1e-10));
    std::cout << "Finished accel_point_mass_01\n";
    return 0;
}

int cheb3d_01() {
    std::cout << "Starting cheb3d_01\n";
    double t = 0;
    int N = 3;
    double Ta = -1;
    double Tb = 1;
    Matrix Cx(3, 1);
    Cx(1,1) = 1; Cx(2,1) = 0.5; Cx(3,1) = 0.2;
    Matrix Cy(3, 1);
    Cy(1,1) = 0; Cy(2,1) = 0.3; Cy(3,1) = 0.1;
    Matrix Cz(3, 1);
    Cz(1,1) = 0; Cz(2,1) = 0.1; Cz(3,1) = 0.05;
    Matrix expected(3, 1);
    expected(1,1) = 0.8; expected(2,1) = -0.1; expected(3,1) = -0.05;
    Matrix result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    _assert(m_equals(result, expected, 1e-10));
    std::cout << "Finished cheb3d_01\n";
    return 0;
}

int ecc_anom_01() {
    std::cout << "Starting ecc_anom_01\n";
    double M = 1.0;
    double e = 0.1;
    double expected = 1.08859775239789;
    std::cout << "Calling EccAnom\n";
    double result = EccAnom(M, e);
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-10);
    std::cout << "Finished ecc_anom_01\n";
    return 0;
}

int frac_01() {
    std::cout << "Starting frac_01\n";
    double x = 3.7;
    double expected = 0.7;
    std::cout << "Calling Frac\n";
    double result = Frac(x);
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-10);
    std::cout << "Finished frac_01\n";
    return 0;
}

int mean_obliquity_01() {
    std::cout << "Starting mean_obliquity_01\n";
    double Mjd_TT = 47;
    double expected = 0.409412778202273;
    std::cout << "Calling MeanObliquity\n";
    double result = MeanObliquity(Mjd_TT);
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-12);
    std::cout << "Finished mean_obliquity_01\n";
    return 0;
}

int mjday_test_01() {
    std::cout << "Starting mjday_test_01\n";
    std::cout << "Calling Mjday\n";
    double result = Mjday(2023, 1, 1, 12, 0, 0);
    double expected = 59945.5;
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-9);
    std::cout << "Finished mjday_test_01\n";
    return 0;
}

int mjday_tdb_test_01() {
    std::cout << "Starting mjday_tdb_test_01\n";
    std::cout << "Calling Mjday_TDB\n";
    double result = Mjday_TDB(58000.0);
    double expected = 57999.9999999833;
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-10);
    std::cout << "Finished mjday_tdb_test_01\n";
    return 0;
}

int position_test_01() {
    std::cout << "Starting position_test_01\n";
    std::cout << "Calling Position\n";
    Matrix result = Position(10.0, 45.0, 1000.0);
    std::cout << "Initializing expected\n";
    Matrix expected(3, 1);
    expected(1,1) = -2818651.27539222;
    expected(2,1) = -1827503.07323191;
    expected(3,1) =  5404810.30792904;
    double tolerance = 1e-8;
    std::cout << "Comparing results\n";
    _assert(fabs(result(1,1) - expected(1,1)) < tolerance);
    _assert(fabs(result(2,1) - expected(2,1)) < tolerance);
    _assert(fabs(result(3,1) - expected(3,1)) < tolerance);
    std::cout << "Finished position_test_01\n";
    return 0;
}

int rx_test_01() {
    std::cout << "Starting rx_test_01\n";
    double angle = 4.0;
    std::cout << "Calling R_x\n";
    Matrix result = R_x(angle);
    std::cout << "Initializing expected\n";
    Matrix expected(3, 3);
    expected(1,1) = 1.0; expected(1,2) = 0.0; expected(1,3) = 0.0;
    expected(2,1) = 0.0; expected(2,2) = -0.653643620863612; expected(2,3) = -0.756802495307928;
    expected(3,1) = 0.0; expected(3,2) = 0.756802495307928; expected(3,3) = -0.653643620863612;
    std::cout << "Comparing results\n";
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            _assert(fabs(result(i,j) - expected(i,j)) < 1e-15);
        }
    }
    std::cout << "Finished rx_test_01\n";
    return 0;
}

int ry_test_01() {
    std::cout << "Starting ry_test_01\n";
    double angle = 4.0;
    std::cout << "Calling R_y\n";
    Matrix result = R_y(angle);
    std::cout << "Initializing expected\n";
    Matrix expected(3,3);
    expected(1,1) = -0.653643620863612; expected(1,2) = 0.0; expected(1,3) = 0.756802495307928;
    expected(2,1) = 0.0; expected(2,2) = 1.0; expected(2,3) = 0.0;
    expected(3,1) = -0.756802495307928; expected(3,2) = 0.0; expected(3,3) = -0.653643620863612;
    std::cout << "Comparing results\n";
    _assert(m_equals(result, expected, 1e-15));
    std::cout << "Finished ry_test_01\n";
    return 0;
}

int rz_test_01() {
    std::cout << "Starting rz_test_01\n";
    double angle = 4.0;
    std::cout << "Calling R_z\n";
    Matrix result = R_z(angle);
    std::cout << "Initializing expected\n";
    Matrix expected(3,3);
    expected(1,1) = -0.653643620863612; expected(1,2) = -0.756802495307928; expected(1,3) = 0.0;
    expected(2,1) = 0.756802495307928; expected(2,2) = -0.653643620863612; expected(2,3) = 0.0;
    expected(3,1) = 0.0; expected(3,2) = 0.0; expected(3,3) = 1.0;
    std::cout << "Comparing results\n";
    _assert(m_equals(result, expected, 1e-15));
    std::cout << "Finished rz_test_01\n";
    return 0;
}

int sign_test_01() {
    std::cout << "Starting sign_test_01\n";
    std::cout << "Calling sign_\n";
    double result = sign_(5.0, -3.0);
    double expected = -5.0;
    std::cout << "Comparing results\n";
    _assert(fabs(result - expected) < 1e-10);
    std::cout << "Finished sign_test_01\n";
    return 0;
}

int timediff_test_01() {
    std::cout << "Starting timediff_test_01\n";
    std::cout << "Calling timediff\n";
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(0.2, 32.0);
    std::cout << "Comparing results\n";
    _assert(fabs(UT1_TAI - (-31.8)) < 1e-10);  
    _assert(fabs(UTC_GPS - (-13.0)) < 1e-10);  
    _assert(fabs(UT1_GPS - (-12.8)) < 1e-10);  
    _assert(fabs(TT_UTC - 64.184) < 1e-10);    
    _assert(fabs(GPS_UTC - 13.0) < 1e-10);     
    std::cout << "Finished timediff_test_01\n";
    return 0;
}

int azelpa_test_01() {
    std::cout << "Starting azelpa_test_01\n";
	
    Matrix s(3, 1);
	
    s(1,1) = 1000.0; s(2,1) = 2000.0; s(3,1) = 3000.0;
	
    auto [Az, El, dAds, dEds] = AzElPa(s);
    double expected_Az = 0.463647609000806;
    double expected_El = 0.930274014115472;
    double rho = sqrt(1000.0*1000.0 + 2000.0*2000.0);
    double rho_sq = rho*rho; 
	
    Matrix expected_dAds(1, 3);
    expected_dAds(1,1) = 2000.0/rho_sq;   
    expected_dAds(1,2) = -1000.0/rho_sq;  
    expected_dAds(1,3) = 0.0;
	
    Matrix expected_dEds(1, 3);
    expected_dEds(1,1) = -9.5831484749991e-05;
    expected_dEds(1,2) = -0.000191662969499982;
    expected_dEds(1,3) = 0.000159719141249985;
	
    double tolerance = 1e-12;
	
    _assert(fabs(Az - expected_Az) < tolerance);
    _assert(fabs(El - expected_El) < tolerance);
	
    for (int j = 1; j <= 3; j++) {
        _assert(fabs(dAds(1,j) - expected_dAds(1,j)) < tolerance);
        _assert(fabs(dEds(1,j) - expected_dEds(1,j)) < 1e-8);
    }
	
    std::cout << "Finished azelpa_test_01\n";
    return 0;
}

int iers_test_01() {
    std::cout << "Starting iers_test_01\n";

    Matrix eop(13, 2);
    eop(4, 1) = 49746.0; 
    eop(5, 1) = -5.59518621231704e-07 * ARCS;
    eop(6, 1) = 2.33458634442529e-06 * ARCS;
    eop(7, 1) = 0.3260677;                   
    eop(8, 1) = 0.0027213;                   
    eop(9, 1) = -1.16864337831454e-07 * ARCS; 
    eop(10, 1) = -2.48709418409192e-08 * ARCS; 
    eop(11, 1) = -8.19335121075116e-10 * ARCS; 
    eop(12, 1) = -1.53201123230613e-09 * ARCS; 
    eop(13, 1) = 29.0;                      
    eop(4, 2) = 49747.0;
    eop(5, 2) = eop(5, 1);
    eop(6, 2) = eop(6, 1);
    eop(7, 2) = eop(7, 1);
    eop(8, 2) = eop(8, 1);
    eop(9, 2) = eop(9, 1);
    eop(10, 2) = eop(10, 1);
    eop(11, 2) = eop(11, 1);
    eop(12, 2) = eop(12, 1);
    eop(13, 2) = eop(13, 1);

    double Mjd_UTC = 49746.0;
    auto result = IERS(eop, Mjd_UTC, 'l');

    double expected_x_pole = -5.59518621231704e-07;
    double expected_y_pole = 2.33458634442529e-06;
    double expected_UT1_UTC = 0.3260677;
    double expected_LOD = 0.0027213;
    double expected_dpsi = -1.16864337831454e-07;
    double expected_deps = -2.48709418409192e-08;
    double expected_dx_pole = -8.19335121075116e-10;
    double expected_dy_pole = -1.53201123230613e-09;
    double expected_TAI_UTC = 29.0;
	
    _assert(std::abs(std::get<0>(result) - expected_x_pole) < 1e-12);
    _assert(std::abs(std::get<1>(result) - expected_y_pole) < 1e-12);
    _assert(std::abs(std::get<2>(result) - expected_UT1_UTC) < 1e-7);
    _assert(std::abs(std::get<3>(result) - expected_LOD) < 1e-7);
    _assert(std::abs(std::get<4>(result) - expected_dpsi) < 1e-12);
    _assert(std::abs(std::get<5>(result) - expected_deps) < 1e-12);
    _assert(std::abs(std::get<6>(result) - expected_dx_pole) < 1e-12);
    _assert(std::abs(std::get<7>(result) - expected_dy_pole) < 1e-12);
    _assert(std::abs(std::get<8>(result) - expected_TAI_UTC) < 1e-7);

    std::cout << "Finished iers_test_01\n";
    return 0;
}

int legendre_test_01() {
    std::cout << "Starting legendre_test_01\n";
    int n = 3;
    int m = 3;
    double fi = 3.141592653589793 / 4.0;

    auto [pnm, dpnm] = Legendre(n, m, fi);

    Matrix expected_pnm(4, 4);
    expected_pnm(1, 1) = 1.0;
    expected_pnm(1, 2) = 0.0;
    expected_pnm(1, 3) = 0.0;
    expected_pnm(1, 4) = 0.0;
    expected_pnm(2, 1) = 1.22474487139159;
    expected_pnm(2, 2) = 1.22474487139159;
    expected_pnm(2, 3) = 0.0;
    expected_pnm(2, 4) = 0.0;
    expected_pnm(3, 1) = 0.559016994374947;
    expected_pnm(3, 2) = 1.93649167310371;
    expected_pnm(3, 3) = 0.968245836551854;
    expected_pnm(3, 4) = 0.0;
    expected_pnm(4, 1) = -0.467707173346743;
    expected_pnm(4, 2) = 1.71846588560844;
    expected_pnm(4, 3) = 1.81142209327368;
    expected_pnm(4, 4) = 0.739509972887452;

    Matrix expected_dpnm(4, 4);
    expected_dpnm(1, 1) = 0.0;
    expected_dpnm(1, 2) = 0.0;
    expected_dpnm(1, 3) = 0.0;
    expected_dpnm(1, 4) = 0.0;
    expected_dpnm(2, 1) = 1.22474487139159;
    expected_dpnm(2, 2) = -1.22474487139159;
    expected_dpnm(2, 3) = 0.0;
    expected_dpnm(2, 4) = 0.0;
    expected_dpnm(3, 1) = 3.35410196624968;
    expected_dpnm(3, 2) = 7.44760245974182e-16;
    expected_dpnm(3, 3) = -1.93649167310371;
    expected_dpnm(3, 4) = 0.0;
    expected_dpnm(4, 1) = 4.20936456012068;
    expected_dpnm(4, 2) = 4.00975373308636;
    expected_dpnm(4, 3) = -1.81142209327368;
    expected_dpnm(4, 4) = -2.21852991866236;

    _assert(m_equals(pnm, expected_pnm, 1e-10));
    _assert(m_equals(dpnm, expected_dpnm, 1e-10));

    std::cout << "Finished legendre_test_01\n";
    return 0;
}

int nutangles_test_01() {
    std::cout << "Starting nutangles_test_01\n";
    
    double Mjd_TT = 60000.0;
    auto [dpsi, deps] = NutAngles(Mjd_TT);
    
    double expected_dpsi = -4.49660202989237e-05;
    double expected_deps = 3.750608080778e-05;
    
    double tolerance = 1e-15;
    
    _assert(fabs(dpsi - expected_dpsi) < tolerance);
    _assert(fabs(deps - expected_deps) < tolerance);
    
    std::cout << "Finished nutangles_test_01\n";
    return 0;
}

int timeupdate_test_01() {
    Matrix P(2, 2);
    P(1,1) = 1.0; P(1,2) = 0.0;
    P(2,1) = 0.0; P(2,2) = 1.0;

    Matrix Phi(2, 2);
    Phi(1,1) = 1.0; Phi(1,2) = 0.1;
    Phi(2,1) = 0.0; Phi(2,2) = 1.0;

    Matrix Qdt(2, 2);
    Qdt(1,1) = 0.01; Qdt(1,2) = 0.0;
    Qdt(2,1) = 0.0; Qdt(2,2) = 0.01;

    TimeUpdate(P, Phi, Qdt);

    Matrix expected(2, 2);
    expected(1,1) = 1.02; expected(1,2) = 0.1;
    expected(2,1) = 0.1; expected(2,2) = 1.01;

    for(int i=1; i<=2; i++) {
        for(int j=1; j<=2; j++) {
            if(fabs(P(i,j)-expected(i,j)) > 1e-10) {
                return 1;
            }
        }
    }
    return 0;
}

int accel_harmonic_test_01() {
    std::cout << "Starting accel_harmonic_test_01\n";
    
    Matrix r(3, 1);
    r(1,1) = 6378136.3; r(2,1) = 0.0; r(3,1) = 0.0;
    
    Matrix E = Matrix::eye(3);
    
    int n_max = 4;
    int m_max = 4;
    
    Matrix expected(3, 1);
    expected(1,1) = -9.81424910913625;
    expected(2,1) = 3.49605441401039e-05;
    expected(3,1) = 0.000106491609085503;
    
    Matrix result = AccelHarmonic(r, E, n_max, m_max);
    
    _assert(m_equals(result, expected, 1e-10));
    
    std::cout << "Finished accel_harmonic_test_01\n";
    return 0;
}

int eqn_equinox_test_01() {
    std::cout << "Starting eqn_equinox_test_01\n";
    double Mjd_TT = 60000.0;
    double expected = -4.12564567660139e-05;
    double tolerance = 1e-12;

    double EqE = EqnEquinox(Mjd_TT);

    _assert(fabs(EqE - expected) < tolerance);

    std::cout << "Finished eqn_equinox_test_01\n";
    return 0;
}

int jpl_eph_de430_test_01() {
    std::cout << "Starting jpl_eph_de430_test_01\n";
	
    Matrix r_Mercury_exp(3, 1);
    r_Mercury_exp(1, 1) = 18364570296.7767;
    r_Mercury_exp(2, 1) = -194240890308.454;
    r_Mercury_exp(3, 1) = -89580185100.4172;
	
    Matrix r_Venus_exp(3, 1);
    r_Venus_exp(1, 1) = 134236480603.468;
    r_Venus_exp(2, 1) = -121788430081.357;
    r_Venus_exp(3, 1) = -59451902859.1386;
	
    Matrix r_Earth_exp(3, 1);
    r_Earth_exp(1, 1) = -26742322957.984;
    r_Earth_exp(2, 1) = 133827327107.372;
    r_Earth_exp(3, 1) = 58018326578.9664;
	
    Matrix r_Mars_exp(3, 1);
    r_Mars_exp(1, 1) = -170687646029.87;
    r_Mars_exp(2, 1) = -255905377504.909;
    r_Mars_exp(3, 1) = -108721468283.956;
	
    Matrix r_Jupiter_exp(3, 1);
    r_Jupiter_exp(1, 1) = 105439249049.906;
    r_Jupiter_exp(2, 1) = -847168695382.084;
    r_Jupiter_exp(3, 1) = -365696885445.354;
	
    Matrix r_Saturn_exp(3, 1);
    r_Saturn_exp(1, 1) = 594596636534.412;
    r_Saturn_exp(2, 1) = -1408094254368.45;
    r_Saturn_exp(3, 1) = -608810768966.399;
	
    Matrix r_Uranus_exp(3, 1);
    r_Uranus_exp(1, 1) = 2453302678658.18;
    r_Uranus_exp(2, 1) = 1439200900737.1;
    r_Uranus_exp(3, 1) = 596605466768.131;
	
    Matrix r_Neptune_exp(3, 1);
    r_Neptune_exp(1, 1) = 4400884479114.87;
    r_Neptune_exp(2, 1) = -974205297187;
    r_Neptune_exp(3, 1) = -510890407765.413;
	
    Matrix r_Pluto_exp(3, 1);
    r_Pluto_exp(1, 1) = 1967686040254.54;
    r_Pluto_exp(2, 1) = -4414806220850.92;
    r_Pluto_exp(3, 1) = -1978780908608.13;
	
    Matrix r_Moon_exp(3, 1);
    r_Moon_exp(1, 1) = 398673022.292818;
    r_Moon_exp(2, 1) = -38480993.8201334;
    r_Moon_exp(3, 1) = -55662807.01252;
	
    Matrix r_Sun_exp(3, 1);
    r_Sun_exp(1, 1) = 26173432883.2097;
    r_Sun_exp(2, 1) = -132807686289.682;
    r_Sun_exp(3, 1) = -57572484281.2736;
	
    auto [r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun] = JPL_Eph_DE430(58849.5);
    
    _assert(m_equals(r_Mercury, r_Mercury_exp, abs(r_Mercury_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Venus, r_Venus_exp, abs(r_Venus_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Earth, r_Earth_exp, abs(r_Earth_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Mars, r_Mars_exp, abs(r_Mars_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Jupiter, r_Jupiter_exp, abs(r_Jupiter_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Saturn, r_Saturn_exp, abs(r_Saturn_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Uranus, r_Uranus_exp, abs(r_Uranus_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Neptune, r_Neptune_exp, abs(r_Neptune_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Pluto, r_Pluto_exp, abs(r_Pluto_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Moon, r_Moon_exp, abs(r_Moon_exp(1, 1)*1e-11)));
    _assert(m_equals(r_Sun, r_Sun_exp, abs(r_Sun_exp(1, 1)*1e-11)));
    
    std::cout << "Finished jpl_eph_de430_test_01\n";
    return 0;
}

int ltc_test_01() {
    std::cout << "Starting ltc_test_01\n";
    double lon = 0.56;
    double lat = 0.76;
    
    Matrix M = LTC(lon, lat);
    
    Matrix expected(3, 3);
    expected(1,1) = -0.531186197920883; 
	expected(1,2) = 0.847255111013416;  
	expected(1,3) = 0.0;
    expected(2,1) = -0.583692215456663; 
	expected(2,2) = -0.365945563094434; 
	expected(2,3) = 0.724836010740905;
    expected(3,1) = 0.614121014746807;  
	expected(3,2) = 0.385022884661602;  
	expected(3,3) = 0.688921445110551;
    
    _assert(M.n_row == expected.n_row);
    _assert(M.n_column == expected.n_column);
    
    for(int i = 1; i <= 3; i++) {
        for(int j = 1; j <= 3; j++) {
            _assert(fabs(M(i,j) - expected(i,j)) < 1e-12);
        }
    }
    
    std::cout << "Finished ltc_test_01\n";
    return 0;
}

int nut_matrix_test_01() {
    std::cout << "Starting nut_matrix_test_01\n";
    double Mjd_TT = 60000.0;
    
    Matrix NutMat = NutMatrix(Mjd_TT);
    
    Matrix expected(3, 3);
    expected(1,1) = 0.999999998989028; 
	expected(1,2) = 4.12564567521108e-05; 
	expected(1,3) = 1.78842879763684e-05;
    expected(2,1) = -4.12557859535431e-05; 
	expected(2,2) = 0.999999998445613; 
	expected(2,3) = -3.75064497141753e-05;
    expected(3,1) = -1.78858353317898e-05; 
	expected(3,2) = 3.75057118458533e-05; 
	expected(3,3) = 0.999999999136709;
    
    _assert(NutMat.n_row == expected.n_row);
    _assert(NutMat.n_column == expected.n_column);
    
    for(int i = 1; i <= 3; i++) {
        for(int j = 1; j <= 3; j++) {
            _assert(fabs(NutMat(i,j) - expected(i,j)) < 1e-12);
        }
    }
    
    std::cout << "Finished nut_matrix_test_01\n";
    return 0;
}

int pole_matrix_test_01() {
    std::cout << "Starting pole_matrix_test_01\n";
    double xp = 0.0001;
    double yp = 0.0002;
    
    Matrix PoleMat = PoleMatrix(xp, yp);
    
    Matrix expected(3, 3);
    expected(1,1) = 0.999999995; 
	expected(1,2) = 1.99999998333333e-08; 
	expected(1,3) = 9.99999978333334e-05;
    expected(2,1) = 0.0; 
	expected(2,2) = 0.99999998; 
	expected(2,3) = -0.000199999998666667;
    expected(3,1) = -9.99999998333333e-05; 
	expected(3,2) = 0.000199999997666667; 
	expected(3,3) = 0.999999975;
    
    _assert(PoleMat.n_row == expected.n_row);
    _assert(PoleMat.n_column == expected.n_column);
    
    for(int i = 1; i <= 3; i++) {
        for(int j = 1; j <= 3; j++) {
            _assert(fabs(PoleMat(i,j) - expected(i,j)) < 1e-12);
        }
    }
    
    std::cout << "Finished pole_matrix_test_01\n";
    return 0;
}

int prec_matrix_test_01() {
    std::cout << "Starting prec_matrix_test_01\n";
    double Mjd_1 = 60000.0;
    double Mjd_2 = 60100.0;
    
    Matrix PrecMat = PrecMatrix(Mjd_1, Mjd_2);
    
    Matrix expected(3, 3);
    expected(1,1) = 0.999999997771519; 
	expected(1,2) = -6.12316906343852e-05; 
	expected(1,3) = -2.66015332954581e-05;
    expected(2,1) = 6.12316906343852e-05; 
	expected(2,2) = 0.99999999812534; 
	expected(2,3) = -8.14428812735772e-10;
    expected(3,1) = 2.66015332954582e-05; 
	expected(3,2) = -8.14428046226956e-10; 
	expected(3,3) = 0.999999999646179;
    
    _assert(PrecMat.n_row == expected.n_row);
    _assert(PrecMat.n_column == expected.n_column);
    
    for(int i = 1; i <= 3; i++) {
        for(int j = 1; j <= 3; j++) {
            _assert(fabs(PrecMat(i,j) - expected(i,j)) < 1e-12);
        }
    }
    
    std::cout << "Finished prec_matrix_test_01\n";
    return 0;
}

int gmst_test_01() {
    std::cout << "Starting gmst_test_01\n";
    double Mjd_UT1 = 60000.0;
    double expected = 2.69831296672573;

    double gmstime = gmst(Mjd_UT1);

    _assert(fabs(gmstime - expected) < 1e-12);

    std::cout << "Finished gmst_test_01\n";
    return 0;
}

int gast_test_01() {
    std::cout << "Starting gast_test_01\n";
    double Mjd_UT1 = 60000.0;
    double expected = 2.69827171026896;

    double result = gast(Mjd_UT1);

    _assert(fabs(result - expected) < 1e-12);

    std::cout << "Finished gast_test_01\n";
    return 0;
}

int meas_update_test_01() {
    std::cout << "Starting meas_update_test_01\n";
    
    Matrix x(3, 1);
    x(1,1) = 1.0; x(2,1) = 1.0; x(3,1) = 1.0;
    
    Matrix z(2, 1);
    z(1,1) = 1.5; z(2,1) = 2.3;
    
    Matrix g(2, 1);
    g(1,1) = 1.2; g(2,1) = 2.1;
    
    Matrix s(2, 1);
    s(1,1) = 0.1; s(2,1) = 0.2;
    
    Matrix G(2, 3);
    G(1,1) = 1.0; G(1,2) = 0.0; G(1,3) = 0.0;
    G(2,1) = 0.0; G(2,2) = 1.0; G(2,3) = 0.0;
    
    Matrix P = Matrix::eye(3);
    int n = 3;
	
    auto [K, x_new, P_new] = MeasUpdate(x, z, g, s, G, P, n);
    
    Matrix expected_K(3, 2);
    expected_K(1,1) = 0.99009900990099; expected_K(1,2) = 0.0;
    expected_K(2,1) = 0.0; expected_K(2,2) = 0.961538461538461;
    expected_K(3,1) = 0.0; expected_K(3,2) = 0.0;
    
    Matrix expected_x(3, 1);
    expected_x(1,1) = 1.2970297029703;
    expected_x(2,1) = 1.19230769230769;
    expected_x(3,1) = 1.0;
    
    Matrix expected_P(3, 3);
    expected_P(1,1) = 0.00990099009900991; expected_P(1,2) = 0.0; expected_P(1,3) = 0.0;
    expected_P(2,1) = 0.0; expected_P(2,2) = 0.0384615384615385; expected_P(2,3) = 0.0;
    expected_P(3,1) = 0.0; expected_P(3,2) = 0.0; expected_P(3,3) = 1.0;
    
    _assert(K.n_row == expected_K.n_row);
    _assert(K.n_column == expected_K.n_column);
    _assert(x_new.n_row == expected_x.n_row);
    _assert(x_new.n_column == expected_x.n_column);
    _assert(P_new.n_row == expected_P.n_row);
    _assert(P_new.n_column == expected_P.n_column);
    
    for(int i = 1; i <= K.n_row; i++) {
        for(int j = 1; j <= K.n_column; j++) {
            _assert(fabs(K(i,j) - expected_K(i,j)) < 1e-12);
        }
    }
    
    for(int i = 1; i <= x_new.n_row; i++) {
        _assert(fabs(x_new(i,1) - expected_x(i,1)) < 1e-12);
    }
    
    for(int i = 1; i <= P_new.n_row; i++) {
        for(int j = 1; j <= P_new.n_column; j++) {
            _assert(fabs(P_new(i,j) - expected_P(i,j)) < 1e-12);
        }
    }
    
    std::cout << "Finished meas_update_test_01\n";
    return 0;
}

int g_accelharmonic_test_01() {
    std::cout << "Starting g_accelharmonic_test_01\n";
    
    Matrix r(3, 1);
    r(1,1) = 7000e3; r(2,1) = 0; r(3,1) = 0;
    
    Matrix U = Matrix::eye(3);
    int n_max = 4;
    int m_max = 4;
    
    Matrix G_expected(3, 3);
    G_expected(1,1) = 2.33048175246608e-06; 
	G_expected(1,2) = -1.88080662155699e-11; 
	G_expected(1,3) = -5.18909359925601e-11;
    G_expected(2,1) = -1.88040301950557e-11; 
	G_expected(2,2) = -1.1636930401132e-06; 
	G_expected(2,3) = -5.90145164317381e-12;
    G_expected(3,1) = -5.18900588865551e-11; 
	G_expected(3,2) = -5.90145164995007e-12; 
	G_expected(3,3) = -1.16678871423606e-06;
    
    Matrix G = G_AccelHarmonic(r, U, n_max, m_max);
    
    _assert(m_equals(G, G_expected, 1e-15));
    
    std::cout << "Finished g_accelharmonic_test_01\n";
    return 0;
}

int gha_matrix_test_01() {
    std::cout << "Starting gha_matrix_test_01\n";
    
    Matrix GHAmat = GHAMatrix(49746.0);
    
    Matrix expected_GHAmat(3, 3);
    expected_GHAmat(1,1) = -0.612591243162461;
    expected_GHAmat(1,2) = 0.790399879048998;
    expected_GHAmat(1,3) = 0.0;
    expected_GHAmat(2,1) = -0.790399879048998;
    expected_GHAmat(2,2) = -0.612591243162461;
    expected_GHAmat(2,3) = 0.0;
    expected_GHAmat(3,1) = 0.0;
    expected_GHAmat(3,2) = 0.0;
    expected_GHAmat(3,3) = 1.0;
    
    _assert(m_equals(GHAmat, expected_GHAmat, 1e-6));
    
    std::cout << "Finished gha_matrix_test_01\n";
    return 0;
}

int accel_test_01() {
    std::cout << "Starting m_Accel_01\n";
    
    Matrix R(6, 1); 
    R(1,1) = 9.0;
    R(2,1) = 7.5;
    R(3,1) = 2.8;
    R(4,1) = 6.43475380190004e+126;
    R(5,1) = 1.39220669308953e+127;
    R(6,1) = 4.84358390766759e+127;
    
    Matrix A(6, 1); 
    A(1,1) = 4.6;
    A(2,1) = 3.8;
    A(3,1) = 0.0;
    A(4,1) = 9.0;
    A(5,1) = 7.5;
    A(6,1) = 2.8;
    A = transpose(A);  
    
    Matrix B = Accel(0.0, A);  
    
    _assert(m_equals(R, B, fabs(R(6,1) * 1e-11))); 
    
    std::cout << "Finished m_Accel_01\n";
    return 0;
}

int var_eqn_test_01() {
    std::cout << "Starting var_eqn_test_01\n";
    
    double x = 4.0;
    Matrix yPhi(42, 1);
    yPhi(1,1) = 1e6; 
	yPhi(2,1) = 2e6; 
	yPhi(3,1) = 3e6;
    yPhi(4,1) = 1e3; 
	yPhi(5,1) = 2e3; 
	yPhi(6,1) = 3e3;
    for (int i = 7; i <= 42; i++) {
        yPhi(i,1) = 0.0;
    }
    
    Matrix expected_yPhip(42, 1);
    expected_yPhip(1,1) = 1000.0;
    expected_yPhip(2,1) = 2000.0;
    expected_yPhip(3,1) = 3000.0;
    expected_yPhip(4,1) = -8.55483473704505;
    expected_yPhip(5,1) = -14.7587726037911;
    expected_yPhip(6,1) = -26.646858402194;
    for (int i = 7; i <= 42; i++) {
        expected_yPhip(i,1) = 0.0;
    }
    
    Matrix result = VarEqn(x, yPhi);
    
    _assert(m_equals(result, expected_yPhip, 1e-6));
    
    std::cout << "Finished var_eqn_test_01\n";
    return 0;
}

int deinteg_test_01() {
    std::cout << "Starting deinteg_test_01\n";

    Matrix expected(6, 1);
	expected(1, 1) = 5542555.89427452;
	expected(2, 1) = 3213514.83814162;
	expected(3, 1) = 3990892.92789074;
	expected(4, 1) = 5394.06894044389;
	expected(5, 1) = -2365.2129057402;
	expected(6, 1) = -7061.8448137347;
          
	Matrix A(6, 1);
	A(1, 1) = 6221397.62857869;
	A(2, 1) = 2867713.77965738;
	A(3, 1) = 3006155.98509949;
	A(4, 1) = 4645.04725161806;
	A(5, 1) = -2752.21591588204;
	A(6, 1) = -7507.99940987031;
	
	A = transpose(A);
	Matrix result = DEInteg(Accel, 0, -134.999991953373, 1e-13, 1e-6, 6, A);

    _assert(m_equals(expected, result, abs(expected(5)*1e-6)));

    std::cout << "Finished deinteg_test_01\n";
    return 0;
}

int all_tests() {
    _verify(m_constructor_01);
    _verify(m_constructor_02);
    _verify(m_zeros_01);
    _verify(m_zeros_02);
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
    _verify(accel_point_mass_01);
    _verify(cheb3d_01);
    _verify(ecc_anom_01);
    _verify(frac_01);
    _verify(mean_obliquity_01);
    _verify(mjday_test_01);
    _verify(mjday_tdb_test_01);
    _verify(position_test_01);
    _verify(rx_test_01);
    _verify(ry_test_01);
    _verify(rz_test_01);
    _verify(sign_test_01);
    _verify(timediff_test_01);
    _verify(azelpa_test_01);
	_verify(iers_test_01);
    _verify(legendre_test_01);
	_verify(nutangles_test_01);
	_verify(timeupdate_test_01);
	_verify(accel_harmonic_test_01);
	_verify(eqn_equinox_test_01);
	_verify(jpl_eph_de430_test_01);
	_verify(ltc_test_01);
	_verify(nut_matrix_test_01);
	_verify(pole_matrix_test_01);
	_verify(prec_matrix_test_01);
	_verify(gmst_test_01);
	_verify(gast_test_01);
	_verify(meas_update_test_01);
	_verify(g_accelharmonic_test_01);
	_verify(gha_matrix_test_01);
	_verify(accel_test_01);
	_verify(var_eqn_test_01);
	_verify(deinteg_test_01);
    return 0;
}

int main() {
    std::cout << "Starting tests\n";
	eop19620101(21413);
	GGM03S(16471);
	DE430Coeff(2285, 1020);
	AuxParamLoad();
	GEOS3(46);
    int result = all_tests();
    if (result == 0)
        printf("PASSED\n");
    else
        printf("FAILED\n");
    printf("Tests run: %d\n", tests_run);
    return result != 0;
}