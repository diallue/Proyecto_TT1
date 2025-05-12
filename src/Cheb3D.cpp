#include "..\include\Cheb3D.hpp"

Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz) {
    if (t < Ta || t > Tb) {
        cout << "ERROR: Time out of range in Cheb3D::Value";
    }

    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1 = Matrix::zeros(3, 1);
    Matrix f2 = Matrix::zeros(3, 1);

    for (int i = N; i >= 2; --i) {
        Matrix old_f1 = f1;
        Matrix coef(3, 1);
        coef(1,1) = Cx(i,1);
        coef(2,1) = Cy(i,1);
        coef(3,1) = Cz(i,1);
        f1 = 2.0 * tau * f1 - f2 + coef;
        f2 = old_f1;
    }

    Matrix coef0(3, 1);
    coef0(1,1) = Cx(1,1);
    coef0(2,1) = Cy(1,1);
    coef0(3,1) = Cz(1,1);
    Matrix ChebApp = tau * f1 - f2 + coef0;

    return ChebApp;
}