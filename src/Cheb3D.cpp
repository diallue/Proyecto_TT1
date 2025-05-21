#include "..\include\Cheb3D.hpp"

/**
 * Evalúa una aproximación de Chebyshev 3D en un tiempo dado.
 * @param t Tiempo de evaluación (debe estar en [Ta, Tb]).
 * @param N Grado del polinomio de Chebyshev.
 * @param Ta Límite inferior del intervalo de tiempo.
 * @param Tb Límite superior del intervalo de tiempo.
 * @param Cx Coeficientes de Chebyshev para la componente x (Nx1).
 * @param Cy Coeficientes de Chebyshev para la componente y (Nx1).
 * @param Cz Coeficientes de Chebyshev para la componente z (Nx1).
 * @return Vector de posición aproximado (3x1).
 * @throw Error si t está fuera de [Ta, Tb] o dimensiones inválidas.
 */
Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz) {
    if (t < Ta || t > Tb) {
        cout << "ERROR: Time out of range in Cheb3D::Value" << endl;
        exit(EXIT_FAILURE);
    }

    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1 = zeros(3, 1);
    Matrix f2 = zeros(3, 1);

    for (int i = N; i >= 2; i--) {
        Matrix old_f1 = f1;
        Matrix coef(3, 1);
        coef(1,1) = Cx(i,1);
        coef(2,1) = Cy(i,1);
        coef(3,1) = Cz(i,1);
        f1 = f1 * (2.0 * tau) - f2 + coef;
        f2 = old_f1;
    }

    Matrix coef0(3, 1);
    coef0(1,1) = Cx(1,1);
    coef0(2,1) = Cy(1,1);
    coef0(3,1) = Cz(1,1);
    Matrix ChebApp = f1 * tau - f2 + coef0;

    return ChebApp;
}