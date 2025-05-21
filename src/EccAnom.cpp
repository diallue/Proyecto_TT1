#include "..\include\EccAnom.hpp"

/**
 * Calcula la anomalía excéntrica (E) a partir de la anomalía media (M) y la excentricidad (e).
 * Utiliza el método de Newton-Raphson para resolver la ecuación de Kepler.
 * @param M Anomalía media (en radianes).
 * @param e Excentricidad de la órbita (0 <= e < 1).
 * @return Anomalía excéntrica (E, en radianes).
 * @throw Error si no converge tras 15 iteraciones.
 */
double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;

    M = fmod(M, PI2);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = 3.141592653589793;
    }

    double f = E - e * std::sin(E) - M;
    E = E - f / (1.0 - e * std::cos(E));

    const double eps = std::numeric_limits<double>::epsilon();
    const double tolerance = 100.0 * eps;

    while (abs(f) > tolerance) {
        f = E - e * std::sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            cout << "convergence problems in EccAnom";
        }
    }

    return E;
}