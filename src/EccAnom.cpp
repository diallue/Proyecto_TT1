#include "..\include\R_x.hpp"

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;

    M = fmod(M, PI2);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = pi;
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