#include "..\include\AzElPa.hpp"

/**
 * Calcula azimut, elevación y sus matrices de derivadas parciales para un vector de posición dado.
 * @param s Vector de posición en coordenadas topocéntricas [x, y, z] (metros)
 * @return Tupla con:
 *         [0] Azimut (radianes) [0, 2pi)
 *         [1] Elevación (radianes)
 *         [2] Matriz de derivadas parciales del azimut respecto a s (1x3)
 *         [3] Matriz de derivadas parciales de la elevación respecto a s (1x3)
 */
tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s) {
    double rho = sqrt(s(1,1)*s(1,1) + s(2,1)*s(2,1));
    
    double Az = atan2(s(1,1), s(2,1));
    if (Az < 0.0) {
        Az += PI2;
    }
    
    double El = atan(s(3,1)/rho);
    
    Matrix dAds(1, 3);
    double rho_squared = rho * rho;
    dAds(1,1) = s(2,1)/rho_squared;
    dAds(1,2) = -s(1,1)/rho_squared;
    dAds(1,3) = 0.0;

    Matrix dEds(1, 3);
    double s_dot_s = s(1,1)*s(1,1) + s(2,1)*s(2,1) + s(3,1)*s(3,1);
    double denom = rho * s_dot_s;
    
    dEds(1,1) = (-s(1,1)*s(3,1)) / denom;
    dEds(1,2) = (-s(2,1)*s(3,1)) / denom;
    dEds(1,3) = rho / s_dot_s;
    
    return make_tuple(Az, El, dAds, dEds);
}