#include "..\include\Frac.hpp"

/**
 * Calcula la parte fraccionaria de un número real.
 * @param x Número real de entrada.
 * @return Parte fraccionaria de x (x - floor(x)).
 */
double Frac(double x) {
    return x - floor(x);
}