#ifndef _SIGN__
#define _SIGN__

#include "..\include\matrix.hpp"
#include <cmath>

/**
 * @details Devuelve el valor absoluto de 'a' con el signo de 'b'.
 * @param a Valor numérico cuyo signo será modificado.
 * @param b Valor numérico que determina el signo del resultado.
 * @return |a| si b >= 0, -|a| si b < 0.
 */
double sign_(double a, double b);

#endif