#ifndef _CHEB3D_
#define _CHEB3D_

#include "..\include\matrix.hpp"
#include <cmath>
#include <iostream>

using namespace std;

/**
 * @details Evalúa una aproximación de Chebyshev 3D en un tiempo dado.
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
Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);

#endif