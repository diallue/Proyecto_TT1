#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include <cmath>
#include <tuple>

/**
 * @details Calcula las posiciones de los planetas, la Luna y el Sol usando las efemérides JPL DE430
 * @param Mjd_TDB Fecha Juliana Modificada en Tiempo Dinámico Bariocéntrico (TDB)
 * @return Tupla con 11 vectores de posición (3x1) en metros:
 *         [0] Mercurio (geocéntrico)
 *         [1] Venus (geocéntrico)
 *         [2] Tierra (baricéntrico)
 *         [3] Marte (geocéntrico)
 *         [4] Júpiter (geocéntrico)
 *         [5] Saturno (geocéntrico)
 *         [6] Urano (geocéntrico)
 *         [7] Neptuno (geocéntrico)
 *         [8] Plutón (geocéntrico)
 *         [9] Luna (geocéntrico)
 *         [10] Sol (geocéntrico)
 */
tuple<Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix> JPL_Eph_DE430(double Mjd_TDB);

#endif