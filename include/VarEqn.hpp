#ifndef _VAREQN_
#define _VAREQN_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\global.hpp"
#include <vector>
#include <tuple>

/**
 * @details Calcula la ecuación variacional para la propagación de órbitas y matrices de transición de estado
 * @param x Tiempo desde la época de referencia [s]
 * @param yPhi Vector de estado extendido (42x1) que contiene:
 *             - Posición (3 elementos)
 *             - Velocidad (3 elementos)
 *             - Matriz de transición de estado (36 elementos)
 * @return Derivada del vector de estado extendido (42x1)
 */
Matrix VarEqn(double x, Matrix yPhi);

#endif