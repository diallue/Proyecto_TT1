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

Matrix VarEqn(double x, Matrix& yPhi);

#endif