#indef _AZELPA_
#define _AZELPA_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>
#include <vector>

tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s);

#endif