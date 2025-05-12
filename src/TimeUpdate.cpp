#include "..\include\TimeUpdate.hpp"

void TimeUpdate(Matrix& P, Matrix& Phi, Matrix& Qdt) {
	P = Phi * P * transpose(Phi) + Qdt;
}