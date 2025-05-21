#include "..\include\TimeUpdate.hpp"

/**
 * Realiza la actualización temporal de la matriz de covarianza en un filtro de Kalman.
 * 
 * @param P Matriz de covarianza del estado (n x n)
 * @param Phi Matriz de transición de estado (n x n)
 * @param Qdt Matriz de covarianza del ruido del proceso integrada (n x n)
 */
void TimeUpdate(Matrix& P, Matrix& Phi, Matrix& Qdt) {
	P = Phi * P * transpose(Phi) + Qdt;
}