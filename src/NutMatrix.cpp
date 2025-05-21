#include "..\include\NutMatrix.hpp"

/**
 * Calcula la matriz de nutación que transforma coordenadas del sistema medio al sistema verdadero.
 *
 * @param Mjd_TT Fecha en Tiempo Terrestre (Modified Julian Date)
 * @return Matriz de transformación 3x3 que representa la nutación
 */
Matrix NutMatrix(double Mjd_TT) {
	double eps = MeanObliquity(Mjd_TT);

    auto [dpsi, deps] = NutAngles(Mjd_TT);

    Matrix NutMat = R_x(-eps - deps) * R_z(-dpsi) * R_x(eps);

    return NutMat;
}