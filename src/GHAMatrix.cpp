#include "..\include\GHAMatrix.hpp"

/**
 * Calcula la matriz de rotación del ángulo horario de Greenwich (GHA), que transforma
 * coordenadas del sistema inercial (ECI) al sistema fijo a la Tierra (ECEF).
 * @param Mjd_UT1 Fecha juliana modificada (UT1)
 * @return Matriz de rotación 3x3 que representa la transformación ECI -> ECEF
 */
Matrix GHAMatrix(double Mjd_UT1) {
    double gast_value = gast(Mjd_UT1);
    Matrix GHAmat = R_z(gast_value);
    return GHAmat;
}