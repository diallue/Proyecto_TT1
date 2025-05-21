#include "..\include\LTC.hpp"

/**
 * Calcula la matriz de transformación Local Tangent Coordinate (LTC) en función
 * de la latitud y longitud geodésica del observador sobre la Tierra.
 *
 * @param lon Longitud geodésica en radianes
 * @param lat Latitud geodésica en radianes
 * @return Matriz de transformación 3x3 del sistema ECEF al sistema local ENZ
 */
Matrix LTC(double lon, double lat) {
    Matrix M = R_y(-1.0 * lat) * R_z(lon);

    for (int j = 1; j <= 3; ++j) {
        double Aux = M(1, j);
        M(1, j) = M(2, j);
        M(2, j) = M(3, j);
        M(3, j) = Aux;
    }

    return M;
}