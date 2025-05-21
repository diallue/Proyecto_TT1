#include "..\include\G_AccelHarmonic.hpp"

/**
 * Calcula numéricamente el gradiente de la aceleración debida al campo gravitatorio armónico.
 * 
 * @param r Vector de posición en coordenadas geocéntricas (ECEF), dimensión 3x1
 * @param U Matriz de coeficientes armónicos normalizados del campo gravitatorio
 * @param n_max Grado máximo del desarrollo armónico
 * @param m_max Orden máximo del desarrollo armónico
 * @return Matriz 3x3 del gradiente de la aceleración (jacobiano)
 */
Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    double d = 1.0;
    
    Matrix G = zeros(3, 3);
    
    Matrix dr(3, 1);
    for (int i = 1; i <= 3; i++) {
        dr(i,1) = 0.0;
    }
    
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dr(j,1) = 0.0;
        }
        dr(i,1) = d;
        
        Matrix dr_half = dr * (0.5);
        Matrix minus_dr_half = dr * (-0.5); 
        
        Matrix r_plus = r + dr_half;
        Matrix r_minus = r + minus_dr_half;
        
        Matrix accel_plus = AccelHarmonic(r_plus, U, n_max, m_max);
        Matrix accel_minus = AccelHarmonic(r_minus, U, n_max, m_max);
        Matrix da = accel_plus - accel_minus;
        
        for (int j = 1; j <= 3; j++) {
            G(j,i) = da(j,1) / d;
        }
    }
    
    return G;
}