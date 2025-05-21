#include "..\include\MeasUpdate.hpp"

/**
 * Realiza la actualización de medida de un filtro de Kalman.
 * 
 * @param x Estado estimado antes de la medida (vector columna de dimensión n)
 * @param z Vector de medida observado (vector columna de dimensión m)
 * @param g Vector de medida estimado (predicción de la medida)
 * @param s Vector con las desviaciones típicas de la medida (dimensión m×1)
 * @param G Matriz de observación (dimensión m×n)
 * @param P Matriz de covarianza del estado estimado (dimensión n×n)
 * @param n Dimensión del vector de estado
 * @return Tupla con:
 *         - K: matriz de ganancia de Kalman (n×m)
 *         - x_new: nuevo vector de estado estimado (n×1)
 *         - P_new: nueva matriz de covarianza del estado (n×n)
 */
tuple<Matrix, Matrix, Matrix> MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n) {
    int m = z.n_row;
    
    Matrix Inv_W = zeros(m, m);
    for (int i = 1; i <= m; i++) {
        Inv_W(i,i) = s(i,1) * s(i,1);
    }
    
    Matrix Gt = transpose(G);
    Matrix P_Gt = P * Gt;              
    Matrix G_P_Gt = G * P_Gt;          
    Matrix sum = Inv_W + G_P_Gt;       
    Matrix inv_sum = inv(sum);
    Matrix K = P_Gt * inv_sum;         
    
    Matrix innovation = z - g;
    Matrix K_innovation = K * innovation; 
    Matrix x_new = x + K_innovation;  
    
    Matrix I = Matrix::eye(n);
    Matrix KG = K * G;
    Matrix I_minus_KG = I - KG;         
    Matrix P_new = I_minus_KG * P;      
    
    return make_tuple(K, x_new, P_new);
}