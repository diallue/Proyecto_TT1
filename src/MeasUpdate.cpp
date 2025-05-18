#include "..\include\MeasUpdate.hpp"

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