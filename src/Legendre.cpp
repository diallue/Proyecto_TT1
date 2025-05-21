#include "..\include\Legendre.hpp"

using namespace std;

/**
 * Calcula los polinomios de Legendre asociados normalizados y sus derivadas.
 * 
 * @param n Grado máximo del polinomio
 * @param m Orden máximo del polinomio
 * @param fi Latitud geocéntrica en radianes
 * @param[out] pnm Matriz (n+1)x(m+1) de polinomios de Legendre normalizados
 * @param[out] dpnm Matriz (n+1)x(m+1) de derivadas de los polinomios
 */
void Legendre(int n, int m, double fi, Matrix &pnm, Matrix &dpnm) {
    pnm = Matrix(n+1, m+1);
    dpnm = Matrix(n+1, m+1);

    pnm(1, 1) = 1.0;
    dpnm(1, 1) = 0.0;
    
    if (m >= 1) {
        pnm(2, 2) = sqrt(3.0) * cos(fi);
        dpnm(2, 2) = -sqrt(3.0) * sin(fi);
    }

    for (int i = 2; i <= n; i++) {
        if (i <= m) {
            pnm(i+1, i+1) = sqrt((2.0*i + 1.0)/(2.0*i)) * cos(fi) * pnm(i, i);
            dpnm(i+1, i+1) = sqrt((2.0*i + 1.0)/(2.0*i)) * 
                            (cos(fi) * dpnm(i, i) - sin(fi) * pnm(i, i));
        }
    }

    for (int i = 1; i <= n; i++) {
        if (i <= m) {
            pnm(i+1, i) = sqrt(2.0*i + 1.0) * sin(fi) * pnm(i, i);
            dpnm(i+1, i) = sqrt(2.0*i + 1.0) * 
                          (cos(fi) * pnm(i, i) + sin(fi) * dpnm(i, i));
        }
    }

    for (int j = 0; j <= m; j++) {
        for (int i = j+1; i <= n; i++) {
            if (j+1 > m) continue;
            
            double factor = sqrt((2.0*i + 1.0)/((i-j)*(i+j)));
            
            double term1 = sqrt(2.0*i - 1.0) * sin(fi) * pnm(i, j+1);
            double term2 = 0.0;
            if (i-1 >= j+1) {
                term2 = sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * pnm(i-1, j+1);
            }
            pnm(i+1, j+1) = factor * (term1 - term2);
            
            double dterm1 = sqrt(2.0*i - 1.0) * sin(fi) * dpnm(i, j+1);
            double dterm2 = sqrt(2.0*i - 1.0) * cos(fi) * pnm(i, j+1);
            double dterm3 = 0.0;
            if (i-1 >= j+1) {
                dterm3 = sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * dpnm(i-1, j+1);
            }
            dpnm(i+1, j+1) = factor * (dterm1 + dterm2 - dterm3);
        }
    }
}