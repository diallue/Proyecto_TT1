#include "..\include\Legendre.hpp"

using namespace std;

tuple<Matrix, Matrix> Legendre(double n, double m, double fi) {
    int n_int = static_cast<int>(n);
    int m_int = static_cast<int>(m);
    
    if (n_int < m_int || n_int < 0 || m_int < 0) {
        throw std::invalid_argument("n must be >= m and both must be non-negative");
    }

    Matrix pnm = zeros(n_int+1, m_int+1);
    Matrix dpnm = zeros(n_int+1, m_int+1);
    
    pnm(0,0) = 1.0;
    dpnm(0,0) = 0.0;
    if (n_int >= 1 && m_int >= 1) {
        pnm(1,1) = sqrt(3.0) * cos(fi);
        dpnm(1,1) = -sqrt(3.0) * sin(fi);
    }
    
    for (int i = 2; i <= n_int; i++) {
        pnm(i,i) = sqrt((2.0*i + 1.0)/(2.0*i)) * cos(fi) * pnm(i-1,i-1);
    }
    
    for (int i = 2; i <= n_int; i++) {
        dpnm(i,i) = sqrt((2.0*i + 1.0)/(2.0*i)) * 
                    (cos(fi) * dpnm(i-1,i-1) - sin(fi) * pnm(i-1,i-1));
    }
    
    for (int i = 1; i <= n_int; i++) {
        pnm(i,i-1) = sqrt(2.0*i + 1.0) * sin(fi) * pnm(i-1,i-1);
    }
    
    for (int i = 1; i <= n_int; i++) {
        dpnm(i,i-1) = sqrt(2.0*i + 1.0) * 
                      (cos(fi) * pnm(i-1,i-1) + sin(fi) * dpnm(i-1,i-1));
    }
    
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n_int; i++) {
            pnm(i,j) = sqrt((2.0*i + 1.0)/((i-j)*(i+j))) * 
                       (sqrt(2.0*i - 1.0) * sin(fi) * pnm(i-1,j) - 
                        sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * pnm(i-2,j));
        }
        j++;
        k++;
        if (j > m_int) break;
    }
    
    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n_int; i++) {
            dpnm(i,j) = sqrt((2.0*i + 1.0)/((i-j)*(i+j))) * 
                        (sqrt(2.0*i - 1.0) * sin(fi) * dpnm(i-1,j) + 
                         sqrt(2.0*i - 1.0) * cos(fi) * pnm(i-1,j) - 
                         sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * dpnm(i-2,j));
        }
        j++;
        k++;
        if (j > m_int) break;
    }
    
    return make_tuple(pnm, dpnm);
}