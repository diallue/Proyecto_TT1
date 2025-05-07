#include "..\include\Legendre.hpp"

using namespace std;

tuple<Matrix, Matrix> Legendre(int n, int m, double fi) {
	Matri pnm = zeros(n+1, m+1);
	Matrix dpnm = zeros(n+1, m+1);
    
    pnm(2,2) = sqrt(3) * cos(fi);
    dpnm(2,2) = sqrt(3) * sin(fi);
	
    for (int i = 2; i <= n; i++) {
        result.pnm[i][i] = sqrt((2.0*i + 1.0)/(2.0*i)) * cos(fi) * result.pnm[i-1][i-1];
    }
    
    for (int i = 2; i <= n; i++) {
        result.dpnm[i][i] = sqrt((2.0*i + 1.0)/(2.0*i)) * 
                           (cos(fi) * result.dpnm[i-1][i-1] - 
                           sin(fi) * result.pnm[i-1][i-1]);
    }
    
    for (int i = 1; i <= n; i++) {
        result.pnm[i][i-1] = sqrt(2.0*i + 1.0) * sin(fi) * result.pnm[i-1][i-1];
    }
    
    for (int i = 1; i <= n; i++) {
        result.dpnm[i][i-1] = sqrt(2.0*i + 1.0) * 
                             (cos(fi) * result.pnm[i-1][i-1] + 
                              sin(fi) * result.dpnm[i-1][i-1]);
    }
    
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            result.pnm[i][j] = sqrt((2.0*i + 1.0)/((i-j)*(i+j))) * 
                              (sqrt(2.0*i - 1.0) * sin(fi) * result.pnm[i-1][j] - 
                              sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * result.pnm[i-2][j]);
        }
        j++;
        k++;
        if (j > m) break;
    }
    
    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            result.dpnm[i][j] = sqrt((2.0*i + 1.0)/((i-j)*(i+j))) * 
                               (sqrt(2.0*i - 1.0) * sin(fi) * result.dpnm[i-1][j] + 
                                sqrt(2.0*i - 1.0) * cos(fi) * result.pnm[i-1][j] - 
                                sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * result.dpnm[i-2][j]);
        }
        j++;
        k++;
        if (j > m) break;
    }
    
    return make_tuple(pnm, dpnm);
}