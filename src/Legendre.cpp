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
std::tuple<Matrix, Matrix> Legendre(int n, int m, double fi) {
    Matrix pnm = zeros(n + 1, m + 1);
    Matrix dpnm = zeros(n + 1, m + 1);

    pnm(1, 1) = 1.0;
    dpnm(1, 1) = 0.0;
    if (n >= 1 && m >= 1) {
        pnm(2, 2) = std::sqrt(3.0) * std::cos(fi);
        dpnm(2, 2) = -std::sqrt(3.0) * std::sin(fi);
    }

    for (int i = 2; i <= n; ++i) {
        pnm(i + 1, i + 1) = std::sqrt((2.0 * i + 1.0) / (2.0 * i)) * std::cos(fi) * pnm(i, i);
    }
    for (int i = 2; i <= n; ++i) {
        dpnm(i + 1, i + 1) = std::sqrt((2.0 * i + 1.0) / (2.0 * i)) * (std::cos(fi) * dpnm(i, i) - std::sin(fi) * pnm(i, i));
    }

    for (int i = 1; i <= n; ++i) {
        pnm(i + 1, i) = std::sqrt(2.0 * i + 1.0) * std::sin(fi) * pnm(i, i);
    }
    for (int i = 1; i <= n; ++i) {
        dpnm(i + 1, i) = std::sqrt(2.0 * i + 1.0) * (std::cos(fi) * pnm(i, i) + std::sin(fi) * dpnm(i, i));
    }

    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; ++i) {
            pnm(i + 1, j + 1) = std::sqrt((2.0 * i + 1.0) / ((i - j) * (i + j))) * (std::sqrt(2.0 * i - 1.0) * std::sin(fi) * pnm(i, j + 1) -std::sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm(i - 1, j + 1));
        }
        j = j + 1;
        k = k + 1;
        if (j > m) {
            break;
        }
    }

    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; ++i) {
            dpnm(i + 1, j + 1) = std::sqrt((2.0 * i + 1.0) / ((i - j) * (i + j))) * (
                std::sqrt(2.0 * i - 1.0) * std::sin(fi) * dpnm(i, j + 1) +
                std::sqrt(2.0 * i - 1.0) * std::cos(fi) * pnm(i, j + 1) -
                std::sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * dpnm(i - 1, j + 1)
            );
        }
        j = j + 1;
        k = k + 1;
        if (j > m) {
            break;
        }
    }

    return std::make_tuple(pnm, dpnm);
}