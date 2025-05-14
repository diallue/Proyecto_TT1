#include "..\include\JPL_Eph_DE430.hpp"

static Matrix PC;

tuple<Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix> JPL_Eph_DE430(double Mjd_TDB) {
    double JD = Mjd_TDB + 2400000.5;

    int i = 0;
    for (int j = 1; j <= PC.n_row; ++j) {
        if (PC(j, 1) <= JD && JD <= PC(j, 2)) {
            i = j;
            break;
        }
    }
    Matrix PCtemp = PC.extract_row(i);

    double t1 = PCtemp(1, 1) - 2400000.5;
    double dt = Mjd_TDB - t1;

    Matrix Cx_Earth(1, 13), Cy_Earth(1, 13), Cz_Earth(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_Earth(1, k) = PCtemp(1, 231 + k - 1);
        Cy_Earth(1, k) = PCtemp(1, 244 + k - 1);
        Cz_Earth(1, k) = PCtemp(1, 257 + k - 1);
    }
    Matrix Cx(1, 13), Cy(1, 13), Cz(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx(1, k) = PCtemp(1, 270 + k - 1);
        Cy(1, k) = PCtemp(1, 283 + k - 1);
        Cz(1, k) = PCtemp(1, 296 + k - 1);
    }
    Cx_Earth = Cx_Earth.union_vector(Cx, true);
    Cy_Earth = Cy_Earth.union_vector(Cy, true);
    Cz_Earth = Cz_Earth.union_vector(Cz, true);
    int j = (0 <= dt && dt <= 16) ? 0 : 1;
    double Mjd0 = t1 + 16 * j;
    Matrix r_Earth = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16,
        Cx_Earth.extract_column(13 * j + 1, 13 * j + 13),
        Cy_Earth.extract_column(13 * j + 1, 13 * j + 13),
        Cz_Earth.extract_column(13 * j + 1, 13 * j + 13));

    Matrix Cx_Moon(1, 13), Cy_Moon(1, 13), Cz_Moon(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_Moon(1, k) = PCtemp(1, 441 + k - 1);
        Cy_Moon(1, k) = PCtemp(1, 454 + k - 1);
        Cz_Moon(1, k) = PCtemp(1, 467 + k - 1);
    }
    for (int i = 1; i <= 7; ++i) {
        Matrix Cx(1, 13), Cy(1, 13), Cz(1, 13);
        for (int k = 1; k <= 13; ++k) {
            Cx(1, k) = PCtemp(1, 441 + 39 * i + k - 1);
            Cy(1, k) = PCtemp(1, 454 + 39 * i + k - 1);
            Cz(1, k) = PCtemp(1, 467 + 39 * i + k - 1);
        }
        Cx_Moon = Cx_Moon.union_vector(Cx, true);
        Cy_Moon = Cy_Moon.union_vector(Cy, true);
        Cz_Moon = Cz_Moon.union_vector(Cz, true);
    }
    if (0 <= dt && dt <= 4) j = 0;
    else if (4 < dt && dt <= 8) j = 1;
    else if (8 < dt && dt <= 12) j = 2;
    else if (12 < dt && dt <= 16) j = 3;
    else if (16 < dt && dt <= 20) j = 4;
    else if (20 < dt && dt <= 24) j = 5;
    else if (24 < dt && dt <= 28) j = 6;
    else j = 7;
    Mjd0 = t1 + 4 * j;
    Matrix r_Moon = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4,
        Cx_Moon.extract_column(13 * j + 1, 13 * j + 13),
        Cy_Moon.extract_column(13 * j + 1, 13 * j + 13),
        Cz_Moon.extract_column(13 * j + 1, 13 * j + 13));

    Matrix Cx_Sun(1, 11), Cy_Sun(1, 11), Cz_Sun(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Sun(1, k) = PCtemp(1, 753 + k - 1);
        Cy_Sun(1, k) = PCtemp(1, 764 + k - 1);
        Cz_Sun(1, k) = PCtemp(1, 775 + k - 1);
    }
    Cx(1, 11); Cy(1, 11); Cz(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx(1, k) = PCtemp(1, 786 + k - 1);
        Cy(1, k) = PCtemp(1, 797 + k - 1);
        Cz(1, k) = PCtemp(1, 808 + k - 1);
    }
    Cx_Sun = Cx_Sun.union_vector(Cx, true);
    Cy_Sun = Cy_Sun.union_vector(Cy, true);
    Cz_Sun = Cz_Sun.union_vector(Cz, true);
    j = (0 <= dt && dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16 * j;
    Matrix r_Sun = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16,
        Cx_Sun.extract_column(11 * j + 1, 11 * j + 11),
        Cy_Sun.extract_column(11 * j + 1, 11 * j + 11),
        Cz_Sun.extract_column(11 * j + 1, 11 * j + 11));

    Matrix Cx_Mercury(1, 14), Cy_Mercury(1, 14), Cz_Mercury(1, 14);
    for (int k = 1; k <= 14; ++k) {
        Cx_Mercury(1, k) = PCtemp(1, 3 + k - 1);
        Cy_Mercury(1, k) = PCtemp(1, 17 + k - 1);
        Cz_Mercury(1, k) = PCtemp(1, 31 + k - 1);
    }
    for (int i = 1; i <= 3; ++i) {
        Matrix Cx(1, 14), Cy(1, 14), Cz(1, 14);
        for (int k = 1; k <= 14; ++k) {
            Cx(1, k) = PCtemp(1, 3 + 42 * i + k - 1);
            Cy(1, k) = PCtemp(1, 17 + 42 * i + k - 1);
            Cz(1, k) = PCtemp(1, 31 + 42 * i + k - 1);
        }
        Cx_Mercury = Cx_Mercury.union_vector(Cx, true);
        Cy_Mercury = Cy_Mercury.union_vector(Cy, true);
        Cz_Mercury = Cz_Mercury.union_vector(Cz, true);
    }
    if (0 <= dt && dt <= 8) j = 0;
    else if (8 < dt && dt <= 16) j = 1;
    else if (16 < dt && dt <= 24) j = 2;
    else j = 3;
    Mjd0 = t1 + 8 * j;
    Matrix r_Mercury = 1e3 * Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8,
        Cx_Mercury.extract_column(14 * j + 1, 14 * j + 14),
        Cy_Mercury.extract_column(14 * j + 1, 14 * j + 14),
        Cz_Mercury.extract_column(14 * j + 1, 14 * j + 14));

    Matrix Cx_Venus(1, 10), Cy_Venus(1, 10), Cz_Venus(1, 10);
    for (int k = 1; k <= 10; ++k) {
        Cx_Venus(1, k) = PCtemp(1, 171 + k - 1);
        Cy_Venus(1, k) = PCtemp(1, 181 + k - 1);
        Cz_Venus(1, k) = PCtemp(1, 191 + k - 1);
    }
    Cx(1, 10); Cy(1, 10); Cz(1, 10);
    for (int k = 1; k <= 10; ++k) {
        Cx(1, k) = PCtemp(1, 201 + k - 1);
        Cy(1, k) = PCtemp(1, 211 + k - 1);
        Cz(1, k) = PCtemp(1, 221 + k - 1);
    }
    Cx_Venus = Cx_Venus.union_vector(Cx, true);
    Cy_Venus = Cy_Venus.union_vector(Cy, true);
    Cz_Venus = Cz_Venus.union_vector(Cz, true);
    j = (0 <= dt && dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16 * j;
    Matrix r_Venus = 1e3 * Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16,
        Cx_Venus.extract_column(10 * j + 1, 10 * j + 10),
        Cy_Venus.extract_column(10 * j + 1, 10 * j + 10),
        Cz_Venus.extract_column(10 * j + 1, 10 * j + 10));

    Matrix Cx_Mars(1, 11), Cy_Mars(1, 11), Cz_Mars(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Mars(1, k) = PCtemp(1, 309 + k - 1);
        Cy_Mars(1, k) = PCtemp(1, 320 + k - 1);
        Cz_Mars(1, k) = PCtemp(1, 331 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Mars = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32,
        Cx_Mars.extract_column(11 * j + 1, 11 * j + 11),
        Cy_Mars.extract_column(11 * j + 1, 11 * j + 11),
        Cz_Mars.extract_column(11 * j + 1, 11 * j + 11));

    Matrix Cx_Jupiter(1, 8), Cy_Jupiter(1, 8), Cz_Jupiter(1, 8);
    for (int k = 1; k <= 8; ++k) {
        Cx_Jupiter(1, k) = PCtemp(1, 342 + k - 1);
        Cy_Jupiter(1, k) = PCtemp(1, 350 + k - 1);
        Cz_Jupiter(1, k) = PCtemp(1, 358 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Jupiter = 1e3 * Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32,
        Cx_Jupiter.extract_column(8 * j + 1, 8 * j + 8),
        Cy_Jupiter.extract_column(8 * j + 1, 8 * j + 8),
        Cz_Jupiter.extract_column(8 * j + 1, 8 * j + 8));

    Matrix Cx_Saturn(1, 7), Cy_Saturn(1, 7), Cz_Saturn(1, 7);
    for (int k = 1; k <= 7; ++k) {
        Cx_Saturn(1, k) = PCtemp(1, 366 + k - 1);
        Cy_Saturn(1, k) = PCtemp(1, 373 + k - 1);
        Cz_Saturn(1, k) = PCtemp(1, 380 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Saturn = 1e3 * Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32,
        Cx_Saturn.extract_column(7 * j + 1, 7 * j + 7),
        Cy_Saturn.extract_column(7 * j + 1, 7 * j + 7),
        Cz_Saturn.extract_column(7 * j + 1, 7 * j + 7));

    Matrix Cx_Uranus(1, 6), Cy_Uranus(1, 6), Cz_Uranus(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Uranus(1, k) = PCtemp(1, 387 + k - 1);
        Cy_Uranus(1, k) = PCtemp(1, 393 + k - 1);
        Cz_Uranus(1, k) = PCtemp(1, 399 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Uranus = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
        Cx_Uranus.extract_column(6 * j + 1, 6 * j + 6),
        Cy_Uranus.extract_column(6 * j + 1, 6 * j + 6),
        Cz_Uranus.extract_column(6 * j + 1, 6 * j + 6));

    Matrix Cx_Neptune(1, 6), Cy_Neptune(1, 6), Cz_Neptune(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Neptune(1, k) = PCtemp(1, 405 + k - 1);
        Cy_Neptune(1, k) = PCtemp(1, 411 + k - 1);
        Cz_Neptune(1, k) = PCtemp(1, 417 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Neptune = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
        Cx_Neptune.extract_column(6 * j + 1, 6 * j + 6),
        Cy_Neptune.extract_column(6 * j + 1, 6 * j + 6),
        Cz_Neptune.extract_column(6 * j + 1, 6 * j + 6));

    Matrix Cx_Pluto(1, 6), Cy_Pluto(1, 6), Cz_Pluto(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Pluto(1, k) = PCtemp(1, 423 + k - 1);
        Cy_Pluto(1, k) = PCtemp(1, 429 + k - 1);
        Cz_Pluto(1, k) = PCtemp(1, 435 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Pluto = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
        Cx_Pluto.extract_column(6 * j + 1, 6 * j + 6),
        Cy_Pluto.extract_column(6 * j + 1, 6 * j + 6),
        Cz_Pluto.extract_column(6 * j + 1, 6 * j + 6));

    const double EMRAT = 81.30056907419062;
    const double EMRAT1 = 1.0 / (1.0 + EMRAT);
    r_Earth = r_Earth - r_Moon * EMRAT1;
    r_Mercury = r_Mercury - r_Earth;
    r_Venus = r_Venus - r_Earth;
    r_Mars = r_Mars - r_Earth;
    r_Jupiter = r_Jupiter - r_Earth;
    r_Saturn = r_Saturn - r_Earth;
    r_Uranus = r_Uranus - r_Earth;
    r_Neptune = r_Neptune - r_Earth;
    r_Pluto = r_Pluto - r_Earth;
    r_Sun = r_Sun - r_Earth;

    return make_tuple(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);
}