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
    Matrix Cx_earth(1, 13), Cy_earth(1, 13), Cz_earth(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_earth(1, k) = PCtemp(1, 270 + k - 1);
        Cy_earth(1, k) = PCtemp(1, 283 + k - 1);
        Cz_earth(1, k) = PCtemp(1, 296 + k - 1);
    }
    Cx_Earth = Cx_Earth.union_vector(Cx_earth, true);
    Cy_Earth = Cy_Earth.union_vector(Cy_earth, true);
    Cz_Earth = Cz_Earth.union_vector(Cz_earth, true);
    int j = (0 <= dt && dt <= 16) ? 0 : 1;
    double Mjd0 = t1 + 16 * j;
    Matrix Cx_Earth_sub(1, 13), Cy_Earth_sub(1, 13), Cz_Earth_sub(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_Earth_sub(1, k) = Cx_Earth(1, 13 * j + k);
        Cy_Earth_sub(1, k) = Cy_Earth(1, 13 * j + k);
        Cz_Earth_sub(1, k) = Cz_Earth(1, 13 * j + k);
    }
    Matrix r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, Cx_Earth_sub, Cy_Earth_sub, Cz_Earth_sub) * 1e3;

    Matrix Cx_Moon(1, 13), Cy_Moon(1, 13), Cz_Moon(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_Moon(1, k) = PCtemp(1, 441 + k - 1);
        Cy_Moon(1, k) = PCtemp(1, 454 + k - 1);
        Cz_Moon(1, k) = PCtemp(1, 467 + k - 1);
    }
    for (int i = 1; i <= 7; ++i) {
        Matrix Cx_moon(1, 13), Cy_moon(1, 13), Cz_moon(1, 13);
        for (int k = 1; k <= 13; ++k) {
            Cx_moon(1, k) = PCtemp(1, 441 + 39 * i + k - 1);
            Cy_moon(1, k) = PCtemp(1, 454 + 39 * i + k - 1);
            Cz_moon(1, k) = PCtemp(1, 467 + 39 * i + k - 1);
        }
        Cx_Moon = Cx_Moon.union_vector(Cx_moon, true);
        Cy_Moon = Cy_Moon.union_vector(Cy_moon, true);
        Cz_Moon = Cz_Moon.union_vector(Cz_moon, true);
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
    Matrix Cx_Moon_sub(1, 13), Cy_Moon_sub(1, 13), Cz_Moon_sub(1, 13);
    for (int k = 1; k <= 13; ++k) {
        Cx_Moon_sub(1, k) = Cx_Moon(1, 13 * j + k);
        Cy_Moon_sub(1, k) = Cy_Moon(1, 13 * j + k);
        Cz_Moon_sub(1, k) = Cz_Moon(1, 13 * j + k);
    }
    Matrix r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, Cx_Moon_sub, Cy_Moon_sub, Cz_Moon_sub) * 1e3;

    Matrix Cx_Sun(1, 11), Cy_Sun(1, 11), Cz_Sun(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Sun(1, k) = PCtemp(1, 753 + k - 1);
        Cy_Sun(1, k) = PCtemp(1, 764 + k - 1);
        Cz_Sun(1, k) = PCtemp(1, 775 + k - 1);
    }
    Matrix Cx_sun(1, 11), Cy_sun(1, 11), Cz_sun(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_sun(1, k) = PCtemp(1, 786 + k - 1);
        Cy_sun(1, k) = PCtemp(1, 797 + k - 1);
        Cz_sun(1, k) = PCtemp(1, 808 + k - 1);
    }
    Cx_Sun = Cx_Sun.union_vector(Cx_sun, true);
    Cy_Sun = Cy_Sun.union_vector(Cy_sun, true);
    Cz_Sun = Cz_Sun.union_vector(Cz_sun, true);
    j = (0 <= dt && dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16 * j;
    Matrix Cx_Sun_sub(1, 11), Cy_Sun_sub(1, 11), Cz_Sun_sub(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Sun_sub(1, k) = Cx_Sun(1, 11 * j + k);
        Cy_Sun_sub(1, k) = Cy_Sun(1, 11 * j + k);
        Cz_Sun_sub(1, k) = Cz_Sun(1, 11 * j + k);
    }
    Matrix r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, Cx_Sun_sub, Cy_Sun_sub, Cz_Sun_sub) * 1e3;

    Matrix Cx_Mercury(1, 14), Cy_Mercury(1, 14), Cz_Mercury(1, 14);
    for (int k = 1; k <= 14; ++k) {
        Cx_Mercury(1, k) = PCtemp(1, 3 + k - 1);
        Cy_Mercury(1, k) = PCtemp(1, 17 + k - 1);
        Cz_Mercury(1, k) = PCtemp(1, 31 + k - 1);
    }
    for (int i = 1; i <= 3; ++i) {
        Matrix Cx_mercury(1, 14), Cy_mercury(1, 14), Cz_mercury(1, 14);
        for (int k = 1; k <= 14; ++k) {
            Cx_mercury(1, k) = PCtemp(1, 3 + 42 * i + k - 1);
            Cy_mercury(1, k) = PCtemp(1, 17 + 42 * i + k - 1);
            Cz_mercury(1, k) = PCtemp(1, 31 + 42 * i + k - 1);
        }
        Cx_Mercury = Cx_Mercury.union_vector(Cx_mercury, true);
        Cy_Mercury = Cy_Mercury.union_vector(Cy_mercury, true);
        Cz_Mercury = Cz_Mercury.union_vector(Cz_mercury, true);
    }
    if (0 <= dt && dt <= 8) j = 0;
    else if (8 < dt && dt <= 16) j = 1;
    else if (16 < dt && dt <= 24) j = 2;
    else j = 3;
    Mjd0 = t1 + 8 * j;
    Matrix Cx_Mercury_sub(1, 14), Cy_Mercury_sub(1, 14), Cz_Mercury_sub(1, 14);
    for (int k = 1; k <= 14; ++k) {
        Cx_Mercury_sub(1, k) = Cx_Mercury(1, 14 * j + k);
        Cy_Mercury_sub(1, k) = Cy_Mercury(1, 14 * j + k);
        Cz_Mercury_sub(1, k) = Cz_Mercury(1, 14 * j + k);
    }
    Matrix r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, Cx_Mercury_sub, Cy_Mercury_sub, Cz_Mercury_sub) * 1e3;

    Matrix Cx_Venus(1, 10), Cy_Venus(1, 10), Cz_Venus(1, 10);
    for (int k = 1; k <= 10; ++k) {
        Cx_Venus(1, k) = PCtemp(1, 171 + k - 1);
        Cy_Venus(1, k) = PCtemp(1, 181 + k - 1);
        Cz_Venus(1, k) = PCtemp(1, 191 + k - 1);
    }
    Matrix Cx_venus(1, 10), Cy_venus(1, 10), Cz_venus(1, 10);
    for (int k = 1; k <= 10; ++k) {
        Cx_venus(1, k) = PCtemp(1, 201 + k - 1);
        Cy_venus(1, k) = PCtemp(1, 211 + k - 1);
        Cz_venus(1, k) = PCtemp(1, 221 + k - 1);
    }
    Cx_Venus = Cx_Venus.union_vector(Cx_venus, true);
    Cy_Venus = Cy_Venus.union_vector(Cy_venus, true);
    Cz_Venus = Cz_Venus.union_vector(Cz_venus, true);
    j = (0 <= dt && dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16 * j;
    Matrix Cx_Venus_sub(1, 10), Cy_Venus_sub(1, 10), Cz_Venus_sub(1, 10);
    for (int k = 1; k <= 10; ++k) {
        Cx_Venus_sub(1, k) = Cx_Venus(1, 10 * j + k);
        Cy_Venus_sub(1, k) = Cy_Venus(1, 10 * j + k);
        Cz_Venus_sub(1, k) = Cz_Venus(1, 10 * j + k);
    }
    Matrix r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, Cx_Venus_sub, Cy_Venus_sub, Cz_Venus_sub) * 1e3;

    Matrix Cx_Mars(1, 11), Cy_Mars(1, 11), Cz_Mars(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Mars(1, k) = PCtemp(1, 309 + k - 1);
        Cy_Mars(1, k) = PCtemp(1, 320 + k - 1);
        Cz_Mars(1, k) = PCtemp(1, 331 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Mars_sub(1, 11), Cy_Mars_sub(1, 11), Cz_Mars_sub(1, 11);
    for (int k = 1; k <= 11; ++k) {
        Cx_Mars_sub(1, k) = Cx_Mars(1, 11 * j + k);
        Cy_Mars_sub(1, k) = Cy_Mars(1, 11 * j + k);
        Cz_Mars_sub(1, k) = Cz_Mars(1, 11 * j + k);
    }
    Matrix r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, Cx_Mars_sub, Cy_Mars_sub, Cz_Mars_sub) * 1e3;

    Matrix Cx_Jupiter(1, 8), Cy_Jupiter(1, 8), Cz_Jupiter(1, 8);
    for (int k = 1; k <= 8; ++k) {
        Cx_Jupiter(1, k) = PCtemp(1, 342 + k - 1);
        Cy_Jupiter(1, k) = PCtemp(1, 350 + k - 1);
        Cz_Jupiter(1, k) = PCtemp(1, 358 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Jupiter_sub(1, 8), Cy_Jupiter_sub(1, 8), Cz_Jupiter_sub(1, 8);
    for (int k = 1; k <= 8; ++k) {
        Cx_Jupiter_sub(1, k) = Cx_Jupiter(1, 8 * j + k);
        Cy_Jupiter_sub(1, k) = Cx_Jupiter(1, 8 * j + k);
        Cz_Jupiter_sub(1, k) = Cz_Jupiter(1, 8 * j + k);
    }
    Matrix r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, Cx_Jupiter_sub, Cy_Jupiter_sub, Cz_Jupiter_sub) * 1e3;

    Matrix Cx_Saturn(1, 7), Cy_Saturn(1, 7), Cz_Saturn(1, 7);
    for (int k = 1; k <= 7; ++k) {
        Cx_Saturn(1, k) = PCtemp(1, 366 + k - 1);
        Cy_Saturn(1, k) = PCtemp(1, 373 + k - 1);
        Cz_Saturn(1, k) = PCtemp(1, 380 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Saturn_sub(1, 7), Cy_Saturn_sub(1, 7), Cz_Saturn_sub(1, 7);
    for (int k = 1; k <= 7; ++k) {
        Cx_Saturn_sub(1, k) = Cx_Saturn(1, 7 * j + k);
        Cy_Saturn_sub(1, k) = Cx_Saturn(1, 7 * j + k);
        Cz_Saturn_sub(1, k) = Cz_Saturn(1, 7 * j + k);
    }
    Matrix r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, Cx_Saturn_sub, Cy_Saturn_sub, Cz_Saturn_sub) * 1e3;

    Matrix Cx_Uranus(1, 6), Cy_Uranus(1, 6), Cz_Uranus(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Uranus(1, k) = PCtemp(1, 387 + k - 1);
        Cy_Uranus(1, k) = PCtemp(1, 393 + k - 1);
        Cz_Uranus(1, k) = PCtemp(1, 399 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Uranus_sub(1, 6), Cy_Uranus_sub(1, 6), Cz_Uranus_sub(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Uranus_sub(1, k) = Cx_Uranus(1, 6 * j + k);
        Cy_Uranus_sub(1, k) = Cx_Uranus(1, 6 * j + k);
        Cz_Uranus_sub(1, k) = Cz_Uranus(1, 6 * j + k);
    }
    Matrix r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Uranus_sub, Cy_Uranus_sub, Cz_Uranus_sub) * 1e3;

    Matrix Cx_Neptune(1, 6), Cy_Neptune(1, 6), Cz_Neptune(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Neptune(1, k) = PCtemp(1, 405 + k - 1);
        Cy_Neptune(1, k) = PCtemp(1, 411 + k - 1);
        Cz_Neptune(1, k) = PCtemp(1, 417 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Neptune_sub(1, 6), Cy_Neptune_sub(1, 6), Cz_Neptune_sub(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Neptune_sub(1, k) = Cx_Neptune(1, 6 * j + k);
        Cy_Neptune_sub(1, k) = Cx_Neptune(1, 6 * j + k);
        Cz_Neptune_sub(1, k) = Cz_Neptune(1, 6 * j + k);
    }
    Matrix r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Neptune_sub, Cy_Neptune_sub, Cz_Neptune_sub) * 1e3;

    Matrix Cx_Pluto(1, 6), Cy_Pluto(1, 6), Cz_Pluto(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Pluto(1, k) = PCtemp(1, 423 + k - 1);
        Cy_Pluto(1, k) = PCtemp(1, 429 + k - 1);
        Cz_Pluto(1, k) = PCtemp(1, 435 + k - 1);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Pluto_sub(1, 6), Cy_Pluto_sub(1, 6), Cz_Pluto_sub(1, 6);
    for (int k = 1; k <= 6; ++k) {
        Cx_Pluto_sub(1, k) = Cx_Pluto(1, 6 * j + k);
        Cy_Pluto_sub(1, k) = Cx_Pluto(1, 6 * j + k);
        Cz_Pluto_sub(1, k) = Cz_Pluto(1, 6 * j + k);
    }
    Matrix r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Pluto_sub, Cy_Pluto_sub, Cz_Pluto_sub) * 1e3;

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