#include "..\include\JPL_Eph_DE430.hpp"

tuple<Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix> JPL_Eph_DE430(double Mjd_TDB) {
    const double EMRAT = 81.30056907419062;
    const double EMRAT1 = 1.0 / (1.0 + EMRAT);

    double JD = Mjd_TDB + 2400000.5;
    int i = 0;
    for (int k = 1; k <= PC.n_row; ++k) {
        if (PC(k, 1) <= JD && JD <= PC(k, 2)) {
            i = k;
            break;
        }
    }
    Matrix PCtemp = PC.extract_row(i);
	std::cout << "PCtemp: " << PCtemp.n_row << "x" << PCtemp.n_column << std::endl;

    double t1 = PCtemp(1, 1) - 2400000.5;
    double dt = Mjd_TDB - t1;

    Matrix Cx_Earth(26, 1), Cy_Earth(26, 1), Cz_Earth(26, 1);
    for (int k = 0; k < 13; ++k) {
        Cx_Earth(k + 1, 1) = PCtemp(1, 231 + k);
        Cy_Earth(k + 1, 1) = PCtemp(1, 244 + k);
        Cz_Earth(k + 1, 1) = PCtemp(1, 257 + k);
    }
    for (int k = 0; k < 13; ++k) {
        Cx_Earth(13 + k + 1, 1) = PCtemp(1, 270 + k);
        Cy_Earth(13 + k + 1, 1) = PCtemp(1, 283 + k);
        Cz_Earth(13 + k + 1, 1) = PCtemp(1, 296 + k);
    }

    int j = (dt <= 16.0) ? 0 : 1;
    double Mjd0 = t1 + 16.0 * j;
    Matrix Cx_Earth_sub(13, 1), Cy_Earth_sub(13, 1), Cz_Earth_sub(13, 1);
    for (int k = 0; k < 13; ++k) {
        Cx_Earth_sub(k + 1, 1) = Cx_Earth(j * 13 + k + 1, 1);
        Cy_Earth_sub(k + 1, 1) = Cy_Earth(j * 13 + k + 1, 1);
        Cz_Earth_sub(k + 1, 1) = Cz_Earth(j * 13 + k + 1, 1);
    }
    Matrix r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16.0, 
                            Cx_Earth_sub, Cy_Earth_sub, Cz_Earth_sub) * 1e3;

    Matrix Cx_Moon(104, 1), Cy_Moon(104, 1), Cz_Moon(104, 1);
    for (int k = 0; k < 13; ++k) {
        Cx_Moon(k + 1, 1) = PCtemp(1, 441 + k);
        Cy_Moon(k + 1, 1) = PCtemp(1, 454 + k);
        Cz_Moon(k + 1, 1) = PCtemp(1, 467 + k);
    }
    for (int m = 1; m <= 7; ++m) {
        for (int k = 0; k < 13; ++k) {
            Cx_Moon(13 * m + k + 1, 1) = PCtemp(1, 441 + 39 * m + k);
            Cy_Moon(13 * m + k + 1, 1) = PCtemp(1, 454 + 39 * m + k);
            Cz_Moon(13 * m + k + 1, 1) = PCtemp(1, 467 + 39 * m + k);
        }
    }
    j = (dt <= 4.0) ? 0 : (dt <= 8.0) ? 1 : (dt <= 12.0) ? 2 : (dt <= 16.0) ? 3 :
        (dt <= 20.0) ? 4 : (dt <= 24.0) ? 5 : (dt <= 28.0) ? 6 : 7;
    Mjd0 = t1 + 4.0 * j;
    Matrix Cx_Moon_sub(13, 1), Cy_Moon_sub(13, 1), Cz_Moon_sub(13, 1);
    for (int k = 0; k < 13; ++k) {
        Cx_Moon_sub(k + 1, 1) = Cx_Moon(j * 13 + k + 1, 1);
        Cy_Moon_sub(k + 1, 1) = Cy_Moon(j * 13 + k + 1, 1);
        Cz_Moon_sub(k + 1, 1) = Cz_Moon(j * 13 + k + 1, 1);
    }
    Matrix r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4.0, 
                           Cx_Moon_sub, Cy_Moon_sub, Cz_Moon_sub) * 1e3;

    Matrix Cx_Sun(22, 1), Cy_Sun(22, 1), Cz_Sun(22, 1);
    for (int k = 0; k < 11; ++k) {
        Cx_Sun(k + 1, 1) = PCtemp(1, 753 + k);
        Cy_Sun(k + 1, 1) = PCtemp(1, 764 + k);
        Cz_Sun(k + 1, 1) = PCtemp(1, 775 + k);
    }
    for (int k = 0; k < 11; ++k) {
        Cx_Sun(11 + k + 1, 1) = PCtemp(1, 786 + k);
        Cy_Sun(11 + k + 1, 1) = PCtemp(1, 797 + k);
        Cz_Sun(11 + k + 1, 1) = PCtemp(1, 808 + k);
    }
    j = (dt <= 16.0) ? 0 : 1;
    Mjd0 = t1 + 16.0 * j;
    Matrix Cx_Sun_sub(11, 1), Cy_Sun_sub(11, 1), Cz_Sun_sub(11, 1);
    for (int k = 0; k < 11; ++k) {
        Cx_Sun_sub(k + 1, 1) = Cx_Sun(j * 11 + k + 1, 1);
        Cy_Sun_sub(k + 1, 1) = Cy_Sun(j * 11 + k + 1, 1);
        Cz_Sun_sub(k + 1, 1) = Cz_Sun(j * 11 + k + 1, 1);
    }
    Matrix r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16.0, 
                          Cx_Sun_sub, Cy_Sun_sub, Cz_Sun_sub) * 1e3;

    Matrix Cx_Mercury(56, 1), Cy_Mercury(56, 1), Cz_Mercury(56, 1);
    for (int k = 0; k < 14; ++k) {
        Cx_Mercury(k + 1, 1) = PCtemp(1, 3 + k);
        Cy_Mercury(k + 1, 1) = PCtemp(1, 17 + k);
        Cz_Mercury(k + 1, 1) = PCtemp(1, 31 + k);
    }
    for (int m = 1; m <= 3; ++m) {
        for (int k = 0; k < 14; ++k) {
            Cx_Mercury(14 * m + k + 1, 1) = PCtemp(1, 3 + 42 * m + k);
            Cy_Mercury(14 * m + k + 1, 1) = PCtemp(1, 17 + 42 * m + k);
            Cz_Mercury(14 * m + k + 1, 1) = PCtemp(1, 31 + 42 * m + k);
        }
    }
    j = (dt <= 8.0) ? 0 : (dt <= 16.0) ? 1 : (dt <= 24.0) ? 2 : 3;
    Mjd0 = t1 + 8.0 * j;
    Matrix Cx_Mercury_sub(14, 1), Cy_Mercury_sub(14, 1), Cz_Mercury_sub(14, 1);
    for (int k = 0; k < 14; ++k) {
        Cx_Mercury_sub(k + 1, 1) = Cx_Mercury(j * 14 + k + 1, 1);
        Cy_Mercury_sub(k + 1, 1) = Cy_Mercury(j * 14 + k + 1, 1);
        Cz_Mercury_sub(k + 1, 1) = Cz_Mercury(j * 14 + k + 1, 1);
    }
    Matrix r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8.0, 
                              Cx_Mercury_sub, Cy_Mercury_sub, Cz_Mercury_sub) * 1e3;

    Matrix Cx_Venus(20, 1), Cy_Venus(20, 1), Cz_Venus(20, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Venus(k + 1, 1) = PCtemp(1, 171 + k);
        Cy_Venus(k + 1, 1) = PCtemp(1, 181 + k);
        Cz_Venus(k + 1, 1) = PCtemp(1, 191 + k);
    }
    for (int k = 0; k < 10; ++k) {
        Cx_Venus(10 + k + 1, 1) = PCtemp(1, 201 + k);
        Cy_Venus(10 + k + 1, 1) = PCtemp(1, 211 + k);
        Cz_Venus(10 + k + 1, 1) = PCtemp(1, 221 + k);
    }
    j = (dt <= 16.0) ? 0 : 1;
    Mjd0 = t1 + 16.0 * j;
    Matrix Cx_Venus_sub(10, 1), Cy_Venus_sub(10, 1), Cz_Venus_sub(10, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Venus_sub(k + 1, 1) = Cx_Venus(j * 10 + k + 1, 1);
        Cy_Venus_sub(k + 1, 1) = Cy_Venus(j * 10 + k + 1, 1);
        Cz_Venus_sub(k + 1, 1) = Cz_Venus(j * 10 + k + 1, 1);
    }
    Matrix r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16.0, 
                            Cx_Venus_sub, Cy_Venus_sub, Cz_Venus_sub) * 1e3;

    Matrix Cx_Mars(11, 1), Cy_Mars(11, 1), Cz_Mars(11, 1);
    for (int k = 0; k < 11; ++k) {
        Cx_Mars(k + 1, 1) = PCtemp(1, 309 + k);
        Cy_Mars(k + 1, 1) = PCtemp(1, 320 + k);
        Cz_Mars(k + 1, 1) = PCtemp(1, 331 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32.0, Cx_Mars, Cy_Mars, Cz_Mars) * 1e3;

    Matrix Cx_Jupiter(8, 1), Cy_Jupiter(8, 1), Cz_Jupiter(8, 1);
    for (int k = 0; k < 8; ++k) {
        Cx_Jupiter(k + 1, 1) = PCtemp(1, 342 + k);
        Cy_Jupiter(k + 1, 1) = PCtemp(1, 350 + k);
        Cz_Jupiter(k + 1, 1) = PCtemp(1, 358 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32.0, Cx_Jupiter, Cy_Jupiter, Cz_Jupiter) * 1e3;

    Matrix Cx_Saturn(7, 1), Cy_Saturn(7, 1), Cz_Saturn(7, 1);
    for (int k = 0; k < 7; ++k) {
        Cx_Saturn(k + 1, 1) = PCtemp(1, 366 + k);
        Cy_Saturn(k + 1, 1) = PCtemp(1, 373 + k);
        Cz_Saturn(k + 1, 1) = PCtemp(1, 380 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32.0, Cx_Saturn, Cy_Saturn, Cz_Saturn) * 1e3;

    Matrix Cx_Uranus(6, 1), Cy_Uranus(6, 1), Cz_Uranus(6, 1);
    for (int k = 0; k < 6; ++k) {
        Cx_Uranus(k + 1, 1) = PCtemp(1, 387 + k);
        Cy_Uranus(k + 1, 1) = PCtemp(1, 393 + k);
        Cz_Uranus(k + 1, 1) = PCtemp(1, 399 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32.0, Cx_Uranus, Cy_Uranus, Cz_Uranus) * 1e3;

    Matrix Cx_Neptune(6, 1), Cy_Neptune(6, 1), Cz_Neptune(6, 1);
    for (int k = 0; k < 6; ++k) {
        Cx_Neptune(k + 1, 1) = PCtemp(1, 405 + k);
        Cy_Neptune(k + 1, 1) = PCtemp(1, 411 + k);
        Cz_Neptune(k + 1, 1) = PCtemp(1, 417 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32.0, Cx_Neptune, Cy_Neptune, Cz_Neptune) * 1e3;

    Matrix Cx_Pluto(6, 1), Cy_Pluto(6, 1), Cz_Pluto(6, 1);
    for (int k = 0; k < 6; ++k) {
        Cx_Pluto(k + 1, 1) = PCtemp(1, 423 + k);
        Cy_Pluto(k + 1, 1) = PCtemp(1, 429 + k);
        Cz_Pluto(k + 1, 1) = PCtemp(1, 435 + k);
    }
    j = 0;
    Mjd0 = t1;
    Matrix r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32.0, Cx_Pluto, Cy_Pluto, Cz_Pluto) * 1e3;

    Matrix Cx_Nutations(40, 1), Cy_Nutations(40, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Nutations(k + 1, 1) = PCtemp(1, 819 + k);
        Cy_Nutations(k + 1, 1) = PCtemp(1, 829 + k);
    }
    for (int m = 1; m <= 3; ++m) {
        for (int k = 0; k < 10; ++k) {
            Cx_Nutations(10 * m + k + 1, 1) = PCtemp(1, 819 + 20 * m + k);
            Cy_Nutations(10 * m + k + 1, 1) = PCtemp(1, 829 + 20 * m + k);
        }
    }
    j = (dt <= 8.0) ? 0 : (dt <= 16.0) ? 1 : (dt <= 24.0) ? 2 : 3;
    Mjd0 = t1 + 8.0 * j;
    Matrix Cx_Nutations_sub(10, 1), Cy_Nutations_sub(10, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Nutations_sub(k + 1, 1) = Cx_Nutations(j * 10 + k + 1, 1);
        Cy_Nutations_sub(k + 1, 1) = Cy_Nutations(j * 10 + k + 1, 1);
    }
    Matrix Cz_Nutations_sub(10, 1); 
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8.0, 
                              Cx_Nutations_sub, Cy_Nutations_sub, Cz_Nutations_sub);

    Matrix Cx_Librations(40, 1), Cy_Librations(40, 1), Cz_Librations(40, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Librations(k + 1, 1) = PCtemp(1, 899 + k);
        Cy_Librations(k + 1, 1) = PCtemp(1, 909 + k);
        Cz_Librations(k + 1, 1) = PCtemp(1, 919 + k);
    }
    for (int m = 1; m <= 3; ++m) {
        for (int k = 0; k < 10; ++k) {
            Cx_Librations(10 * m + k + 1, 1) = PCtemp(1, 899 + 30 * m + k);
            Cy_Librations(10 * m + k + 1, 1) = PCtemp(1, 909 + 30 * m + k);
            Cz_Librations(10 * m + k + 1, 1) = PCtemp(1, 919 + 30 * m + k);
        }
    }
    j = (dt <= 8.0) ? 0 : (dt <= 16.0) ? 1 : (dt <= 24.0) ? 2 : 3;
    Mjd0 = t1 + 8.0 * j;
    Matrix Cx_Librations_sub(10, 1), Cy_Librations_sub(10, 1), Cz_Librations_sub(10, 1);
    for (int k = 0; k < 10; ++k) {
        Cx_Librations_sub(k + 1, 1) = Cx_Librations(j * 10 + k + 1, 1);
        Cy_Librations_sub(k + 1, 1) = Cy_Librations(j * 10 + k + 1, 1);
        Cz_Librations_sub(k + 1, 1) = Cz_Librations(j * 10 + k + 1, 1);
    }
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8.0, 
                               Cx_Librations_sub, Cy_Librations_sub, Cz_Librations_sub);

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
	
	std::cout << "r_Mercury: " << r_Mercury.n_row << "x" << r_Mercury.n_column << std::endl;
    std::cout << "r_Venus: " << r_Venus.n_row << "x" << r_Venus.n_column << std::endl;
    std::cout << "r_Earth: " << r_Earth.n_row << "x" << r_Earth.n_column << std::endl;
    std::cout << "r_Mars: " << r_Mars.n_row << "x" << r_Mars.n_column << std::endl;
    std::cout << "r_Jupiter: " << r_Jupiter.n_row << "x" << r_Jupiter.n_column << std::endl;
    std::cout << "r_Saturn: " << r_Saturn.n_row << "x" << r_Saturn.n_column << std::endl;
    std::cout << "r_Uranus: " << r_Uranus.n_row << "x" << r_Uranus.n_column << std::endl;
    std::cout << "r_Neptune: " << r_Neptune.n_row << "x" << r_Neptune.n_column << std::endl;
    std::cout << "r_Pluto: " << r_Pluto.n_row << "x" << r_Pluto.n_column << std::endl;
    std::cout << "r_Moon: " << r_Moon.n_row << "x" << r_Moon.n_column << std::endl;
    std::cout << "r_Sun: " << r_Sun.n_row << "x" << r_Sun.n_column << std::endl;

    return make_tuple(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);
}