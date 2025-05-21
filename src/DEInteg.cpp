#include "..\include\DEInteg.hpp"

/**
 * Estados posibles del integrador de ecuaciones diferenciales.
 */
enum class DE_STATE {
    DE_INIT = 1,      
    DE_DONE = 2,     
    DE_BADACC = 3,   
    DE_NUMSTEPS = 4, 
    DE_STIFF = 5,    
    DE_INVPARAM = 6   
};

/**
 * Integra ecuaciones diferenciales ordinarias (ODE) usando un método de diferencias hacia atrás de paso variable (método de Gear).
 * @param f Función que define la ODE (dy/dt = f(t,y)).
 * @param t Tiempo inicial.
 * @param tout Tiempo final deseado.
 * @param relerr Tolerancia relativa de error.
 * @param abserr Tolerancia absoluta de error.
 * @param n_eqn Número de ecuaciones.
 * @param y Vector de estado inicial (n_eqn x 1).
 * @return Vector de estado en el tiempo tout.
 * @throw Error si los parámetros son inválidos o se accede a índices fuera de rango.
 */
Matrix DEInteg(Matrix f(double t, Matrix z), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y) {
	if(y.n_row<y.n_column){
		y=transpose(y);
	}
    const double twou = 2.0 * std::numeric_limits<double>::epsilon();
    const double fouru = 4.0 * std::numeric_limits<double>::epsilon();

    DE_STATE State_ = DE_STATE::DE_INIT;
    bool PermitTOUT = true;
    double told = 0.0;

    std::array<double, 14> two = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
    std::array<double, 14> gstr = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468};

    Matrix yy(n_eqn, 1), wt(n_eqn, 1), p(n_eqn, 1), yp(n_eqn, 1), phi(n_eqn, 17);
    Matrix g(14, 1), sig(14, 1), rho(14, 1), w(13, 1), alpha(13, 1);
    Matrix beta(13, 1), v(13, 1), psi_(13, 1);

    if (t == tout) {
        return y;
    }

    double epsilon = std::max(relerr, abserr);
    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 ||
        static_cast<int>(State_) > static_cast<int>(DE_STATE::DE_INVPARAM) ||
        (State_ != DE_STATE::DE_INIT && t != told)) {
        State_ = DE_STATE::DE_INVPARAM;
        return y;
    }

    double del = tout - t, absdel = std::abs(del);
    double tend = t + 100.0 * del;
    if (!PermitTOUT) tend = tout;

    int nostep = 0, kle4 = 0;
    bool stiff = false, start = false, phase1 = false, nornd = false, crash = false;
    double releps = relerr / epsilon, abseps = abserr / epsilon;
    double x = 0.0, h = 0.0, hold = 0.0, hnew = 0.0, delsgn = 0.0;
    double absh = 0.0;
    int k = 0, kold = 0, ns = 0, ifail = 0;
    int kp1 = 0, kp2 = 0, km1 = 0, km2 = 0, knew = 0;
    double erk = 0.0, erkm1 = 0.0, erkm2 = 0.0, erkp1 = 0.0, err = 0.0;
    bool OldPermit = false;

    if (State_ == DE_STATE::DE_INIT || !OldPermit || delsgn * del <= 0.0) {
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(std::max(fouru * std::abs(x), std::abs(tout - x)), tout - x);
    }

    while (true) {
        if (std::abs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1), ypout(n_eqn, 1);
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            double hi = tout - x;
            int ki = kold + 1;

            for (int i = 1; i <= ki; ++i) w(i + 1, 1) = 1.0 / i;
            double term = 0.0;
            for (int j = 2; j <= ki; ++j) {
                double psijm1 = psi_(j, 1);
                double gamma = (hi + term) / psijm1, eta = hi / psijm1;
                for (int i = 1; i <= ki + 1 - j; ++i)
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }

            for (int j = 1; j <= ki; ++j) {
                int i = ki + 1 - j;
                for (int l = 1; l <= n_eqn; ++l) {
                    yout(l, 1) += g(i + 1, 1) * phi(l, i + 1);
                    ypout(l, 1) += rho(i + 1, 1) * phi(l, i + 1);
                }
            }
            yout = y + yout * hi;
            y = yout;
            State_ = DE_STATE::DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        if (!PermitTOUT && std::abs(tout - x) < fouru * std::abs(x)) {
            h = tout - x;
            yp = f(x, yy);
            y = yy + yp * h;
            State_ = DE_STATE::DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        h = sign_(std::min(std::abs(h), std::abs(tend - x)), h);
        for (int l = 1; l <= n_eqn; ++l)
            wt(l, 1) = releps * std::abs(yy(l, 1)) + abseps;

        if (std::abs(h) < fouru * std::abs(x)) {
            h = sign_(fouru * std::abs(x), h);
            crash = true;
            return y;
        }

        double p5eps = 0.5 * epsilon;
        crash = false;
        g(2, 1) = 1.0;
        g(3, 1) = 0.5;
        sig(2, 1) = 1.0;

        double round = 0.0;
        for (int l = 1; l <= n_eqn; ++l)
            round += (y(l, 1) * y(l, 1)) / (wt(l, 1) * wt(l, 1));
        round = twou * std::sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return y;
        }

        if (start) {
            yp = f(x, y);
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum += (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));
            }
            sum = std::sqrt(sum);
            absh = std::abs(h);
            if (epsilon < 16.0 * sum * h * h) absh = 0.25 * std::sqrt(epsilon / sum);
            h = sign_(std::max(absh, fouru * std::abs(x)), h);
            hold = 0.0;
            hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; ++l) phi(l, 16) = 0.0;
            }
        }

        while (true) {
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            if (h != hold) ns = 0;
            if (ns <= kold) ns++;
            int nsp1 = ns + 1;

            if (k >= ns) {
                beta(ns + 1, 1) = 1.0;
                alpha(ns + 1, 1) = 1.0 / ns;
                double temp1 = h * ns;
                sig(nsp1 + 1, 1) = 1.0;
                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; ++i) {
                        int im1 = i - 1;
                        double temp2 = psi_(im1 + 1, 1);
                        psi_(im1 + 1, 1) = temp1;
                        beta(i + 1, 1) = beta(im1 + 1, 1) * psi_(im1 + 1, 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i + 1, 1) = h / temp1;
                        sig(i + 2, 1) = i * alpha(i + 1, 1) * sig(i + 1, 1);
                    }
                }
                psi_(k + 1, 1) = temp1;

                if (ns > 1) {
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v(k + 1, 1) = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 1; j <= nsm2; ++j) {
                            int i = k - j;
                            v(i + 1, 1) -= alpha(j + 2, 1) * v(i + 2, 1);
                        }
                    }
                    int limit1 = kp1 - ns;
                    double temp5 = alpha(ns + 1, 1);
                    for (int iq = 1; iq <= limit1; ++iq) {
                        v(iq + 1, 1) -= temp5 * v(iq + 2, 1);
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }
                    g(nsp1 + 1, 1) = w(2, 1);
                } else {
                    for (int iq = 1; iq <= k; ++iq) {
                        double temp3 = iq * (iq + 1);
                        v(iq + 1, 1) = 1.0 / temp3;
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }
                }

                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; ++i) {
                        int limit2 = kp2 - i;
                        double temp6 = alpha(i, 1);
                        for (int iq = 1; iq <= limit2; ++iq)
                            w(iq + 1, 1) -= temp6 * w(iq + 2, 1);
                        g(i + 1, 1) = w(2, 1);
                    }
                }
            }

            if (k >= nsp1) {
                for (int i = nsp1; i <= k; ++i) {
                    double temp1 = beta(i + 1, 1);
                    for (int l = 1; l <= n_eqn; ++l)
                        phi(l, i + 1) = temp1 * phi(l, i + 1);
                }
            }

            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, kp2 + 1) = phi(l, kp1 + 1);
                phi(l, kp1 + 1) = 0.0;
                p(l, 1) = 0.0;
            }
            for (int j = 1; j <= k; ++j) {
                int i = kp1 - j, ip1 = i + 1;
                double temp2 = g(i + 1, 1);
                for (int l = 1; l <= n_eqn; ++l) {
                    p(l, 1) += temp2 * phi(l, i + 1);
                    phi(l, i + 1) += phi(l, ip1 + 1);
                }
            }
            if (nornd) {
                p = y + p * h;
            } else {
                for (int l = 1; l <= n_eqn; ++l) {
                    double tau = h * p(l, 1) - phi(l, 16);
                    p(l, 1) = y(l, 1) + tau;
                    phi(l, 17) = (p(l, 1) - y(l, 1)) - tau;
                }
            }
            double xold = x;
            x += h;
            absh = std::abs(h);
            yp = f(x, p);

            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                double temp3 = 1.0 / wt(l, 1), temp4 = yp(l, 1) - phi(l, 2);
                if (km2 > 0)
                    erkm2 += std::pow((phi(l, km1 + 1) + temp4) * temp3, 2);
                if (km2 >= 0)
                    erkm1 += std::pow((phi(l, k + 1) + temp4) * temp3, 2);
                erk += std::pow(temp4 * temp3, 2);
            }

            if (km2 > 0) {
                if (static_cast<size_t>(km2) >= gstr.size())
                    throw std::out_of_range("Índice km2 fuera de rango en gstr");
                erkm2 = absh * sig(km1 + 1, 1) * gstr[km2] * std::sqrt(erkm2);
            }
            if (km2 >= 0) {
                if (static_cast<size_t>(km1) >= gstr.size())
                    throw std::out_of_range("Índice km1 fuera de rango en gstr");
                erkm1 = absh * sig(k + 1, 1) * gstr[km1] * std::sqrt(erkm1);
            }

            double temp5 = absh * std::sqrt(erk);
            err = temp5 * (g(k + 1, 1) - g(kp1 + 1, 1)); 
            if (static_cast<size_t>(k) >= gstr.size())
                throw std::out_of_range("Índice k fuera de rango en gstr");
            erk = temp5 * sig(kp1 + 1, 1) * gstr[k];
            knew = k;

            if (km2 > 0 && std::max(erkm1, erkm2) <= erk) knew = km1;
            if (km2 == 0 && erkm1 <= 0.5 * erk) knew = km1;

            bool success = err <= epsilon;

            if (!success) {
                phase1 = false;
                x = xold;
                for (int i = 1; i <= k; ++i) {
                    double temp1 = 1.0 / beta(i + 1, 1);
                    int ip1 = i + 1;
                    for (int l = 1; l <= n_eqn; ++l)
                        phi(l, i + 1) = temp1 * (phi(l, i + 1) - phi(l, ip1 + 1));
                }
                if (k >= 2) {
                    for (int i = 2; i <= k; ++i)
                        psi_(i, 1) = psi_(i + 1, 1) - h;
                }
                ifail++;
                double temp2 = 0.5;
                if (ifail > 3 && p5eps < 0.25 * erk) temp2 = std::sqrt(p5eps / erk);
                if (ifail >= 3) knew = 1;
                h *= temp2;
                k = knew;
                if (std::abs(h) < fouru * std::abs(x)) {
                    crash = true;
                    h = sign_(fouru * std::abs(x), h);
                    epsilon *= 2.0;
                    return y;
                }
            } else {
                break;
            }
        }

        kold = k;
        hold = h;

        double temp1 = h * g(kp1 + 1, 1);
        if (nornd) {
            for (int l = 1; l <= n_eqn; ++l)
                y(l, 1) = p(l, 1) + temp1 * (yp(l, 1) - phi(l, 2));
        } else {
            for (int l = 1; l <= n_eqn; ++l) {
                double rho = temp1 * (yp(l, 1) - phi(l, 2)) - phi(l, 17);
                y(l, 1) = p(l, 1) + rho;
                phi(l, 16) = (y(l, 1) - p(l, 1)) - rho;
            }
        }
        yp = f(x, y);

        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, kp1 + 1) = yp(l, 1) - phi(l, 2);
            phi(l, kp2 + 1) = phi(l, kp1 + 1) - phi(l, kp2 + 1);
        }
        for (int i = 1; i <= k; ++i) {
            for (int l = 1; l <= n_eqn; ++l)
                phi(l, i + 1) += phi(l, kp1 + 1);
        }

        erkp1 = 0.0;
        if (knew == km1 || k == 12) phase1 = false;

        if (!phase1 && kp1 <= ns) {
            for (int l = 1; l <= n_eqn; ++l)
                erkp1 += std::pow(phi(l, kp2 + 1) / wt(l, 1), 2);
            if (static_cast<size_t>(kp1) >= gstr.size())
                throw std::out_of_range("Índice kp1 fuera de rango en gstr");
            erkp1 = absh * gstr[kp1] * std::sqrt(erkp1);
            if (k > 1) {
                if (erkm1 <= std::min(erk, erkp1)) {
                    k = km1;
                    erk = erkm1;
                } else if (erkp1 < erk && k != 12) {
                    k = kp1;
                    erk = erkp1;
                }
            } else if (erkp1 < 0.5 * erk) {
                k = kp1;
                erk = erkp1;
            }
        }

        if (phase1 || p5eps >= erk * two[k]) {
            hnew = 2.0 * h;
        } else if (p5eps < erk) {
            double temp2 = k + 1;
            double r = std::pow(p5eps / erk, 1.0 / temp2);
            hnew = absh * std::max(0.5, std::min(0.9, r));
            hnew = sign_(std::max(hnew, fouru * std::abs(x)), h);
        } else {
            hnew = h;
        }
        h = hnew;

        if (crash) {
            State_ = DE_STATE::DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            y = yy;
            t = x;
            told = t;
            OldPermit = true;
            return y;
        }

        nostep++;
        kle4++;
        if (kold > 4) kle4 = 0;
        if (kle4 >= 50) stiff = true;
    }

    if (State_ == DE_STATE::DE_INVPARAM) {
        std::cerr << "Parámetros inválidos en DEInteg\n";
    } else if (State_ == DE_STATE::DE_BADACC) {
        std::cerr << "No se logró la precisión requerida en DEInteg\n";
    } else if (State_ == DE_STATE::DE_STIFF) {
        std::cerr << "Se sospecha un problema rígido en DEInteg\n";
    }

    return y;
}