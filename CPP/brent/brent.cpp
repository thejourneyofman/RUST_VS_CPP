// brent.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <array>
#include <random>
#include <chrono>

class Brent {
public:
    template <class F>
    double solveImpl(const F& f, double v, double x1, double x2, double tol) const {

        double a = x1;
        double b = x2;
        double c = x2;

        double fa = f(a, v);
        double fb = f(b, v);
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
            throw "Root must be bracketed in zbrent";
        }

        double d;
        double e;
        double q;
        double r;
        double p;

        double fc = fb;
        double index = 0;

        while (evaluationNumber_ <= maxEvaluations_) {

            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                // rename a, b, c and adjust bounding interval d
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
            if (std::fabs(fc) < std::fabs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            // convergence check
            double tol1 = 2.0 * std::numeric_limits<double>::epsilon() * std::fabs(b) + 0.5 * tol;
            double xm = 0.5 * (c - b);
            if (std::fabs(xm) <= tol1 || fb == 0.0) {
                return b;
            }
            if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
                // attempt inverse quadratic interpolation
                double s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                // check whether in bounds
                if (p > 0.0) {
                    q = -q;
                }
                p = std::fabs(p);
                double min1 = 3.0 * xm * q * std::fabs((tol1 - q));
                double min2 = std::fabs((e * q));
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    // accept interpolation
                    e = d;
                    d = p / q;
                }
                else {
                    // interpolation failed, use bisection
                    d = xm;
                    e = d;
                }
            }
            else {
                // bounds decreasing too slowly, use bisection
                d = xm;
                e = d;
            }

            // move last best guess to a
            a = b;
            fa = fb;
            // evaluate new trial root
            if (std::fabs(d) > tol1) {
                b += d;
            }
            else {
                b += sign(tol1, xm);
            }

            fb = f(b, v);

            evaluationNumber_++;
        }
        
        std::cout << "maximum number of function evaluations ("
            << maxEvaluations_ << ") exceeded ";
    }
private:
    double sign(double a, double b) const {
        return b >= 0.0 ? std::fabs(a) : -std::fabs(a);
    }


protected:
    mutable double root_, xMin_, xMax_, fxMin_, fxMax_;
    size_t maxEvaluations_ = 100;
    mutable size_t evaluationNumber_ = 0;
};

double func1(double x, double v) {
    return sin(x) - v;
}

void brent_arcsin() {

    std::array<double, 5> samples = { 0.0, 0.01, 0.02, 0.03, 0.04 };
    double min = 0.0;
    double max = M_PI * 0.5;
    double tol = 1e-10;
    double max_iter = 100;

    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);

    const int nSimulations = 1000000;
    for (int i = 0; i < nSimulations; i++)
    {
        double v = unif(rng);
        double y = Brent().solveImpl(func1, v, min, max, tol);
        //double expected = asin(v);
        //std::cout << "result=" << v << " expected=" << expected << std::endl;
    }

    //for (auto v : samples) {
    //    double y = Brent().solveImpl(func1, v, min, max, tol);
    //    double expected = asin(v);
    //    printf("result = %.18f expected= %.18f \n", v, expected);
    //}
}

int main()
{
    auto beginTime = std::chrono::high_resolution_clock::now();
    brent_arcsin();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);
    std::cout << "Time elapsed in brent_arcsin Solver with C++ is: " << elapsedTime.count() << " milliseconds" << std::endl;   
   
}




