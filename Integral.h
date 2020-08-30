#pragma once
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

namespace Integral {
    //(i: 0- cos, 1- sin), N- accuracy, a- lower limit, b- upper limit, p- trigonometric function's parameter, f- function
    double Calc(int i, int N, double a, double b, int p, double(*f)(double)) {
        if ((N % 2) != 0) {
            throw (string)"N%2 must be equal to 0";
        }
        if (N <= 0) {
            throw (string)"N must be more than 0";
        }
        if (i != 0 && i != 1) {
            throw (string)"Incorrect value of i";
        }
        double J;
        double h = (b - a) / N;
        double lambda = 4 * h * pow(p*h, -2)*(pow(p*h, -1)*sin(p*h) - cos(p*h));
        double mu = h * pow(p*h, -2)*(1 + pow(cos(p*h), 2) - pow(p*h, -1)*sin(2 * p*h));
        double nu = (1 + pow(2 * p*h, -1)*sin(2 * p*h) - pow(p*h, -2) * 2 * pow(sin(p*h), 2)) / p;
        //
        if (!i) {
            J = mu * cos(p*a) - nu * sin(p*a);
            J *= f(a);
            for (int m = 1; m <= N / 2 - 1; m++) {
                J += 2 * mu*(f(a + 2 * m*h)*(cos(p*(a + 2 * m*h))));
            }
            for (int m = 0; m <= N / 2 - 1; m++) {
                J += lambda * f(a + (2 * m + 1)*h)*cos(p*(a + (2 * m + 1)*h));
            }
            J += (mu*cos(p*b) + nu * sin(p*b))*f(b);
        }
        else {
            J = mu * sin(p*a) + nu * cos(p*a);
            J *= f(a);
            for (int m = 1; m <= N / 2 - 1; m++) {
                J += 2 * mu*(f(a + 2 * m*h)*(sin(p*(a + 2 * m*h))));
            }
            for (int m = 0; m <= N / 2 - 1; m++) {
                J += lambda * f(a + (2 * m + 1)*h)*sin(p*(a + (2 * m + 1)*h));
            }
            J += (mu*sin(p*b) - nu * cos(p*b))*f(b);
        }
        return J;
    }
}