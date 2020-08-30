#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include "Integral.h"
#include "fGNV.h"


using namespace std;

# define M_PI 3.141592653589793238462643383279502884L

//Глобальные переменные для функции
static double g_ld1, g_ld2, g_bi, g_phi, g_nd, g_nmax, g_r, g_w2, g_lastn;
//Интегранд
double integrand(double teta) {
    short error;
    double lambdaTeta = sqrt(g_ld1 * pow(cos(teta), 2) + g_ld2 * pow(sin(teta), 2));
    double muS = 3 * (1 + g_bi);
    double aTeta = sqrt(muS) / lambdaTeta;
    return (1 / pow(lambdaTeta, 2))*fGNV(g_lastn, g_lastn, aTeta*g_r, error);
}
//Интегранд при n = 0
double integrand2(double teta) {
    short error;
    double lambdaTeta = sqrt(g_ld1 * pow(cos(teta), 2) + g_ld2 * pow(sin(teta), 2));
    double muS = 3 * (1 + g_bi);
    double aTeta = sqrt(muS) / lambdaTeta;
    return (1 / pow(lambdaTeta, 2)*fGNV(0, 0, aTeta*g_r, error)) / cos(2 * teta);
}
//Вычисление температуры в точке с заданными параметрами в полярной системе координат
double Calc(double ld1, double ld2, double bi, double phi, int nd, int nmax, double r, double w2 = 1) {
    g_ld1 = ld1;
    g_ld2 = ld2;
    g_bi = bi;
    g_phi = phi;
    g_nd = nd;
    g_nmax = nmax;
    g_r = r;
    g_w2 = w2;
    double result = (3 * g_w2) / pow(M_PI, 2), sum = Integral::Calc(0, g_nd, 0, M_PI / 2, 2, integrand2);
    for (int n = 1; n < g_nmax; n++) {
        g_lastn = n;
        sum += 2 * cos(2 * n*g_phi)*Integral::Calc(0, g_nd, 0, M_PI / 2, 2 * n, integrand);
    }
    return result * sum;
}
struct Point {
    double x, y;
};
//Вычисление температуры в заданном квадрате с заданными параметрами и вывод результата в текстовый файл
void CalculateTemperature(double a, int n, int count, string path, double ld1, double ld2, double bi) {
    double r, phi;
    ofstream fout(path + "Temp.txt");
    ofstream foutA(path + "A.txt");
    if (fout.is_open() && foutA.is_open()) {
        double** result = new double*[n + 1];
        for (int i = 0; i < n + 1; i++) {
            result[i] = new double[n + 1];
            for (int j = 0; j < n + 1; j++) {
                result[i][j] = 0;
            }
        }
        Point A;
        double w2;
        int nd = 24, nmax = 11;
        for (int i = 0; i < count; i++) {
            cout << "Точка " << i + 1 << "\n\tЦентр источника тепла\n\t\tx = ";
            cin >> A.x;
            cout << "\t\ty = ";
            cin >> A.y;
            cout << "\tИнтенсивность источника тепла\n";
            cout << "\t\tW2 = ";
            cin >> w2;
            double x = -a / 2 - A.x, y = -a / 2 - A.y;
            for (int j = 0; j < n + 1; j++) {
                y = -a / 2 - A.y;
                for (int k = 0; k < n + 1; k++) {
                    r = sqrt(pow(x, 2) + pow(y, 2));
                    if (x == 0) {
                        if (y < 0) {
                            phi = (3 * M_PI) / 2;
                        } if (y > 0) {
                            phi = M_PI / 2;
                        }
                        else {
                            r = 0.00001;
                            phi = atan(0) + M_PI;
                        }
                    }
                    else if (x < 0) {
                        phi = atan(y / x) + M_PI;
                    }
                    else phi = atan(y / x);
                    result[j][k] += Calc(ld1, ld2, bi, phi, nd, nmax, r, w2);
                    y = y + a / n;
                }
                x = x + a / n;
            }
        }
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < n + 1; j++) {
                fout << setw(9) << setiosflags(ios::fixed) << setprecision(5) << result[i][j];
            }
            fout << "\n";
        }
        fout.close();
        foutA << a << "\n" << n;
        foutA.close();
    }
    else {
        cout << "Ошибка доступа к файлу" << endl;
    }
}
int main() {
    setlocale(LC_ALL, "Rus");
    int count, n;
    double a, ld1, ld2, bi;
    cout << "Введите сторону квадрата a: ";
    cin >> a;
    cout << "Введите количество отрезков разбиения n: ";
    cin >> n;
    cout << "Введите количество источников тепла: ";
    cin >> count;
    cout << "Коэфициент теплопроводности Ld1 = ";
    cin >> ld1;
    cout << "Коэфициент теплопроводности Ld2 = ";
    cin >> ld2;
    cout << "Коэфициент теплообмена с внешней средой Bi = ";
    cin >> bi;
    CalculateTemperature(a, n, count, "", ld1, ld2, bi);
}