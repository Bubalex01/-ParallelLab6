#include <tbb/tbb.h>
#include <chrono>
#include <iostream>
using namespace std;

int Rrand(int min, int max) {
    return 1 + rand() % (max - min + 1) + min;
}

double* GaussTBB(double** A, double* Y, int n) {
    double* X = new double[n];

    for (int k = 0; k < n; k++)
    {
        tbb::parallel_for(tbb::blocked_range<int>(k + 1, n), [&](const tbb::blocked_range<int>& range) {
            for (int j = range.begin(); j < range.end(); j++) {
                double m = A[j][k] / A[k][k];
                for (int i = k; i < n; i++) {
                    A[j][i] = A[j][i] - m * A[k][i];
                }
                Y[j] = Y[j] - m * Y[k];
            }
            });
    }
    for (int k = n - 1; k >= 0; k--) {
        X[k] = Y[k];
        for (int i = k + 1; i < n; i++) {
            X[k] = X[k] - A[k][i] * X[i];

        }
        X[k] = X[k] / A[k][k];
    }
    return X;
}

double* FindY(double** A, double* X, int n) {
    double* Y = new double[n];
    for (int i = 0; i < n; i++) {
        Y[i] = 0;
        for (int j = 0; j < n; j++) {
            Y[i] += A[i][j] * X[j];
        }
    }
    return Y;
}

int ArrEquial(double* X1, double* X2, int n, double eps) {
    int c = 0;
    for (int i = 0; i < n; i++) {
        if (std::abs(X1[i] - X2[i]) > eps) {
            c++;
        }
    }
    return c;
}

double* Gauss(double** A, double* Y, int n) {
    double* X = new double[n];
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++) {
            double m = A[j][k] / A[k][k];
            for (int i = k; i < n; i++) {
                A[j][i] = A[j][i] - m * A[k][i];
            }
            Y[j] = Y[j] - m * Y[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        X[k] = Y[k];
        for (int i = k + 1; i < n; i++) {
            X[k] = X[k] - A[k][i] * X[i];

        }
        X[k] = X[k] / A[k][k];
    }

    return X;
}

int main()
{
    setlocale(LC_ALL, "Russian");

int nX=3200;
int f = 0;

double** A = new double* [nX];
for (int i = 0; i < nX; i++) {
    A[i] = new double[nX];
}
for (int i = 0; i < nX; i++) {
    for (int j = 0; j < nX; j++) {
        A[i][j] = 0;
    }
}

double* X = new double[nX];
double* X2 = new double[nX];
double* Y = new double[nX];

for (int i = 0; i < nX; i++) { //заполнение матрицы
    for (int j = 0; j < nX; j++) {
        A[i][j] = Rrand(-100, 100);
    }
}

for (int i = 0; i < nX; i++) { //заполнение правой части
    X2[i] = Rrand(-100, 100);
}

Y = FindY(A, X2, nX);

std::cout << "OneTBB решение" << std::endl;
auto begin = std::chrono::steady_clock::now();
X = GaussTBB(A, Y, nX);
auto end = std::chrono::steady_clock::now();
auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
std::cout << "Время: " << elapsed_ms.count() << " ms\n";
f = ArrEquial(X, X2, nX, 0.001);
std::cout << "Количество несовпадений:" << f << std::endl;

std::cout << "Прямое решение" << std::endl;
begin = std::chrono::steady_clock::now();
X = Gauss(A, Y, nX);
end = std::chrono::steady_clock::now();
elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
std::cout << "Время: " << elapsed_ms.count() << " ms\n";
f = ArrEquial(X, X2, nX, 0.001);
std::cout << "Количество несовпадений:" << f << std::endl;

for (int i = 0; i < nX; i++) {
    delete[] A[i];
}
delete[] A;
delete[] X;
delete[] Y;
}
