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

int main()
{
    std::cout << "Hello World!\n";
}
