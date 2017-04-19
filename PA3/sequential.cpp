
#include <iostream>
#include <algorithm>
#include "stdio.h"
#include "stdlib.h"
#include <fstream>

#define N 4
using namespace std;

int main() {
    double **b = new double *[4];
    double **c = new double *[4];
    double **a = new double *[4];

    for (int i = 0; i < 4; i++) {
        c[i] = new double[4];
        b[i] = new double[4];
        a[i] = new double[4];
    }
    ifstream matrixfile("matrix");
    ifstream matrixfile2("matrix2");
    if (!(matrixfile.is_open())) {
        cout << "Error: file not found" << endl;
        return 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrixfile >> b[i][j];
            matrixfile2 >> a[i][j];
        }
    }
    matrixfile.close();
    matrixfile2.close();
    int i, j, k;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            c[i][j] = 0;
            for (k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }


    for (int i = 0; i < 4; i++) {
        cout << '\n';
        for (int j = 0; j < 4; j++) {
            printf("%lf\t", c[i][j]);
        }
    }

    for(int i = 0; i < 4; i++)
    {
        delete [] a[i];
        delete [] b[i];
        delete [] c[i];
    }
    return 0;
}