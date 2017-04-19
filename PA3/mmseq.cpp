#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/time.h>

using namespace std;

int main(int argc, char **argv) {

    timeval start, end;
    double elapsedTime;
	
    int row = atoi(argv[1]);
    int column = atoi(argv[2]);
    if(row != column) return 0;
    double **a = new double *[row];
    double **b = new double *[row];
    double **c = new double *[row];
    for (int i = 0; i < row; i++) {
        a[i] = new double[column];
        b[i] = new double[column];
        c[i] = new double[column];
    }

    ifstream matrixfile("matrix");
    ifstream matrixfile2("matrix2");

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            matrixfile >> a[i][j];
            matrixfile2 >> b[i][j];
        }
    }

    matrixfile.close();
    matrixfile2.close();
    gettimeofday(&start, NULL);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            c[i][j] = 0.;
            for (int k = 0; k < row; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

//    gettimeofday(&end, NULL);


    for (int i = 0; i < row; i++) {
        cout << '\n';
        for (int j = 0; j < column; j++) { printf("%lf\t", c[i][j]); }
    }
    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
    printf("\nTime: %lf ms.\n", elapsedTime);

    for (int i = 0; i < row; i++) {
        delete [] a[i];
        delete [] b[i];
        delete [] c[i];
    }

    return 0;
}
