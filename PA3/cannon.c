#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <sys/time.h>
#include <time.h>
#include <string.h>

static void matrix_malloc(double **, int, int, int);

static void MatrixMultiplication(int, int, int, double *, double *, double *, MPI_Comm comm);

static void MM(double *, double *, double *, int, int, int);

int main(int argc, char **argv) {

    double elapsedTime;
    struct timeval start, end;
    FILE *matrix, *matrix2;
    matrix = fopen("matrix", "r");
    matrix2 = fopen("matrix2", "r");
    if (matrix == NULL || matrix2 == NULL) {
        printf("Error Reading File\n");
        exit(0);
    }
    int rank, numprocs;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int sqr = sqrt(numprocs - 1);

    MPI_Comm new_commuicator;

    // int color = rank / sqrt(numprocs);

    MPI_Comm_split(MPI_COMM_WORLD, rank != 0 ? 1 : MPI_UNDEFINED, rank, &new_commuicator);

    // if(&new_commuicator == NULL) printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

    char *p, *q;
    double *a, *b, *c;
    int row, column;
    if (rank == 0) {
//        row = atoi(argv[1]);
//        column = atoi(argv[2]);
//
        if (argc < 3 || argc > 4) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (argc == 4) {
            row = (int) strtol(argv[2], &p, 10);
            column = (int) strtol(argv[3], &p, 10);
            char *random_flag = argv[1];

            int size = row * column;
            if (row != column || strcmp(random_flag, "-t")) {
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }

            a = (double *) malloc(sizeof(double) * size);
            b = (double *) malloc(sizeof(double) * size);


            for (int i = 0; i < size; i++) {
                fscanf(matrix, "%lf,", &a[i]);
                fscanf(matrix2, "%lf,", &b[i]);
            }
        } else if (argc == 3) {
            row = (int) strtol(argv[1], &p, 10);
            column = (int) strtol(argv[2], &p, 10);

            if (row != column) {
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }

            matrix_malloc(&a, row, column, 1);
            matrix_malloc(&b, row, column, 1);
        }

        matrix_malloc(&c, row, column, 0);
    }

//    (*a) = (double *) malloc(row * column * sizeof(double));
//    (*b) = (double *) malloc(row * column * sizeof(double));
//    (*c) = (double *) malloc(row * column * sizeof(double));
//
//
//    for (int i = 0; i < row * column; i++) {
//        (*a)[i] = (double) rand();
//        (*b)[i] = (double) rand();
//        (*c)[i] = 0.0;
//    }

    gettimeofday(&start, NULL);
    if (rank == 0) {
        //       int sqr = sqrt(numprocs - 1);
        int A_i_offset = 0, B_i_offset = 0, A_j_offset = 0, B_j_offset = 0;
        for (int i = 1; i < numprocs; i++) {
            int Pi = (i - 1) / sqr;
            int Pj = (i - 1) % sqr;

            //  printf("Pi is %d and Pj is %d\n", Pi, Pj);
//            int my_m = Pi
            int my_m = Pi < row % sqr ? row / sqr + 1 : row / sqr;

            int my_n = Pj < column % sqr ? column / sqr + 1 : column / sqr;

            int my_l_a = Pj < column % sqr ? column / sqr + 1 : column / sqr;

            int my_l_b = Pi < column % sqr ? column / sqr + 1 : column / sqr;


            int my_m_max = 0 < row % sqr ? row / sqr + 1 : row / sqr;
            int my_l_max = 0 < column % sqr ? column / sqr + 1 : column / sqr;
            int my_n_max = 0 < column % sqr ? column / sqr + 1 : column / sqr;

            double *subA, *subB;

            matrix_malloc(&subA, my_m_max, my_l_max, 0);
            matrix_malloc(&subB, my_l_max, my_n_max, 0);

//            int sizeA = sizeof(subA);
//            int sizeB = sizeof(sizeB);
//            printf("subA has size %d and subB has size %d\n", sizeA, sizeB);

//            (*subA) = (double *) malloc(my_m_max * my_l_max * sizeof(double));
//            (*subB) = (double *) malloc(my_l_max * my_m_max * sizeof(double));
//
//            for (int i = 0; i < my_m_max * my_l_max; i++) {
//                (*subA)[i] = 0.0;
//                (*subB)[i] = 0.0;
//            }

            for (int i = 0; i < my_m; i++) {
                for (int j = 0; j < my_l_a; j++) {
                    subA[i * my_l_max + j] = a[(i + A_i_offset) * column + (j + A_j_offset)];
                }
//                printf("subA has i is %d and value is %lf and a at index %d\n", i * my_l_max + j,
//                       subA[i * my_l_max + j], (i + A_i_offset) * column + (j + A_j_offset));
            }
            for (int i = 0; i < my_l_b; i++) {
                for (int j = 0; j < my_n; j++) {
                    subB[i * my_n_max + j] = b[(i + B_i_offset) * column + (j + B_j_offset)];
                }
//                printf("subB has i is %d and value is %lf and a at index %d\n", i * my_n_max + j,
//                       subB[i * my_n_max + j], (i + B_i_offset) * column + (j + B_j_offset));
            }

            //A_j_offset = B_j_offset = 0;
            if (i % sqr == 0) {
                A_j_offset = B_j_offset = 0;
                A_i_offset += my_m;
                B_i_offset += my_l_b;
            } else {
                //  A_j_offset = B_j_offset = 0;
                A_j_offset += my_l_a;
                B_j_offset += my_n;

            }
//            printf("a_i_off is %d and a_j_off is %d and b_i_off is %d and b_j_off is %d\n", A_i_offset, A_j_offset,
//                   B_i_offset, B_j_offset);



            MPI_Send(&my_m_max, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&my_l_max, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&my_n_max, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            MPI_Send(subA, my_m_max * my_l_max, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(subB, my_l_max * my_n_max, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

            free(subA);
            free(subB);
        }

        /*******************receive*****************/

        int C_i_offset = 0, C_j_offset = 0;

        MPI_Status temp_status;
        for (int i = 1; i < numprocs; i++) {

            int Pi = (i - 1) / sqr;
            int Pj = (i - 1) % sqr;
            //  printf("Pi is %d and Pj is %d\n", Pi, Pj);

            //int my_m = Pi < column % sqr -1  ? column / sqr + 1 : column / sqr - 1;


            int my_m = Pi < row % sqr ? row / sqr + 1 : row / sqr;
            int my_n = Pj < column % sqr ? column / sqr + 1 : column / sqr;

            int my_m_max = 0 < row % sqr ? row / sqr + 1 : row / sqr;
            int my_n_max = 0 < column % sqr ? column / sqr + 1 : column / sqr;


            double *subC;

            matrix_malloc(&subC, my_m_max, my_n_max, 0);

            //    printf("!!!!!!!!!!!!!size of subC is %d\n", sizeof(subC));

//                (*subC) = (double *) malloc(my_m_max * my_n_max * sizeof(double));
//
//                for (int i = 0; i < my_m_max * my_n_max; i++) {
//                    (*subC)[i] = 0.0;
//                }

            MPI_Recv(subC, my_m_max * my_n_max, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &temp_status);

//            for(int i = 0; i < my_m_max * my_n_max; i++)
//                printf("receive subC has %lf\n", subC[i]);



            for (int i = 0; i < my_m; i++) {
                for (int j = 0; j < my_n; j++) {

                    c[(i + C_i_offset) * column + (j + C_j_offset)] += subC[i * my_n_max + j];
                }

//                printf("subC has i is %d and value is %lf and c at index %d\n", i * my_n_max + j,
//                       subC[i * my_n_max + j], (i + C_i_offset) * column + (j + C_j_offset));
            }

            if (i % sqr == 0) {
                C_i_offset += my_m;
                C_j_offset = 0;
            } else {
                //    C_i_offset = 0;
                C_j_offset += my_n;
            }

            //         printf("c_i_off is %d and c_j_off is %d\n", C_i_offset, C_j_offset);

            free(subC);

        }

    } else {
        int my_m, my_l, my_n;

        MPI_Recv(&my_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&my_l, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&my_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        //    printf("!!!!!!Receive row for m1 is %d, column for m1 is %d, column , column for m2 is %d !!!!!!!\n", my_m, my_n, my_l);

        double *a_temp, *b_temp, *c_temp;

        matrix_malloc(&a_temp, my_m, my_l, 0);
        matrix_malloc(&b_temp, my_l, my_n, 0);
        matrix_malloc(&c_temp, my_m, my_n, 0);

//        (*a_temp) = (double *) malloc(my_m * my_l * sizeof(double));
//        (*b_temp) = (double *) malloc(my_l * my_n * sizeof(double));
//        (*c_temp) = (double *) malloc(my_m * my_n * sizeof(double));
//
//        for (int i = 0; i < my_m_max * my_l_max; i++) {
//            (*subA)[i] = 0.0;
//            (*subB)[i] = 0.0;
//        }

        MPI_Recv(a_temp, my_m * my_l, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(b_temp, my_l * my_n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//        for (int i = 0; i < my_m * my_l; i++)
//            printf("receive a_temp has %lf       b_temp has %lf\n", a_temp[i], b_temp[i]);

        MatrixMultiplication(my_m, my_l, my_n, a_temp, b_temp, c_temp, new_commuicator);


//        for (int i = 0; i < my_m * my_l; i++)
//            printf("now c_temp has %lf\n", c_temp[i]);

        MPI_Send(c_temp, my_m * my_n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        //    if(rank == 0) {
        free(a_temp);
        free(b_temp);
        free(c_temp);
        //   }

    }


    int size = sizeof(c);
    if (rank == 0) {
        int i, j, count = 0;
        double **arr = (double **) malloc(row * sizeof(double *));
        for (i = 0; i < row; i++)
            arr[i] = (double *) malloc(column * sizeof(double));

        for (i = 0; i < row; i++) {
            for (j = 0; j < column; j++) {
                arr[i][j] = c[count];
                count++;
            }
        }

        for (i = 0; i < row; i++) {
            printf("\n");
            for (j = 0; j < column; j++) {
                printf("%lf\t", arr[i][j]);
            }
        }
        gettimeofday(&end, NULL);
        elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
        elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
        printf("\nTime: %lf ms.\n", elapsedTime);
        free(a);
        free(b);
        free(c);

        for (i = 0; i < row; i++)
            free(arr[i]);
        free(arr);
    }

    MPI_Finalize();
    return 0;
}

static void matrix_malloc(double **matrix, int row, int column, int flag) {
    (*matrix) = (double *) malloc(row * column * sizeof(double));
    if (flag == 0) {
        for (int i = 0; i < row * column; i++) (*matrix)[i] = 0.0;
    } else {
        srand((unsigned int) time(NULL));
        for (int i = 0; i < row * column; i++) (*matrix)[i] = (double) rand();
    }
}

static void MatrixMultiplication(int my_m, int my_l, int my_n, double *a, double *b, double *c, MPI_Comm comm) {

    int i;
    int npes;
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank;
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &myrank);


    int dims[2] = {sqrt(npes), sqrt(npes)};
    int periods[2] = {1, 1};

    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);


    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, my_m * my_l, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, my_l * my_n, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);


    for (i = 0; i < dims[0]; i++) {
        MM(a, b, c, my_m, my_l, my_n);

        MPI_Sendrecv_replace(a, my_m * my_l, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, &status);
        MPI_Sendrecv_replace(b, my_l * my_n, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);
    }

    MPI_Comm_free(&comm_2d);
}

static void MM(double *A, double *B, double *C, int m, int l, int n) {
    double *B_temp;
    int i, j, k;
    B_temp = (double *) malloc(l * n * sizeof(double));
    for (i = 0; i < l; i++) {
        for (j = 0; j < n; j++) {
            B_temp[j * l + i] = B[i * n + j];
        }
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            double dot = 0.0;
            for (k = 0; k < l; k++) {
                dot += A[i * l + k] * B_temp[j * l + k];
            }
            C[i * n + j] += dot;
        }
    }
    free(B_temp);
}

//
// Created by cqx820 on 4/4/17.
//

