#include <iostream>
#include "/usr/include/mpi/mpi.h"
#include <algorithm>
#include "stdio.h"
#include "stdlib.h"
#include <fstream>

using namespace std;


void printArray(double **A, double **B, int rank) {

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 4; j++)
            printf("Rank: %d The number in array A is: %lf and the number in array B is: %lf\n", rank, A[j][i],
                   B[i][j]);
}

int main(int argc, char **argv) {
    FILE *matrix, *matrix2;
    matrix = fopen("matrix", "r");
    matrix2 = fopen("matrix2", "r");
    if (matrix == NULL || matrix2 == NULL) {
        printf("Error Reading File\n");
        exit(0);
    }
    //  double **b = new double *[4];
    double **c = new double *[4];
    double **myB = new double *[2];
    double **myA = new double *[4];
    //  double **a = new double *[4];


    double *temp_A = (double *) malloc(sizeof(double) * 8);
    double *temp_B = (double *) malloc(sizeof(double) * 8);


    double *a_lin = (double *) malloc(sizeof(double) * 16);
    double *b_lin = (double *) malloc(sizeof(double) * 16);
    double *a_tem = (double *) malloc(sizeof(double) * 16);


    for (int i = 0; i < 4; i++) {
        c[i] = new double[4];
        myA[i] = new double[2];
        //    b[i] = new double[4];
        //    a[i] = new double[4];
    }

    for (int i = 0; i < 2; i++)
        myB[i] = new double[4];

    for (int i = 0; i < 16; i++) {
        fscanf(matrix, "%lf,", &a_tem[i]);
        fscanf(matrix2, "%lf,", &b_lin[i]);
    }

    int num = 0;
    for (int i = 0; i < 16; i += 4) {
        a_lin[num] = a_tem[i];
        a_lin[num + 1] = a_tem[i + 1];
        num += 2;
    }

    for (int i = 2; i < 16; i += 4) {
        a_lin[num] = a_tem[i];
        a_lin[num + 1] = a_tem[i + 1];
        num += 2;
    }

//    //   free(a_tem);
//    for (int j = 0; j < 16; j++) {
//        printf("we have %f\n", a_lin[j]);
//    }

//	for(int i = 0; i < 4; i++)
//		for(int j = 0; j < 4; j++)
//			printf("we have %f\n", a[i][j]);

//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 2; j++) {
//            //buffer_A[i][j] = a[i][j];
//        }
//    }
//
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 4; j++) {
//            //buffer_B[i][j] = b[i][j];
//        }
//    }


    int rank, numprocs;
    // int namelen;
    // char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//    MPI_Get_processor_name(processor_name,&namelen);

    //  fprintf(stderr,"Hello World! Process %d of %d on %s\n",
    //        rank, numprocs, processor_name);


    MPI_Scatter(a_lin, 8, MPI_DOUBLE, temp_A, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b_lin, 8, MPI_DOUBLE, temp_B, 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // for(int i = 0; i < 8; i++)
////        for(int j = 0; j < 2; j++ )
    //           printf("%d rank: %f\n", rank,  temp_B[i]);

    //   printArray(myA, myB, rank);
    MPI_Request request;
    MPI_Status status;


//
//    for(int i = 0; i < 4; i++)
//    {
//        for(int j = 2; j < 4; j++)
//        {
//            buffer_A[i][j-2] = a[i][j];
//        }
//    }
//
//    for(int i = 2; i < 4; i++)
//    {
//        for(int j = 0; j < 4; j++)
//        {
//            buffer_B[i-2][j] = b[i][j];
//        }
//    }

    // if (rank == 0) {



//        for (int i = 0; i < 4; i++) {
//            cout << '\n';
//            for (int j = 0; j < 2; j++) {
//                printf("%lf\t", myA[i][j]);
//            }
//        }

//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++) {
//            c[i][j] = 0;
//            for (int k = 0; k < 4; k++) {
//                //   printf("This is %lf * %lf\n", myB[i][k], myA[k][j]);
//                c[i][j] += myB[i][k] * myA[k][j];
//            }
//        }
//    }
    double *abuf = (double *) malloc(sizeof(double) * 8);

    MPI_Isend(temp_A, 8, MPI_DOUBLE, (rank + 1) % numprocs, 0, MPI_COMM_WORLD, &request);
    MPI_Irecv(abuf, 8, MPI_DOUBLE, (rank + 1) % numprocs, 0, MPI_COMM_WORLD, &request);
    //numprocs - 1
    double **buffer_A = new double *[4];

    for (int i = 0; i < 4; i++) {
        buffer_A[i] = new double[2];
    }

    int count = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            myB[i][j] = temp_B[count];
            count++;
        }
    }

    count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            myA[i][j] = temp_A[count];
            //     cout<<"send buffer i is " <<i<<" "<<myA[i][j]<<" in rank "<<rank<<endl;
            count++;
        }
    }


    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            //c[0][]=0;

            c[i + rank * 2][j + rank * 2] = 0;
            for (int k = 0; k < 4; k++) {
                c[i + rank * 2][j + rank * 2] += myB[i][k] * myA[k][j];

            }
        }
    }

    MPI_Wait(&request, &status);

    count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            buffer_A[i][j] = abuf[count];
            count++;
        }
    }

//        for(int i = 0; i < 4; i++)
//        {
//            delete [] buffer_A[i];
//        }

    // } else {
//
//    double **buffer_B = new double *[4];
//
//    for (int i = 0; i < 4; i++) {
//        buffer_B[i] = new double[2];
//    }
//
//    count = 0;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 4; j++) {
//            myB[i][j] = temp_B[count];
//            count++;
//        }
//    }
//
//    count = 0;
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 2; j++) {
//            myA[i][j] = temp_A[count];
//            count++;
//        }
//    }
//

//
//    free(temp_A);
//    free(temp_B);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            c[i + ((rank) % numprocs) * 2][j + ((rank + 1) % numprocs) * 2] = 0;

            for (int k = 0; k < 4; k++) {

                c[i + ((rank) % numprocs) * 2][j + ((rank + 1) % numprocs) * 2] += myB[i][k] * buffer_A[k][j];
            }
        }
    }
//
//    if(rank == 0) {
//
//        for (int i = 0; i < 4; i++) {
//            cout << '\n';
//            for (int j = 0; j < 2; j++) {
//                printf("%lf\t", buffer_A[i][j]);
//            }
//        }
//    }



    //    MPI_Irecv(&myA[0][0], 8, MPI_DOUBLE, (rank + numprocs - 1) % numprocs, 1, MPI_COMM_WORLD, request);
    //  MPI_Isend(&myA[0][0], 8, MPI_DOUBLE, (rank + numprocs - 1) % numprocs, 2, MPI_COMM_WORLD, request);//all 0 tag
    //mpi wait
//    for (int i = 2; i < 4; i++) {
//        for (int j = 0; j < 2; j++) {
//            c[i][j] = 0;
//            for (int k = 0; k < 4; k++) {
//                c[i][j] += myB[i - 2][k] * myA[k][j];
//            }
//        }
//    }


//        for(int i = 0; i < 4; i++)
//        {
//            delete [] buffer_B[i];
//        }
    //   }

    // MPI_Waitall(8, request, status);

    //for(i = 0; i < 4; i++)
    // {
    //  for(j = 0; j < 4; j++)
    //   {
    //c[i][j] = 0;
    //for(k = 0; k < 4; k++)
    //{
    //      c[i][j] += b[i][k] * a[k][j];
    //    }
    //  }
    //}


    //myc[i+rank*N/P][j], myc[i+((rank+1)%p)*N/P][j]

//	prinf("PAR:c[%d][%d] %f\n", i, j+rank*N/P, myc[i][j];
//	MPI_Request* request = malloc(sizeof(MPI_Request) * 2);
//	MPI_Status* status = malloc(sizeof(MPI_Status) * 2);
//	if(myid == 0)
//	{
//		MPI_Isend(&b[0][0], 16/2, MPI_DOUBLE, myid +1 % p, 0, MPI_COMM_WORLD, &request[0]);
//		MPI_Wait(request, status);
//
//	}
//	else if(myid == 1)
//	{
//		MPI_Irecv(void* buf, int count, MPI_DOUBLE, int source, 0, MPI_COMM_WORLD, request[1]);
//		MPI_Wait(request, status);
//	}

//	int send_count  = recv_count = (4 * 4) / numprocs;
//	MPI_Comm communicator = MPI_COMM_WORLD;
//	int root = 0;
//	MPI_Waitall;





    MPI_Finalize();

    //  for (int i = 0; i < 4; i++) {
    //cout << '\n';
    //      for (int j = 0; j < 4; j++) { printf("%lf\t", c[i][j]); }
    // }
    return 0;
}



//MPI_Scatter(&a, int , MPI_DOUBLE, void* recv_data, int recv_count, MPI_DOUBLE, RootProcess, MPI_Comm communicator);
//Scatter()
//sendcount  N * N / P
//MPI_FLOAT

//double myA[N/P][N]
//double myB[N][N/P]

//mpiScatter
//mpi_Isend
//mpi_Irecv
//cal1

//sync

//cal2

//tag = 0  send

//MPI wait



//
// Created by cqx820 on 3/27/17.
//
