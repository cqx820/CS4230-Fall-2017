
/*
 * Authors: Qixiang Chao, Qiaofeng Wang, Xuanyang Luo
 *
 * Compile:  mpicc -g -Wall -o dij Dijkstra_MPI
 * Run:      mpiexec -n <p> ./dij
 *           Only shows run time
 * OR
 *           mpiexec -n <p> ./dij -v
 *           Shows both result and run time
 * (row number must equal to column number in the first line of the file and the number of processes p, should evenly divide n)
 *
 * Sequential code refered from:
 * http://www.thecrazyprogrammer.com/2014/03/dijkstra-algorithm-for-finding-shortest-path-of-a-graph.html
 *
 * Input:
 * mat(the matrix file)
 *
 * Result:
 * The path found from start node to each vertex.
 * The length of a shortest path from start vertex to every vertex
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>
//#include "/usr/include/mpi/mpi.h"
#include <string.h>

#define INFINITY 1000000
#define MAXLINE 60

static MPI_Datatype build_col_block_type(int, int);

static void Dijkstra(int *, int *, int *, int, int, int, MPI_Comm);

static void print_results(int *, int *, int, int, int, MPI_Comm);

int n;
//int startnode;

int main(int argc, char **argv) {
    FILE *graph_matrix;
    double elapsedTime;
    struct timeval start, end;

    char line[MAXLINE];
    int p, rank, *matrix, *sink, *pred, *splitted_matrix, i, row, column, flag = 0;

    MPI_Datatype block_column_ty;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n = p;

    // printf("%d", n);

//    if (!strcmp(argv[1], "-t")) //This for taking matrix files
//    {
    graph_matrix = fopen("mat", "r");
    if (graph_matrix == NULL) {
        printf("Error Reading File\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }

    if (argc == 2 && strcmp(argv[1], "-v") == 0) flag = 1;

//    printf("!!!!!!!!!!%d", flag);

    fgets(line, 30, graph_matrix);
//    if (argc != 3) {
//        MPI_Abort(MPI_COMM_WORLD, 1);
//        return -1;
//    }
    if (rank == 0) {
        sscanf(line, "%d %d\n", &row, &column);
//        printf("!!!!!!!!!!!1%d", row);
//        row = atoi(argv[1]);
//        column = atoi(argv[2]);
        if (row != column) {
            printf("Not a valid adjacency matrix\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
    }

    MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);

    sink = (int *) malloc(sizeof(int) * (row / p));
    pred = (int *) malloc(sizeof(int) * (row / p));
    splitted_matrix = (int *) malloc(sizeof(int) * (row / p) * row);

    block_column_ty = build_col_block_type(row, row / p);
    if (rank == 0) {
        matrix = (int *) malloc(sizeof(int) * row * column);
        for (i = 0; i < row * column; i++) {
            fscanf(graph_matrix, "%d", &matrix[i]);
        }
    }
//    } else {
//        if (rank == 0)
//            scanf("%d", &row);
//        MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//        matrix = (int *) malloc(sizeof(int) * row * row);
//        sink = (int *) malloc(sizeof(int) * (row / p));
//        pred = (int *) malloc(sizeof(int) * (row / p));
//        splitted_matrix = (int *) malloc(sizeof(int) * (row / p) * row);
//
//        block_column_ty = build_col_block_type(row, row / p);
//
//        if (rank == 0) {
//            for (i = 0; i < row; i++)
//                for (j = 0; j < row; j++)
//                    scanf("%d", &matrix[i * row + j]);
//        }
//    }

    MPI_Scatter(matrix, 1, block_column_ty, splitted_matrix, row * (row / p), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) free(matrix);

    gettimeofday(&start, NULL);
    Dijkstra(splitted_matrix, sink, pred, row, row / p, rank, MPI_COMM_WORLD);
    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;

    // Print_dists(sink, row, row / p, rank, MPI_COMM_WORLD);
    if (flag == 1 && rank == 0) {
        print_results(pred, sink, row, row / p, rank, MPI_COMM_WORLD);
        printf("\nTime: %lf ms.\n", elapsedTime);
    } else if (flag == 0 && rank == 0)
        printf("\nTime: %lf ms.\n", elapsedTime);

// Print_matrix(loc_mat, n, loc_n, blk_col_mpi_t, rank, comm);
//

    free(splitted_matrix);
    free(sink);
    free(pred);


    /* When you're done with the MPI_Datatype, free it */
    MPI_Type_free(&block_column_ty);


    MPI_Finalize();
    return 0;
}

static MPI_Datatype build_col_block_type(int n, int splittedN) {
    //MPI_Datatype lower_bound, ex;
    MPI_Aint lb, extent;
    MPI_Datatype lb_extent_ty, block_column_ty, block_f_scatter;

    //splittedN = n/p
    MPI_Type_contiguous(splittedN, MPI_INT, &lb_extent_ty);//constructor for extent and lower bound

//    MPI_Type_lb(lower_bound, &lb);//this routine is deprecated as of MPI-2
//    MPI_Type_extent(ex, &extent); extent = upper bound - lower bound


//                        (min(j) disp(j)                          if no entry
//        lb(Typemap) =   (                                        basic type lb
//                        (min(j) {disp(j) such that type(j) = lb} otherwise

    MPI_Type_get_extent(lb_extent_ty, &lb, &extent);

    MPI_Type_vector(n, splittedN, n, MPI_INT,
                    &block_column_ty);//n blocks with n/p copies each of the MPI_INT, with a stride of n elements between the blocks.

    MPI_Type_create_resized(block_column_ty, lb, extent, &block_f_scatter);

//    MPI_Type_create_resized(block_column_ty, 0, 2 * sizeof(int), &block_f_scatter); //ok to do this?


    MPI_Type_commit(&block_f_scatter);

    MPI_Type_free(&lb_extent_ty);
    MPI_Type_free(&block_column_ty);

    return block_f_scatter;
}

static void Dijkstra(int *matrix, int *distance, int *precedence, int n, int splittedN, int rank, MPI_Comm comm) {
    int i, u, temp_u, temp_dis;
    int *local_min_dis = (int *) malloc(sizeof(int) * 2);
    int *min_distance = (int *) malloc(sizeof(int) * 2);
    int *visited = (int *) malloc(sizeof(int) * splittedN);


    // startnode = 0;

    for (i = 0; i < splittedN; i++) {
//        distance[i]=cost[startnode][i];
//        pred[i]=startnode;
//        visited[i]=0;
        distance[i] = matrix[i];
        precedence[i] = 0;
        visited[i] = 0;
    }
    if (rank == 0) visited[0] = 1;

    for (i = 1; i < n; i++) {
        int j, next_node = -1, mindistance = INFINITY;
        //      omp_set_num_threads(n);
//#pragma omp parallel for schedule(static, 1024)
        for (j = 0; j < splittedN; j++) {
            if (!visited[j])
                if (distance[j] < mindistance) {
                    next_node = j;
                    mindistance = distance[j];
                }
        }

        if (next_node >= 0) {
            local_min_dis[0] = distance[next_node];
            local_min_dis[1] = next_node + rank * splittedN;
        } else {
            local_min_dis[0] = INFINITY;
            local_min_dis[1] = -1;
        }

        //compute a global minimum of a pair of int
        MPI_Allreduce(local_min_dis, min_distance, 1, MPI_2INT, MPI_MINLOC, comm);


        temp_u = min_distance[0];
        u = min_distance[1];

        if (u / splittedN == rank) visited[next_node] = 1;

//        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%d\n", n);
        //       omp_set_num_threads(n);
//#pragma omp parallel for schedule(static, 1024)
        for (j = 0; j < splittedN; j++) {
            if (!visited[j]) {

                temp_dis = temp_u + matrix[u * splittedN + j];

                if (temp_dis < distance[j]) {
                    distance[j] = temp_dis;
                    precedence[j] = u;
                }
            }
        }
    }
    if (rank == 0) {
        free(local_min_dis);
        free(min_distance);
        free(visited);
    }

}

static void print_results(int *precedence, int *distance, int n, int splittedN, int rank, MPI_Comm comm) {
    int *pred, v, w, *path, count, i, *dist;

    if (rank == 0) {
        pred = (int *) malloc(n * sizeof(int));
        dist = (int *) malloc(n * sizeof(int));
    }

    MPI_Gather(precedence, splittedN, MPI_INT, pred, splittedN, MPI_INT, 0, comm);
    MPI_Gather(distance, splittedN, MPI_INT, dist, splittedN, MPI_INT, 0, comm);

    if (rank == 0) {

        path = (int *) malloc(n * sizeof(int));

        printf("  v     Path 0->v\n");
        printf("----    ---------\n");
        for (v = 1; v < n; v++) {
            printf("%3d:    ", v);
            count = 0;
            w = v;
            while (w != 0) {
                path[count] = w;
                count++;
                w = pred[w];
            }
            printf("0 ");
            for (i = count - 1; i >= 0; i--)
                printf("%d ", path[i]);
            printf("\n");
        }
        printf("\n");
        printf("  v    dist 0->v\n");
        printf("----   ---------\n");

        for (v = 1; v < n; v++)
            printf("%3d       %4d\n", v, dist[v]);
        printf("\n");

        free(path);
        free(pred);
        free(dist);
    }
}

