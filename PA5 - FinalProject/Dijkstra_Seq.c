/*
 * Authors: Qixiang Chao, Qiaofeng Wang, Xuanyang Luo
 *
 * Compile:  gcc -g -Wall -o dij_seq Dijkstra_Seq.c
 * Run:      ./dij_seq
 *           Only shows run time
 * OR
 *           ./dij_seq -v
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
#include <sys/time.h>
#include <string.h>

#define INFINITY 1000000
#define MAXLINE 60

void dijkstra(int **G, int n, int startnode, int, int, int);

int main(int argc, char **argv) {

    double elapsedTime;
    struct timeval start, end;
    FILE *graph_matrix;
    char line[MAXLINE];
    int row, column, i, j, flag = 0;
    int **matrix;
    graph_matrix = fopen("mat", "r");
    if (graph_matrix == NULL) {
        printf("Error Reading File\n");
        return -1;
    }
    if (argc >= 2 && strcmp(argv[1], "-v") == 0) flag = 1;

    fgets(line, 30, graph_matrix);
    sscanf(line, "%d %d\n", &row, &column);

    if (row != column) {
        printf("Not a valid adjacency matrix\n");
        return -1;
    }
    matrix = (int **) malloc(sizeof(int *) * row);

    for (i = 0; i < row; i++)
        matrix[i] = (int *) malloc(sizeof(int) * column);


    for (i = 0; i < row; i++)
        for (j = 0; j < column; j++) {
            fscanf(graph_matrix, "%d", &matrix[i][j]);
            //  printf("right!!!!!\n");
        }
    //  printf("right!!!!! %d %d\n", row, column);

    gettimeofday(&start, NULL);
    dijkstra(matrix, row, 0, row, column, flag);
    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
    printf("\nTime: %lf ms.\n", elapsedTime);

    return 0;
}

void dijkstra(int **matrix, int n, int startnode, int row, int column, int flag) {

    int cost[row][column], distance[row], pred[row];
    int visited[row], count, mindistance, nextnode, i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (matrix[i][j] == 0)
                cost[i][j] = INFINITY;
            else
                cost[i][j] = matrix[i][j];

    //initialize pred[],distance[] and visited[]
    for (i = 0; i < n; i++) {
        distance[i] = cost[startnode][i];
        pred[i] = startnode;
        visited[i] = 0;
    }

    distance[startnode] = 0;
    visited[startnode] = 1;
    count = 1;

    while (count < n - 1) {
        mindistance = INFINITY;

        //nextnode gives the node at minimum distance
        for (i = 0; i < n; i++)
            if (distance[i] < mindistance && !visited[i]) {
                mindistance = distance[i];
                nextnode = i;
            }

        //check if a better path exists through nextnode
        visited[nextnode] = 1;
        for (i = 0; i < n; i++)
            if (!visited[i])
                if (mindistance + cost[nextnode][i] < distance[i]) {
                    distance[i] = mindistance + cost[nextnode][i];
                    pred[i] = nextnode;
                }
        count++;
    }

    if (flag == 1) {
        int v, w, *path;

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
            printf("%3d       %4d\n", v, distance[v]);
        printf("\n");

        free(path);
    }
}
