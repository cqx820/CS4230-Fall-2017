/*
* Authors: Qixiang Chao, Qiaofeng Wang, Xuanyang Luo
* Compile: gcc -o generator adj_mat_gen.c
* Run: ./generator 5 5
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    FILE *file;
    int i, j, row, column;

    row = atoi(argv[1]);
    column = atoi((argv[2]));

    if (row != column) return -1;

    if ((file = fopen("mat", "w")) == NULL) puts("Errrr");

    fprintf(file, "%d %d\n", row, column);
    srand(time(NULL));
    int count = 0;
    for (i = 0; i < row; i++) {
        for (j = 0; j < column; j++) {
            if (j == count) fprintf(file, "%d ", 0);
            else fprintf(file, "%d ", rand() % 21);
        }
        fprintf(file, "\n");
        count++;
    }

    return 0;
}
