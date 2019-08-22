#include <stdio.h>
#include "stdlib.h"

#include "libs/matrix.h"

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("\nInvalid parameters\n");
        return 0;
    }

    printf("\n## Solving Linear Systems with Replacement Technique ##\n");
    int n = atof(argv[1]);
    float *b = malloc(sizeof(float) * n);

    float *array_result = malloc(sizeof(float) * n);

    // Allocate matrix
    float **matrix = allocateMatrix(n);

    printf("\nN = %s\n", argv[1]);

    printf("\nLoading matrix A\n");
    loadMatrix(matrix, n);

    printMatrix(matrix, n);

    printf("\nLoading matrix B\n");
    for (int i = 0; i < n; i++) {
        scanf("%f", &b[i]);
    }

    printArray(b, n);

    replacementTechnique(matrix, b, array_result, n);

    printf("\nResult matrix:\n");
    printArray(array_result, n);

    return 0;
}


