#include <stdio.h>
#include "libs/matrix.h"

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("\nInvalid parameters\n");
        return 0;
    }
    printf("\n## Gauss Elimination ##\n");
    int n = atof(argv[1]);
    float *b = malloc(sizeof(float) * n);

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

    gaussElimination(matrix, b, n);

    printf("\nResulting matrix A:\n");
    printMatrix(matrix, n);

    printf("\nResulting matrix B:\n");
    printArray(b, n);

    float *array_result = malloc(sizeof(float) * n);

    printf("\nRunning replacement technique...\n");
    replacementTechnique(matrix, b, array_result, n);

    printf("\nResulting matrix:\n");
    printArray(array_result, n);

    return 0;
}