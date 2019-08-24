#include "matrix.h"
#include <stdio.h>
#include "stdlib.h"

void replacementTechnique(float **matrix, float *b, float *x, int n) {
    x[n-1] =  b[n-1] / matrix[n-1][n-1];

    float sum = 0;

    for ( int i = n-2; i >= 0; i-- ) {
        sum = b[i];
        for ( int j = i+1; j < n; j++ ) {
            sum -= matrix[i][j] * x[j];
        }
        x[i] = sum / matrix[i][i];
    }
}

void gaussElimination(float **matrix, float *b, int n) {
    float m = 0;
    for (int k = 0; k < n-1; k++) {
        for (int i = 1 + k; i < n; i++) {
            m = matrix[i][k] / matrix[0][k];
            for (int j = k; j < n; j++) {
                if (i ==  n-1) {
                    // tá na borda
                    matrix[i][j] = matrix[i][j] - matrix[0][j] * m;
                } else {
                    // não tá na borda
                    matrix[i + k][j + k] = matrix[i + k][j + k] - matrix[0][j + k] * m;
                }
            }
            b[i] = b[i] - b[0] * m;
        }
    }
}

void printMatrix(float **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %f ", matrix[i][j]);
            if (j == n-1) {
                printf("\n");
            }
        }
    }
}

float **allocateMatrix(int n) {
    float **matrix = malloc(sizeof(float) * n * n);
    if (matrix != NULL) {
        /*Allocating space for each row*/
        for (int i = 0; i < n; i++) {
            matrix[i] = malloc(sizeof(float) * n);
        }
    }
    return matrix;
}

void printArray(float *arr, int n) {
    printf("[");
    for (int i = 0; i < n; i++) {
        printf(" %f ", arr[i]);
    }
    printf("]\n");
}

void loadMatrix(float **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%f", &matrix[i][j]);
        }
    }
}