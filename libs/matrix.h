#ifndef __MATRIXES
#define __MATRIXES

#include <stdio.h>
#include "stdlib.h"

void replacementTechnique(float **matrix, float *b, float *x, int n);
void printMatrix(float **matrix, int n);
float **allocateMatrix(int n);
void printArray(float *arr, int n);
void loadMatrix(float **matrix, int n);
#endif