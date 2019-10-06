// implementa as funções de Gauss Seidel, diferenças finitas e utilitários
#include "pdelib.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// Função f(xi, yi)
t_float f(t_float xi, t_float yj) {
    return 4 * pow(pi, 2) * (
        sin(2 * pi * xi) *
        sinh(pi * yj) +
        sin(2 * pi * (pi - xi)) *
        sinh(pi * (pi - yj))
    );
}

void allocatePentadiagonalLinearSystem(t_LS5Diag *SL, int nx, int ny)
{
    SL->away_upper_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 2));
    SL->upper_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 1));
    SL->main_diagonal = malloc(sizeof(t_float) * (nx * ny));
    SL->bottom_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 1));
    SL->away_bottom_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 2));
    // Termos independentes
    SL->b = malloc(sizeof(t_float) * (nx * ny));
    // Número de equações do Sistema Linear
    SL->n = nx*ny;
}

void gaussSeidel(t_LS5Diag *SL, t_float *x, t_float error, int max_iterations)
{
    double norma, diff, xk; int k = 1, i;
    do {
        // primeira equação fora do laço
        i = 0;
        xk = (SL->b[i] - (SL->upper_diagonal[i] * x[i+1]) -
            SL->away_upper_diagonal[i] * x[i+2]
        ) / SL->main_diagonal[i];
        norma = fabs(xk - x[0]);
        x[i] = xk;

        // Segunda equação fora do laço
        i++;
        xk = (SL->b[i] - (
            SL->bottom_diagonal[i-1] * x[i-1]) -
            SL->upper_diagonal[i] * x[i+1] -
            SL->away_upper_diagonal[i] * x[i+2]
        ) / SL->main_diagonal[i];
        norma = fabs(xk - x[1]);
        x[i] = xk;

        for (i = 2; i < (SL->n) - 2; i++) {
            xk = (SL->b[i] - SL->away_bottom_diagonal[i-2] * x[i-2] -
                SL->bottom_diagonal[i-1] * x[i-1] -
                SL->upper_diagonal[i] * x[i+1] -
                SL->away_upper_diagonal[i] * x[i+2]
            ) / SL->main_diagonal[i];
            // Calcula norma || x(k) - x(k - 1) ||
            diff = fabs(xk - x[i]);
            norma = (diff > norma) ? (diff) : (norma);
            x[i] = xk;
        }

        // Penúltima equação fora do laço
        i++;
        xk = (SL->b[i] - SL->away_bottom_diagonal[i-2]*x[i-2] -
            SL->bottom_diagonal[i-1]*x[i-1] -
            SL->upper_diagonal[i]*x[i+1]
        ) / SL->main_diagonal[i];
        diff = fabs(xk - x[i]);
        norma = (diff > norma) ? (diff) : (norma);
        x[i] = xk;

        // ultima equação fora do laço
        i++;
        xk = (SL->b[i] - SL->away_bottom_diagonal[i-2]*x[i-2] -
            SL->bottom_diagonal[i-1] * x[i-1]
        ) / SL->main_diagonal[i];
        diff = fabs(xk - x[i]);
        norma = (diff > norma) ? (diff) : (norma);
        x[i] = xk;

        k++;
    } while (norma > error && k < max_iterations) ;
}


void show_help(char *name) {
    fprintf(stderr, "\
        [uso] %s <opções>\n\
        --nx INTEGER     Número de pontos da malha na direção de x.\n\
        --ny INTEGER     Número de pontos da malha na direção de y.\n\
        --i INTEGER      Número máximo de iterações do método de Gauss-Seidel.\n\
        --o STRING       Nome do arquivo de saída.\n\
        ", name);
    exit(-1);
}