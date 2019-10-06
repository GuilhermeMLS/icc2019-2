#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdelib.h"

// Argumentos do programa
const struct option stopcoes[] = {
    {"nx", required_argument, 0, 'a'},
    {"ny", required_argument , 0, 'b'},
    {"i", required_argument , 0, 'i'},
    {"o", required_argument, 0, 'o'},
    {0, 0, 0, 0},
};

int main(int argc, char *argv[]) {
    // Quantidade de pontos na discretização
    int nx; // M
    int ny; // N

    int max_iterations;
    t_float error;
    int opt;
    char *output_file;

    // Leitura dos argumentos da linha de comando
    if (argc < 2) {
        show_help(argv[0]);
    }
    while((opt = getopt_long(argc, argv, "ha:b:i:o:", stopcoes, NULL)) > 0) {
        switch (opt) {
            case 'h': /* -h ou --help */
                show_help(argv[0]);
                break ;
            case 'a': /* -a ou --nx */
                nx = atoi(optarg);
                break ;
            case 'b': /* -b ou --nx */
                ny = atoi(optarg);
                break ;
            case 'i': /* -i ou --i */
                max_iterations = atoi(optarg);
                break ;
            case 'o': /* -o ou --o */
                output_file = optarg;
                break ;
            default:
                fprintf(stderr, "Opção inválida ou argumento faltando: `%c'\n", optopt) ;
                return -1 ;
        }
    }

    // Printa as opões que recebeu
    printf("\nnx = %d\n", nx);
    printf("\nny = %d\n", ny);
    printf("\ni = %d\n", max_iterations);
    printf("\no = %s\n", output_file);

    /******************************************/

    // Comprimento dos lados
    t_float lx = pi;
    t_float ly = pi;
    // Espaço entre dois pontos da malha
    t_float hx = lx / nx + 1;
    t_float hy = ly / ny + 1;
    // Sistema Linear
    t_LS5Diag SL;
    allocatePentadiagonalLinearSystem(&SL, nx, ny);

    /* Calcula os valores das diagonais do Sistema Linear */

    // Diagonal superior afastada
    for (int k = 0; k < (nx*ny)-2; k++) {
        SL.away_upper_diagonal[k] = (-1 * pow(hx, 2)) * (2 - hy);
    }
    // Diagonal superior
    for (int k = 0; k < (nx*ny)-1; k++) {
        SL.upper_diagonal[k] = (-1 * pow(hy, 2)) * (2 - hx);
    }
    // Diagonal principal
    for (int k = 0; k < nx*ny; k++) {
        SL.main_diagonal[k] = 4 * (pow(hx, 2) + pow(hy, 2) + 2 * pow(pi, 2) * pow(hx, 2) * pow(hy, 2));
    }
    // Diagonal inferior
    for (int k = 0; k < (nx*ny)-1; k++) {
        SL.bottom_diagonal[k] = (-1 * pow(hy, 2)) * (2 + hx);
    }
    // Diagonal inferior afastada
    for (int k = 0; k < (nx*ny)-2; k++) {
        SL.away_bottom_diagonal[k] = (-1 * pow(hx, 2)) * (2 + hy);
    }

    /* Calcula o valor dos termos independentes */

    int k = 0;

    for (int j = 0; j < nx; j++) {
        for (int i = 0; i < ny; i++) {
            k++;
            SL.b[k] = pow(hx, 2) * pow(hy, 2) * f(i * hx, j * hy);
            if (i == 1) {
                SL.b[k] -= pow(hy, 2) * 0;
                if (j > 1) {
                    SL.bottom_diagonal[k-1] = 0;
                }
            }
            if (i == ny) {
                SL.b[k] -= pow(hy, 2) * 0;
                if (j < nx) {
                    SL.upper_diagonal[k+1] = 0;
                }
            }
            if (j == 1) {
                SL.b[k] -= pow(hx, 2) * sin(2 * pi * (pi - (i * hx)) * sinh(pow(pi, 2)));
            }
            if (j == nx) {
                SL.b[k] -= pow(hx, 2) * sin(2 * pi * (i * hx) * sinh(pow(pi, 2)));
            }
        }
    }

    /* Printar essa caralha toda agora */

    printf("\n");
    for(int i = 0; i < (nx * ny) - 2; i++) {
        printf("%lf ", SL.away_upper_diagonal[i]);
    }
    printf("\n");
    for(int i = 0; i < (nx * ny) - 1; i++) {
        printf("%lf ", SL.upper_diagonal[i]);
    }
    printf("\n");
    for(int i = 0; i < (nx * ny); i++) {
        printf("%lf ", SL.main_diagonal[i]);
    }
    printf("\n");
    for(int i = 0; i < (nx * ny) - 1; i++) {
        printf("%lf ", SL.bottom_diagonal[i]);
    }
    printf("\n");
    for(int i = 0; i < (nx * ny) - 2; i++) {
        printf("%lf ", SL.away_bottom_diagonal[i]);
    }
    printf("\n\n## Vetor B ##\n\n");
    printf("B = ");
    for(int i = 0; i < SL.n; i++) {
        printf("%lf ", SL.b[i]);
    }
    printf("\n");

    /* Jogar o Sistema Linear no Gauss-Seidel */
    t_float *x = calloc(SL.n, sizeof(t_float)); // Vetor solução
    gaussSeidel(&SL, x, error, max_iterations);

    printf("\n ## Vetor Solução ## \n");
    printf("X = ");
    for(int i = 0; i < SL.n; i++) {
        printf("%lf ", x[i]);
    }
    printf("\n\n");

    return 0;
}