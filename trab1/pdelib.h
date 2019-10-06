#include <unistd.h>
#include <getopt.h>
// Cabeçalhos das funções da pdelib.c e structs

#define pi 3.14159265359
typedef double t_float;

// Linear System Pentadiagonal Type
typedef struct t_LS5Diag {
    double *main_diagonal,
           *bottom_diagonal,
           *upper_diagonal,
           *away_bottom_diagonal,
           *away_upper_diagonal;
    // Termos independentes
    double *b;
    // Número de equações do SL
    unsigned int n;
} t_LS5Diag;

t_float f(t_float xi, t_float yj);

void allocatePentadiagonalLinearSystem(t_LS5Diag *SL, int nx, int ny);

void gaussSeidel(t_LS5Diag *SL, t_float *x, t_float error, int max_iterations);

void show_help(char *name);