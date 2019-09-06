#include <stdio.h>
#include "SistemasLineares.h"
#include "utils.h"

int main() {
    int n = 3;
    SistLinear_t *SL = alocaSistLinear(n);
    real_t *x = malloc(sizeof(real_t)*n);
    inicializaSistLinear(SL, comSolucao, 100);
    eliminacaoGauss(SL, x, 0);
    prnSistLinear(SL);

    return 0;
}