#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    // inicializa gerador de nr aleatoreos
    srand(20192);
    int n = 3;
    SistLinear_t *SL = alocaSistLinear(n);
    real_t *x = malloc(sizeof(real_t)*n);
    inicializaSistLinear(SL, comSolucao, 100);
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++)  {
//            scanf("%lf", &SL->A[i*n+j]);
//        }
//    }
//    for (int i = 0; i < n; i++) {
//        scanf("%lf", &SL->b[i]);
//    }
    printf("\nSistema Linear gerado\n");
    prnSistLinear(SL);
    if (eliminacaoGauss(SL, x, 1) != -1) {
        printf("\nSistema Linear após eliminação de Gauss\n");
        prnSistLinear(SL);
        printf("\nResolução do sistema: \n");
        prnVetor(x, SL->n);
    }
//    replacement(SL, x);
    printf("\nResolução do sistema: \n");
    prnVetor(x, SL->n);


}

