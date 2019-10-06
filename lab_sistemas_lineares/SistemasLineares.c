#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Verifica se existe equação nula no sistema SL

  \param SL Ponteiro para o sistema linear
*/

int temEquacaoNula(SistLinear_t *SL) {
    int n = SL->n;
    real_t *matrix = SL->A;
    real_t *b = SL->b;
    int eq_nula = 0;

    for (int i = 0; i < n; i++) {
        if (matrix[i*n] == 0) {
            eq_nula = 1;
            for (int j = 0; j < n; j++) {
                if (matrix[i*n+j] != 0) {
                    eq_nula = 0;
                }
            }
        }
        if (eq_nula && b[i] == 0) {
            return 1;
        }
    }

    return 0;
}

/*!
  \brief Essa função executa a retrosubstituição

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
void replacement(SistLinear_t *SL, real_t *x) {
    int n = SL->n;
    real_t *matrix = SL->A;
    real_t *b = SL->b;

    x[n-1] =  b[n-1] / matrix[(n-1)*n+(n-1)];

    real_t sum = 0;

    for (int i = n-2; i >= 0; i--) {
        sum = b[i];
        for (int j = i+1; j < n; j++) {
            sum -= matrix[i*n+j] * x[j];
        }
        x[i] = sum / matrix[i*n+i];
    }
}

/*!
  \brief Essa função executa o pivotamento da eleiminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param k Iteração atual da eliminação de Gauss
*/
void pivotar(SistLinear_t *SL, int k)
{
    int n = SL->n;
    real_t *matrix = SL->A;
    real_t *b = SL->b;

    // 1) Encontrar maior elemento da coluna abaixo do pivô
    int column = k;
    int greater_line_index = k;
    real_t greater_element = matrix[k*n+k];
    for(int i = k; i < n; i++) {
        if (fabs(matrix[i*n+column]) > fabs(greater_element)) {
            greater_line_index = i;
        }
    }

    // 2) Trocar linha a com linha b (b seria a do maior elemento)
    int linha1 = k;
    int linha2 = greater_line_index;

    // 2.1) Salva a linha1 em um array auxiliar
    real_t *aux = malloc(sizeof(real_t)*n);
    real_t auxb;
    for (int j = 0; j < n; j++) {
        aux[j] = matrix[linha1*n+j];
    }
    // 2.2) Insere linha2 em linha1
    for (int j = 0; j < n; j++) {
        matrix[linha1*n+j] = matrix[linha2*n+j];

    }
    // 2.3) Insere linha aux (linha1) em linha2
    for (int j = 0; j < n; j++) {
        matrix[linha2*n+j] = aux[j];

    }
    // Troca também as posições em B
    auxb =  b[linha1];
    b[linha1] = b[linha2];
    b[linha2] = auxb;

    free(aux);
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
double normaL2Residuo(SistLinear_t *SL, real_t *x)
{

}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param pivotamento flag para indicar se o pivotamento parcial deve ser feito (!=0)

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, int pivotamento)
{
    // Verifica se há equações nulas
    if (temEquacaoNula(SL)) {
        fprintf(stderr, "Equação nula detectada, os sitema possui infinitas soluções.\n");
        return -1;
    }

    // Acessar matriz por vetor: i*n+j
    int n = SL->n;
    real_t *matrix = SL->A;
    real_t *b = SL->b;
    real_t m = 0;

    // Eliminação de Gauss
    for (int k = 0; k < n-1; k++) {
        // Pivotamento
        if (pivotamento) {
            pivotar(SL, k);
        }
        for (int i = 1 + k; i < n; i++) {
            m = matrix[i*n+k] / matrix[k*n+k];
            for (int j = k; j < n; j++) {
                if (j == k) {
                    matrix[i*n+j] = 0;
                } else {
                    matrix[i*n+j] = matrix[i*n+j] - matrix[k*n+j] * m;
                }
            }
            b[i] = b[i] - b[k] * m;
            printf("\ndebug\n");
            prnSistLinear(SL);
            printf("\n------\n");
        }
    }

    // Resolução do sistema
    replacement(SL, x);

    //TODO: tentar saber se o sistema foi resolvido ou não
    // talvez com os esquemas de erros, ler a ativade no Moodle

    return 0;
}

/*!
  \brief Método de Gauss-Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, real_t erro)
{
    int n = SL->n;
    real_t *b = SL->b;
    real_t *A = SL->b;

    // 1) Isola os elementos das diagonais principais
    // 2) Chuta um vetor inicial de solução
    // 2) Enquanto o erro não estiver satisfatoriamente pequeno
        // Calcular vetor solução a partir do vetor anterior
    real_t *kick = malloc(sizeof(real_t)*n);
    for (int i = 0; i < n; i++) {
        kick[i] = 0;
    }
    // TODO: Ao final desse for eu vou ter a primeira iteração de X em cima do chute inicial
    //     Achar um jeito de salvar ela e usar na próxima iteração
    //     Até que o erro entre a atual e a anterior seja bem pequeno
    for (int i = 0; i < n; i++) {
        kick[i] = b[i] / A[i*n+i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                kick[i] -= (A[i*n+j] * kick[j]) / A[i*n+i];
            }

        }
    }

}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro)
{


}


// Alocaçao de memória
SistLinear_t* alocaSistLinear (unsigned int tam)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  if ( SL ) {
    SL->A = (real_t *) malloc(tam * tam * sizeof(real_t));
    SL->b = (real_t *) malloc(tam * sizeof(real_t));

    if (!(SL->A) || !(SL->b))
      liberaSistLinear(SL);
  }
  
  SL->n = tam;

  return (SL);
}

// Liberacao de memória
void liberaSistLinear (SistLinear_t *SL)
{
  free(SL->A);
  free(SL->b);
  free(SL);
}

/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: comSolucao,
  eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void inicializaSistLinear (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max)
{
  unsigned int tam = SL->n;
  // para gerar valores no intervalo [0,coef_max[
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // inicializa vetor b
  for (unsigned int i=0; i<tam; ++i) {
    SL->b[i] = (real_t)rand() * invRandMax;
  }
    
  if (tipo == hilbert) {
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[i*tam+j] = 1.0 / (real_t)(i+j+1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[i*tam+j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[nula*tam+j] = 0.0;
      }
      SL->b[nula] = 0.0;
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % tam;
      unsigned int propSrc = (propDst + 1) % tam;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[propDst*tam+j] = SL->A[propSrc*tam+j] * mult;
      }
      SL->b[propDst] = SL->b[propSrc] * mult;
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % tam;
      unsigned int combSrc1 = (combDst + 1) % tam;
      unsigned int combSrc2 = (combDst + 2) % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[combDst*tam+j] = SL->A[combSrc1*tam+j] + SL->A[combSrc2*tam+j];
      }
      SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
    }
    else if (tipo == diagDominante) {
      // aumenta o expoente dos termos da diagonal principal
      for (unsigned int i=0; i<tam; ++i) {
        SL->A[i*tam+i] *= (real_t)tam;
      }
    }

  }
}


SistLinear_t *lerSistLinear ()
{
  unsigned int n;
  SistLinear_t *SL;
  
  scanf("%d",&n);

  SL = alocaSistLinear (n);
  
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      scanf ("%lg", &SL->A[i*n+j]);

  for(int i=0; i < n; ++i)
    scanf ("%lg", &SL->b[i]);
  
  return SL;
}


void prnSistLinear (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n\t");
    for(int j=0; j < n; ++j)
      printf ("%10.5lg", SL->A[i*n+j]);
    printf ("   |   %.8lg", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10.10lg ", v[i]);
  printf ("\n\n");

}




