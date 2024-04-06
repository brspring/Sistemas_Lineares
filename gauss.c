#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXIT 50
typedef double rtime_t;

rtime_t timestamp (void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return ( (rtime_t) tp.tv_sec*1.0e3 + (rtime_t) tp.tv_nsec*1.0e-6 );
}

double **alocaMatriz(int ordem)
{
    double **matriz;
    matriz = (double **)calloc(ordem, sizeof(double *));
    for (int i = 0; i < ordem; i++)
        matriz[i] = (double *)calloc(ordem, sizeof(double));
    return matriz;
}

void desalocaMatriz(double **matriz, int ordem)
{
    for (int i = 0; i < ordem; i++)
        free(matriz[i]);
    free(matriz);
}

void separaTridiagonais(double **Matriz, double *a, double *d, double *c, int n){
    for(int i = 0; i < n; ++i){
            d[i] = Matriz[i][i];
            if(i < n-1){
                a[i] = Matriz[i][i+1];
                c[i] = Matriz[i+1][i];
            }   
    }
}

int encontraMax(double **A, int i, uint n) {
    double maxValor = fabs(A[i][0]);
    int indiceMax = 0;

    for (int j = 1; j < n; j++) {
        if (fabs(A[i][j]) > maxValor) {
            maxValor = fabs(A[i][j]);
            indiceMax = j;
        }
    }

    return indiceMax;
}

void residuoEliminacaoDeGaussTriDiagonais(double *a, double *d, double *c, double *b, double *x, double *residuo, int n) {
    //calcula Ax
    for (int i = 0; i < n; ++i) {
        residuo[i] = d[i] * x[i];
        if (i > 0) {
            residuo[i] += a[i - 1] * x[i - 1];
        }
        if (i < n - 1) {
            residuo[i] += c[i] * x[i + 1];
        }
    }

    //subtrai b
    for (int i = 0; i < n; ++i) {
        residuo[i] -= b[i];
    }
}

double residuoMatriz(double **A, double *x, double *b, double *residuo, int n) {
    double *Ax = (double *)malloc(n * sizeof(double));

    //calcula Ax
    for (int i = 0; i < n; ++i) {
        Ax[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            Ax[i] += A[i][j] * x[j];
        }
        residuo[i] = Ax[i] - b[i];
    }

     //libera memoria alocada
    free(Ax);
}

void trocaLinhas(double **A, double *b, int i, int iPivo){
    double *aux = A[i];
    A[i] = A[iPivo];
    A[iPivo] = aux;

    double auxB = b[i];
    b[i] = b[iPivo];
    b[iPivo] = auxB;
}

void eliminacaoDeGauss(double **A, double *b, double *x, uint n){
    //triangularizacao
    for(int i=0; i<n; ++i){
        uint iPivo = encontraMax(A, i, n);
        if(i != iPivo)
            trocaLinhas(A, b, i, iPivo);

        for(int k=i+1; k < n; ++k){
            double m = A[k][i]/A[i][i];
            A[k][i] = 0.00;
            for(int j=i+1; j<n; ++j){
                A[k][j] -= A[i][j] * m;
            }
            b[k] -= m * b[i];
        }        
    //zerando o vetor x
    for(int i = 0; i < n; ++i)
        x[i] = 0.0;

    //retrosubstituicao
    for(int i = n-1; i >= 0; --i){
        x[i] = b [i];
        for(int j = i+1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    }
}

double max(double a, double b) {
    return (a > b) ? a : b;
}

void gaussSeidel(double **A, double *b, double *x, uint n, double tol){
    double erro = 1.0 + tol;
    int j, s, count = 0;
    double *x_antigo = malloc(n * sizeof(double));

    while(erro > tol && count < MAXIT){
        x_antigo = x;

        for(int i = 0; i < n; i++){
            x_antigo[i] = x[i];
            for(j=0, s = 0; j < n; ++j){
                if(i != j)
                    s += A[i][j] * x[j];
            }
            x[i] = (b[i] - s)/A[i][i];
            count++;
            //calcula o erro
            //erro = max(fabs(x[i] - x_antigo[i]));

        }
        //erro = max(fabs(x[i] - x_antigo[i]));
        //print vetor x
        for(int i = 0; i < n; ++i)
            printf("x[%d] = %.2f\n", i, x[i]);
    }
    
    //free(x_antigo);
}
// void gaussSeidelTriDiagonais(double *d, double *a, double *c, double *b, double *x, uint n, double tol){
//     double erro = 1.0 + tol;
//     int j, s;
//     while(erro < tol){
//         for(int i = 0; i < n; i++)
//             for(s=0, j=0; j < n; ++j)
//                 if(i != j)
//                     s += A[i][j] * x[j];

//     x[i] = (b[i] - s)/A[i][i];
//     }
// }
void eliminacaoDeGaussTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n){
    //triangularizacao
    for(int i = 0; i < n-1; ++i){
            int k = i+1;
            double m = a[i]/d[i];
            a[i] = 0;
            d[i+1] -= c[i] * m;
            b[i+1] -= b[i] * m;     
    }
    //retrosubstituicao
    x[n-1] = b[n-1] / d[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (b[i] - c[i] * x[i+1]) / d [i];
    }
}

int main(){
    double *a, *b, *c, *d, *x;
    double *residuo;
    int n;
    double tol;

    //leitura dimensoes da mateiz
    scanf("%d", &n);

   // alocando memoria para a matriz
    double **Matriz = alocaMatriz(n);
    
    //alocando memoria para o vetor b
    b = (double *)malloc(n * sizeof(double));
    //alocando memoria vetor residuo
    residuo = (double *)malloc(n * sizeof(double));
    //alocando memória para os vetores da matriz tridiagonal
    d = (double *)malloc(n * sizeof(double));
    a = (double *)malloc(n * sizeof(double)); 
    c = (double *)malloc(n*  sizeof(double));
    x = (double *)malloc(n*  sizeof(double));

    //lendo a matriz e o vetor b
    for(int i=0; i < n; ++i){
        Matriz[i] = (double *)malloc((n + 1) * sizeof(double));
        for(int j=0; j <= n; ++j){
            if(j==n){
                scanf("%lf", &b[i]);
            }else{
                scanf("%lf", &Matriz[i][j]);
            }
        }
    }

    // vetores matriz tridiagonais 
    separaTridiagonais(Matriz, a, d, c, n);

    double tempo;
    tol = 0.0001;
    tempo = timestamp();
    //eliminacaoDeGauss(Matriz, b, x, n);
    //residuoMatriz(Matriz, x, b, residuo, n);
    gaussSeidel(Matriz, b, x, n, tol);
    //eliminacaoDeGaussTriDiagonais(a, d, c, b, x, n);
    //residuoEliminacaoDeGaussTriDiagonais(a, d, c, b, x, residuo, n);
    tempo = timestamp() - tempo;

    /*printf("EG clássico:\n");
    printf("%.8f ms\n", tempo);
    for (int i = 0; i < n; ++i)
        printf("%lf   ", x[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuo[i]);

    printf("\n\n");*/
    //mostra a matriz lida fazer funcao
    printf("------- Matriz resultante -------\n");
    for(int i=0; i<n; ++i){
        for(int j=0; j<=n; ++j){
            if (j == n )
                printf("| %.2f", b[i]);
            else
                printf("%.2f ", Matriz[i][j]);
        }
        printf("\n");
    }
    for(int i=0; i<n; ++i)
        printf("x[%d] = %.2f\n", i, x[i]);

    desalocaMatriz(Matriz, n);
    free(a);
    free(b);
    free(d);
    free(c);
    free(x);
    free(residuo);

    return 0;
}