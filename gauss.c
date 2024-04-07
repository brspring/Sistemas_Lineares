#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "utils.h"

#define MAXIT 50

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

void eliminacaoDeGauss(double **A, double *b, double *x, int n){
    double tempo = timestamp();
    double *residuo = malloc(n * sizeof(double));

    //triangularizacao
    for(int i=0; i<n; ++i){
        int iPivo = encontraMax(A, i, n);
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
    tempo = timestamp() - tempo;
    residuoMatriz(A, x, b, residuo, n);

    printf("EG clássico:\n");
    printf("%.8f ms\n", tempo);
    for (int i = 0; i < n; ++i)
        printf("%.12f   ", x[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuo[i]);
    printf("\n\n");

    free(residuo);
}

void gaussSeidel(double **A, double *b, double *x, int n, double tol){
    double tempo = timestamp();
    double erro = 1.0;
    double s;
    int j, count = 0;
    double *x_antigo = malloc(n * sizeof(double));
    double *residuo = malloc(n * sizeof(double));

    while(erro > tol && count < MAXIT ){
        for(int i = 0; i < n; i++){
            x_antigo[i] = x[i];
            s = 0.0;
            for(j=0; j < n; ++j){
                if(i != j)
                    s = (A[i][j] * x[j]) + s;
            }
            x[i] = (b[i] - s)/A[i][i];
        }

        if(count == 0)
            erro = 1.0;
        else
            erro = encontrarMaiorSubtracao(x, x_antigo, n);
        count++;
    }
    tempo = timestamp() - tempo;
    residuoMatriz(A, x, b, residuo, n);

    printf("GS Clássico [%d]:\n", count);
    printf("%.8f ms\n", tempo);
    for (int i = 0; i < n; ++i) {
        printf("%.12f  ", x[i]);
    }
    printf("\n");
    for (int i = 0; i < n; ++i) {
        printf("%.12f  ", residuo[i]);
    }
    printf("\n\n");

    free(x_antigo);
    free(residuo);
}

void gaussSeidelTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n, double tol){
    double tempo = timestamp();
    double erro = 1.0;
    int j, s;
    double *residuo = malloc(n * sizeof(double));
    double *x_antigo = malloc(n * sizeof(double));
    int count = 0;

    while(erro > tol && count < MAXIT ){
            s = 0.0;

            x[0] = (b[0] - (c[0] * x[1])) / d[0];

            for(int i = 0; i < n; ++i){
                x_antigo[i] = x[i];
            }
            for (int i=1; i < n-1; ++i){
                x[i] = (b[i] - a[i-1] * x[i-1] - c[i] * x[i+1]) / d[i];
            }

            x[n-1] = (b[n-1] - a[n-2] * x[n-2]) / d[n-1];

        if(count == 0)
            erro = 1.0;
        else
            erro = encontrarMaiorSubtracao(x, x_antigo, n);
        count++;
    }
    tempo = timestamp() - tempo;
    residuoEliminacaoDeGaussTriDiagonais(a, d, c, b, x, residuo, n);

    printf("GS 3-diagonal[%d]:\n", count);
    printf("%.8f ms\n", tempo);
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", x[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuo[i]);
    printf("\n\n");
}

void eliminacaoDeGaussTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n){
    double tempo = timestamp();
    double *residuo = alocaVetor(n);
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
    tempo = timestamp() - tempo;
    residuoEliminacaoDeGaussTriDiagonais(a, d, c, b, x, residuo, n);

    printf("EG 3-diagonal:\n");
    printf("%.8f ms\n", tempo);
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", x[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuo[i]);
    printf("\n\n");

    free(residuo);
}

int main(){
    double *a, *b, *c, *d,  
    *resultadoEG, *resultadoGS, *resultado3EG, *resultado3GS,
    *dS, *aS, *cS;
    double *residuo;
    int n;
    double tolerancia;

    //leitura dimensoes da matriz
    scanf("%d", &n);

   // alocando memoria para as matrizes
    double **Matriz = alocaMatriz(n);
    double **MatrizGS = alocaMatriz(n);
    double **MatrizEG = alocaMatriz(n);

    //alocando memoria para o vetor b
    b = alocaVetor(n);

    //alocando memoria vetor residuo
    residuo = alocaVetor(n);

    //alocando memória para os vetores da matriz tridiagonal
    d = alocaVetor(n);
    a = alocaVetor(n); 
    c = alocaVetor(n);
    dS = alocaVetor(n);
    aS = alocaVetor(n); 
    cS = alocaVetor(n);

    //alocando vetores resultado de cada metodo
    resultadoEG = alocaVetor(n);
    resultadoGS = alocaVetor(n);
    resultado3EG = alocaVetor(n);
    resultado3GS = alocaVetor(n);

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

    copiaMatriz(Matriz, MatrizEG, n);
    copiaMatriz(Matriz, MatrizGS, n);

    // vetores matriz tridiagonais 
    separaTridiagonais(Matriz, a, d, c, n);
    separaTridiagonais(Matriz, aS, dS, cS, n);

    tolerancia = 0.0001;
    
    eliminacaoDeGauss(MatrizEG, b, resultadoEG, n);
    gaussSeidel(MatrizGS, b, resultadoGS, n, tolerancia);
    eliminacaoDeGaussTriDiagonais(a, d, c, b, resultado3EG, n);
    gaussSeidelTriDiagonais(aS, dS, cS, b, resultado3GS, n, tolerancia);

    desalocaMatriz(Matriz, n);
    desalocaMatriz(MatrizGS, n);
    desalocaMatriz(MatrizEG, n);
    free(a);
    free(resultado3EG);
    free(resultado3GS);
    free(resultadoEG);
    free(resultadoGS);
    free(dS);
    free(aS);
    free(cS);
    free(b);
    free(d);
    free(c);
    free(residuo);

    return 0;
}