#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "resolucoes_SL.h"

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
    //triangularizacao
    for(int i=0; i<n; ++i){
        //int iPivo = encontraMax(A, i, n);
        //if(i != iPivo)
            //trocaLinhas(A, b, i, iPivo);

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

void gaussSeidel(double **A, double *b, double *x, int n, double tol, int *count){
    double erro = 1.0;
    double s;
    int j;
    double *x_antigo = malloc(n * sizeof(double));

    *count = 0;

    while(erro > tol && *count < MAXIT ){
        for(int i = 0; i < n; i++){
            x_antigo[i] = x[i];
            s = 0.0;
            for(j=0; j < n; ++j){
                if(i != j)
                    s = (A[i][j] * x[j]) + s;
            }
            x[i] = (b[i] - s)/A[i][i];
        }

        if(*count == 0)
            erro = 1.0;
        else
            erro = encontrarMaiorSubtracao(x, x_antigo, n);
        (*count)++;
    }
    free(x_antigo);
}

void gaussSeidelTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n, double tol, int *countGS3){
    double erro = 1.0;
    int j, s;
    double *x_antigo = malloc(n * sizeof(double));

    *countGS3 = 0;

    while(erro > tol && *countGS3 < MAXIT ){
            s = 0.0;

            x[0] = (b[0] - (c[0] * x[1])) / d[0];

            for(int i = 0; i < n; ++i){
                x_antigo[i] = x[i];
            }
            for (int i=1; i < n-1; ++i){
                x[i] = (b[i] - a[i-1] * x[i-1] - c[i] * x[i+1]) / d[i];
            }

            x[n-1] = (b[n-1] - a[n-2] * x[n-2]) / d[n-1];

        if(*countGS3 == 0)
            erro = 1.0;
        else
            erro = encontrarMaiorSubtracao(x, x_antigo, n);
        (*countGS3)++;
    }
}

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

int encontraMax(double **A, int i, int n) {
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

double encontrarMaiorSubtracao(double *a, double *b, int n){
    double maior = 0.0;
    for(int i = 0; i < n; ++i){
        if(fabs(a[i] - b[i]) > maior){
            maior = fabs(a[i] - b[i]);
        }
    }
    //printa vetores
    return maior;
}

double *alocaVetor(int qntPontos)
{
    double *vetor;
    vetor = (double *)calloc(qntPontos, sizeof(double));
    return vetor;
}

void trocaLinhas(double **A, double *b, int i, int iPivo){
    double *aux = A[i];
    A[i] = A[iPivo];
    A[iPivo] = aux;

    double auxB = b[i];
    b[i] = b[iPivo];
    b[iPivo] = auxB;
}

void copiaVetorResultado(double *b, double *copia, int qntPontos)
{
    for (int i = 0; i < qntPontos; i++)
    {
        copia[i] = b[i];
    }
}

void copiaMatriz(double **A, double **B, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            B[i][j] = A[i][j];
    }
}
