#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"

#define MAXIT 50

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