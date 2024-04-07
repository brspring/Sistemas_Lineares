#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <likwid.h>

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

int main(){
    double *a, *b, *c, *d,  
    *resultadoEG, *resultadoGS, *resultadoEG3, *resultadoGS3,
    *dS, *aS, *cS;
    double *residuo;
    int n;
    int countGS = 0;
    int countGS3 = 0;
    double tolerancia;

    LIKWID_MARKER_INIT;

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
    resultadoEG3 = alocaVetor(n);
    resultadoGS3 = alocaVetor(n);

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

    // ELIMINACAO DE GAUSS CLASSICA

    rtime_t tempoEG = timestamp();
    double *residuoEG = alocaVetor(n);
    LIKWID_MARKER_START ("eliminacaoDeGauss");
    eliminacaoDeGauss(MatrizEG, b, resultadoEG, n);
    LIKWID_MARKER_STOP ("eliminacaoDeGauss");
    tempoEG = timestamp() - tempoEG;
    residuoMatriz(MatrizEG, resultadoEG, b, residuoEG, n);

    printf("EG clássico:\n");
    printf("%.8f ms\n", tempoEG);
    for (int i = 0; i < n; ++i)
        printf("%.12f   ", resultadoEG[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuoEG[i]);
    printf("\n\n");

    // GAUSS SEIDEL CLASSICO
    rtime_t  tempoGS = timestamp();
    double *residuoGS = alocaVetor(n);
    LIKWID_MARKER_START ("gaussSeidel");
    gaussSeidel(MatrizGS, b, resultadoGS, n, tolerancia, &countGS);
    LIKWID_MARKER_STOP ("gaussSeidel");
    tempoGS = timestamp() - tempoGS;
    residuoMatriz(MatrizGS, resultadoGS, b, residuoGS, n);

    printf("GS Clássico [%d]:\n", countGS);
    printf("%.8f ms\n", tempoGS);
    for (int i = 0; i < n; ++i) {
        printf("%.12f  ", resultadoGS[i]);
    }
    printf("\n");
    for (int i = 0; i < n; ++i) {
        printf("%.12f  ", residuoGS[i]);
    }
    printf("\n\n");

    // GAUSS SEIDEL 3 DIAGONAIS
    rtime_t  tempoEG3 = timestamp();
    double *residuoEG3 = alocaVetor(n);
    LIKWID_MARKER_START ("eliminacaoDeGaussTriDiagonais");
    eliminacaoDeGaussTriDiagonais(a, d, c, b, resultadoEG3, n);
    LIKWID_MARKER_STOP ("eliminacaoDeGaussTriDiagonais");
    tempoEG3 = timestamp() - tempoEG3;
    residuoEliminacaoDeGaussTriDiagonais(a, d, c, b, resultadoEG3, residuoEG3, n);

    printf("EG 3-diagonal:\n");
    printf("%.8f ms\n", tempoEG3);
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", resultadoEG3[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuoEG3[i]);
    printf("\n\n");

    // GAUSS SEIDEL 3 DIAGONAIS
    rtime_t  tempoGS3 = timestamp();
    double *residuoGS3 = malloc(n * sizeof(double));
    LIKWID_MARKER_START ("gaussSeidelTriDiagonais");
    gaussSeidelTriDiagonais(aS, dS, cS, b, resultadoGS3, n, tolerancia, &countGS3);
    LIKWID_MARKER_STOP ("gaussSeidelTriDiagonais");
    tempoGS3 = timestamp() - tempoGS3;
    residuoEliminacaoDeGaussTriDiagonais(aS, dS, cS, b, resultadoGS3, residuoGS3, n);

    printf("GS 3-diagonal[%d]:\n", countGS3);
    printf("%.8f ms\n", tempoGS3);
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", resultadoGS3[i]);
    printf("\n");
    for (int i = 0; i < n; ++i)
        printf("%.12f  ", residuoGS3[i]);
    printf("\n\n");

    desalocaMatriz(Matriz, n);
    desalocaMatriz(MatrizGS, n);
    desalocaMatriz(MatrizEG, n);
    free(a);
    free(residuoGS);
    free(residuoEG);
    free(residuoEG3);
    free(residuoGS3);
    free(resultadoEG3);
    free(resultadoGS3);
    free(resultadoEG);
    free(resultadoGS);
    free(dS);
    free(aS);
    free(cS);
    free(b);
    free(d);
    free(c);
    free(residuo);

    LIKWID_MARKER_CLOSE;
    return 0;
}