#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <likwid.h>
//#include <likwid-marker.h>
#include "utils.h"
#include "resolucoes_SL.h"

#define MAXIT 50

int main(){
    double *a, *b, *c, *d,  
    *resultadoEG, *resultadoGS, *resultadoEG3, *resultadoGS3,
    *dS, *aS, *cS,
    *bEG, *bGS, *bEG3, *bGS3;
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
    bEG = alocaVetor(n);
    bGS = alocaVetor(n);
    bEG3 = alocaVetor(n);
    bGS3 = alocaVetor(n);

    //alocando vetores resultado de cada metodo
    resultadoEG = alocaVetor(n);
    resultadoGS = alocaVetor(n);
    resultadoEG3 = alocaVetor(n);
    resultadoGS3 = alocaVetor(n);

    //lendo a matriz e o vetor b
    for(int i=0; i < n; ++i){
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
    copiaVetorResultado(b, bEG, n);
    copiaVetorResultado(b, bGS, n);
    copiaVetorResultado(b, bEG3, n);
    copiaVetorResultado(b, bGS3, n);
    
    // vetores matriz tridiagonais 
    separaTridiagonais(Matriz, a, d, c, n);
    separaTridiagonais(Matriz, aS, dS, cS, n);

    tolerancia = 0.0001;

    // ELIMINACAO DE GAUSS CLASSICA

    rtime_t tempoEG = timestamp();
    double *residuoEG = alocaVetor(n);
    LIKWID_MARKER_START ("eliminacaoDeGauss");
    eliminacaoDeGauss(MatrizEG, bEG, resultadoEG, n);
    LIKWID_MARKER_STOP ("eliminacaoDeGauss");
    tempoEG = timestamp() - tempoEG;
    residuoMatriz(MatrizEG, resultadoEG, bEG, residuoEG, n);

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
    gaussSeidel(MatrizGS, bGS, resultadoGS, n, tolerancia, &countGS);
    LIKWID_MARKER_STOP ("gaussSeidel");
    tempoGS = timestamp() - tempoGS;
    residuoMatriz(MatrizGS, resultadoGS, bGS, residuoGS, n);

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
    eliminacaoDeGaussTriDiagonais(a, d, c, bEG3, resultadoEG3, n);
    LIKWID_MARKER_STOP ("eliminacaoDeGaussTriDiagonais");
    tempoEG3 = timestamp() - tempoEG3;
    residuoEliminacaoDeGaussTriDiagonais(a, d, c, bEG3, resultadoEG3, residuoEG3, n);

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
    gaussSeidelTriDiagonais(aS, dS, cS, bGS3, resultadoGS3, n, tolerancia, &countGS3);
    LIKWID_MARKER_STOP ("gaussSeidelTriDiagonais");
    tempoGS3 = timestamp() - tempoGS3;
    residuoEliminacaoDeGaussTriDiagonais(aS, dS, cS, bGS3, resultadoGS3, residuoGS3, n);

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
    free(bEG);
    free(bGS);
    free(bEG3);

    LIKWID_MARKER_CLOSE;
    return 0;
}
