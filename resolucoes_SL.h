#ifndef __RESOLUCOES_SL_H__
#define __RESOLUCOES_SL_H__

void residuoEliminacaoDeGaussTriDiagonais(double *a, double *d, double *c, double *b, double *x, double *residuo, int n);

double residuoMatriz(double **A, double *x, double *b, double *residuo, int n);

void eliminacaoDeGauss(double **A, double *b, double *x, int n);

void gaussSeidel(double **A, double *b, double *x, int n, double tol, int *count);

void gaussSeidelTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n, double tol, int *countGS3);

void eliminacaoDeGaussTriDiagonais(double *a, double *d, double *c, double *b, double *x, int n);

double **alocaMatriz(int ordem);

void desalocaMatriz(double **matriz, int ordem);

void separaTridiagonais(double **Matriz, double *a, double *d, double *c, int n);

int encontraMax(double **A, int i, int n);

double encontrarMaiorSubtracao(double *a, double *b, int n);

double *alocaVetor(int qntPontos);

void trocaLinhas(double **A, double *b, int i, int iPivo);

void copiaVetorResultado(double *b, double *copia, int qntPontos);

void copiaMatriz(double **A, double **B, int n);
#endif // __RESOLUCOES_SL_H__
