#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// rtime_t: tipo usado para representar valores de tempo em ponto flutuante
typedef double rtime_t;

// Funções
rtime_t timestamp(void);

double **alocaMatriz(int ordem);

void desalocaMatriz(double **matriz, int ordem);

void separaTridiagonais(double **Matriz, double *a, double *d, double *c, int n);

int encontraMax(double **A, int i, int n);

double encontrarMaiorSubtracao(double *a, double *b, int n);

double *alocaVetor(int qntPontos);

void trocaLinhas(double **A, double *b, int i, int iPivo);

void copiaVetorResultado(double *b, double *copia, int qntPontos);

void copiaMatriz(double **A, double **B, int n);

#endif // __UTILS_H__