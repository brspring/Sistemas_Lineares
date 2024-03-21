#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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

void trocaLinhas(double **A, double *b, int i, int iPivo){
    double *aux = A[i];
    A[i] = A[iPivo];
    A[iPivo] = aux;

    double auxB = b[i];
    b[i] = b[iPivo];
    b[iPivo] = auxB;
}

void eliminacaoDeGauss(double **A, double *b, uint n){
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
            
    }
}
int main(){
    double **Matriz;
    double *b;
    uint n;
    //leitura de Matriz e b
    scanf("%d", &n);

   // Alocando memória para a matriz
    Matriz = (double **)malloc(n * sizeof(double *));
    
    // Alocando memória para o vetor b
    b = (double *)malloc(n * sizeof(double));

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

    eliminacaoDeGauss(Matriz, b, n);

    // mostra a matriz lida
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
    
    for(int i = 0; i < n; ++i){
        free(Matriz[i]);
    }
    free(Matriz);
    free(b);
}

/*
EXEMPLO 

3
1 -3 2 11
-2 8 -1 -15
4 -6 5 29
*/