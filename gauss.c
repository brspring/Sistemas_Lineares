#include <stdio.h>
#include <stdlib.h>


int encontraMax(double **A, int i, int n){
    int max = A[i][0];
    for(int j=1; j<n; ++j){
        if(A[i][j] > max)
            max = j;
    }
    return max;
}

void trocaLinhas(double **A, double *b, int i, int j){
    double *aux = A[i];
    A[i] = A[j];
    A[j] = aux;

    double auxB = b[i];
    b[i] = b[j];
    b[j] = auxB;
}

void eliminacaoDeGauss(double **A, double *b, int n){
    for(int i=0; i<n; ++i){
        int iPivo = encontraMax(A, i, n);
        if(i != iPivo)
            trocaLinhas(A, b, i, iPivo);

        for(int k=i+1; k < n; ++k){
            double m = A[k][i]/A[i][i];
            A[k][i] = 0;
            for(int j=i; j<n; ++j){
                A[k][j] -= m*A[i][j];
            }
            b[k] -= m * b[i];
        }        
            
    }
}
int main(){
    double **Matriz;
    double *b;
    int n;
    //leitura de Matriz e b
    scanf("%d", &n);
    //Matriz = (double **)malloc(n*sizeof(double *));

    for(int i=0; i<n; ++i){
        //Matriz[i] = (double *)malloc((n+1)*sizeof(double));
        for(int j=0; j<=n; ++j){
            if(j==n){
                scanf("%lf", &b[i]);
            }else{
                scanf("%lf", &Matriz[i][j]);
            }
        }
    }
    // mostra a matriz lida
    for(int i=0; i<=n; ++i){
        for(int j=0; j<=n; ++j){
            printf("%lf ", Matriz[i][j]);
        }
        printf("\n");
    }
    //eliminacaoDeGauss(Matriz, b, n);
}