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

double residuoMatriz(double **A, double *x, double *b, double *residuo, int n) {
    double *Ax = (double *)malloc(n * sizeof(double));

    //Calcula Ax
    for (int i = 0; i < n; ++i) {
        Ax[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    //Calcula o vetor de residuo
    for (int i = 0; i < n; ++i) {
        residuo[i] = Ax[i] - b[i];
        printf("valor Ax[%d] = %.2f/n Valor b[%d] = %.2f\n", i, Ax[i], i, b[i]);
    }

     //Libera memoria alocada
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
    //retrosubstituicao
    for(int i = n-1; i >= 0; --i){
        x[i] = b [i];
        for(int j = i+1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    }
}

// void gaussSeidel(double A**, double *b, double *x, uint n, double tol){
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
    //Triangularizacao
    for(int i = 0; i < n-1; ++i){
            int k = i+1;
            double m = a[i]/d[i];
            a[i] = 0;
            d[i+1] -= c[i] * m;
            b[i+1] -= b[i] * m;     
    }
    //Retrosubstituicao
    x[n-1] = b[n-1] / d[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (b[i] - c[i] * x[i+1]) / d [i];
    }
}

int main(){
    double **Matriz;
    double *a, *b, *c, *d, *x;
    double *residuo;
    int n;

    //leitura dimensoes da mateiz
    scanf("%d", &n);

   // Alocando memoria para a matriz
    Matriz = (double **)malloc(n * sizeof(double *));
    
    // Alocando memoria para o vetor b
    b = (double *)malloc(n * sizeof(double));
    // Alocando memoria vetor residuo
    residuo = (double *)malloc(n * sizeof(double));
    // Alocando memória para os vetores da matriz tridiagonal
    d = (double *)malloc(n * sizeof(double));
    a = (double *)malloc(n * sizeof(double)); 
    c = (double *)malloc(n*  sizeof(double));
    x = (double *)malloc(n*  sizeof(double));

    // Lendo a matriz e o vetor b
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
    for(int i = 0; i < n; ++i){
            d[i] = Matriz[i][i];
            if(i < n-1){
                a[i] = Matriz[i][i+1];
                c[i] = Matriz[i+1][i];
            }   
    }

    eliminacaoDeGauss(Matriz, b, x, n);
    residuoMatriz(Matriz, x, b, residuo, n);
    //tempo = timestemp()
    //eliminacaoDeGaussTriDiagonais(a, d, c, b, x, n);
    //timestemp() - tempo

    printf("Solução do sistema:\n");
    for (int i = 0; i < n; ++i) {
        printf("x[%d] = %.2f\n", i, x[i]);
    }

    printf("Residuo do sistema:\n");
    for (int i = 0; i < n; ++i) {
        printf("residuo[%d] = %.2f\n", i, residuo[i]);
    }

    // mostra a matriz lida fazer funcao
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
    //printando os 3 vetores 
    for(int i=0; i<n; ++i)
        printf("d[%d] = %.2f\n", i, d[i]);
    printf("\n");
    for(int i=0; i<n-1; ++i)
        printf("a[%d] = %.2f\n", i, a[i]);
    printf("\n");
    for(int i=0; i<n-1; ++i)
        printf("c[%d] = %.2f\n", i, c[i]);
    printf("\n");

    //liberando memoria
    for(int i = 0; i < n; ++i){
        free(Matriz[i]);
    }
    free(Matriz);
    free(a);
    free(b);
    free(d);
    free(c);
    free(x);
    free(residuo);

    return 0;
}
/*
EXEMPLO 

3
1 -3 2 11
-2 8 -1 -15
4 -6 5 29

4
4 -1 0 0 3
-1 4 -1 0 6
0 -1 4 -1 1
0 0 -1 4 12

4
2 1 0 0 5
3 4 2 0 6
0 2 5 1 7
0 0 3 6 8
*/