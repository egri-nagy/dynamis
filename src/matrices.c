#include "matrices.h"
#include "utils.h"

#include <stdlib.h>

/**
 * Removes a row and a column from an m x n matrix. No memory is freed.  
 */
void shrinkMatrix(double **a, int m, int n, int row, int column){
  int i, j;
  //printf("hrink called %d %d, %d, %d, %d",a,m,n,row,column);
  fflush(stdout);
  //left bottom
  for (i = row+1;  i < m; i++){
    for (j = 0; j < column; j++){
      a[i-1][j]  = a[i][j];  
   }
  }
  //right top
  for (i = 0;  i < row; i++){
    for (j = column+1; j < n; j++){
      a[i][j-1]  = a[i][j];  
    }
  }
  //right bottom
  for (i = row+1;  i < m; i++){
    for (j = column+1; j < n; j++){
      a[i-1][j-1]  = a[i][j];  
   }
  }
 
}

void extendMatrix(double **a, int m, int n, double *row, double *column){
  int i;
  //filling the row
  for (i = 0; i < n+1 ; i++){
    a[m][i] = row[i];
  }
  //filling the column
  for (i = 0; i < m; i++){
    a[i][n] = column[i];
  }
}

void fillBlock(double **a, int ax, int ay, double **b, int bm, int bn){
  int i, j;
  for (i = 0; i < bm; i++){
    for (j = 0; j < bn; j++){
      a[i+ax][j+ay] = b[i][j];
    }    
  } 
}

void constant_multiplication(double **a, int m, int n, double c){
  int i,j; 
 for (i=0; i < m; i++){
    for (j=0; j < n; j++){
      a[i][j] *= c;
    }
  }
}

double **safeMatrixAlloc(int m, int n){
  int i;
  double **a = (double**) calloc(m, sizeof(double*));
  if (!a){ 
    fprintf(stderr, "Fail in safeMatrixAlloc()!\n");
  }
  for (i = 0; i < m; i++){
    a[i] =  (double*) calloc(n, sizeof(double));
    if (!a[i]){ 
      fprintf(stderr, "Fail in safeMatrixAlloc()!\n");
    }
  }
  return a;
}

void safeMatrixFree(double **a,int m, int n){
  int i;
  for (i = 0; i < m; i++){
    free(a[i]);
  }
  free(a);
}

double **hypercycleAlloc(int n){
  int i;
  double **a = safeMatrixAlloc(n,n);
  for (i=0; i < n; i++){
    a[i][(i+1)%n] = 1.0;
  }  
  return a;
}

double **randomHypercycleAlloc(int n){
  int i;
  double **a = safeMatrixAlloc(n,n);
  for (i=0; i < n; i++){
    a[i][(i+1)%n] = random_double();
  }  
  return a;
}

double **blockHypercycleAlloc(int n, int blocks){
  int i;
  double **tmp;
  double **a = safeMatrixAlloc(n,n);
  int blocksize = n / blocks;
  printf("%d\n", blocksize);
  if (n % blocks){ fprintf(stderr, "Blocksize not correct!");return 0;}
  for (i = 0; i < blocks; i++){
    tmp = uniformRandomAlloc(blocksize,1);
    fillBlock(a, blocksize*((i+1)%blocks),blocksize*i, tmp, blocksize, blocksize);
    safeMatrixFree(tmp, blocksize, blocksize);
  }
  return a;
}

double **uniformRandomAlloc(int n, double density){
  int i,j;
  double **a = safeMatrixAlloc(n,n);
  for (i=0; i < n; i++){
    for (j=0; j < n; j++){
      if (random_double() < density){
	a[i][j] = 1.0 - (2 * random_double());
      }
    }
  }  
  return a;
}


void printMatrix(FILE *out, double **a, int m, int n){
  int i,j;
  for (i = 0; i < m; i++){
    fprintf(out,"#");
    for (j = 0; j < n; j++){
  
    fprintf(out, "%f ", a[i][j]);
    }  
    fprintf(out, "\n");
  }
}

void readMatrix(FILE *in, double **a, int m, int n){
  int i;
  int j;
  float tmp;
  for (i = 0; i < m; i++){
    for (j = 0; j < n; j++){
      fscanf(in, "%f", &tmp);
      a[i][j] = tmp;
    }
  }
}
