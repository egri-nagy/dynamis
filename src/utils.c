#include "utils.h"

#include <stdio.h>

double random_double(){
  return rand() / (double)RAND_MAX;
}

void normalize(int size, double v[]){
  int  i;
  double sum = 0.0;
  for (i=0;  i  < size;  i++){
    sum  += v[i];
  }
  if (sum != 0.0){
    for (i=0;  i  < size;  i++){
      v[i]  /= sum;
    }
  }
}


void simplex(int size, double v[]){
  int  i;
  for (i=0;  i  < size;  i++){
    v[i]  = random_double();
  }
  normalize(size, v);
}

void generalized_simplex(int size, double sum, double v[]){
  int  i;
  simplex(size, v);
  for (i=0;  i  < size;  i++){
    v[i] *= sum;
  }
}
void multiply_vector(int size, double c, double *dst, double *src){
  int i;
  for (i = 0; i < size; i++){
    dst[i] = c * src[i];
  }  
}

void add_vectors(int size,double *c, double *a, double *b){
  int i;
  for (i = 0; i < size; i++){
    c[i] = a[i] + b[i];
  }
}

void copy_vector(int size, double dst[], double src[]){
  int i;
  for (i = 0; i < size; i++){
    dst[i] = src[i];
  }
}

void print_vector(int size, double y[], char *msg){
  int i;
  printf("%s ", msg);
  for (i = 0; i < size; i++){
    printf("%f ", y[i]);
  }
  printf("\n");
}
