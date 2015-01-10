#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <math.h>

///returns  a random number [0,1]
double random_double(); 
///normalizes a vector
void normalize(int size, double v[]);
///creates  a  simplex into the preallocated vector 
void simplex(int size, double v[]);
///creates  a  simplex into the preallocated vector 
void generalized_simplex(int size, double sum, double v[]);
///c = a + b
void add_vectors(int size, double c[], double a[], double b[]);
///dst = c * src
void multiply_vector(int size, double c, double dst[], double src[]);

void copy_vector(int size, double dst[], double src[]);
void print_vector(int size, double y[], char *msg);
#endif
