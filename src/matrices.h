#ifndef MATRICES_H
#define MATRICES_H

#include <stdio.h>

/**
 * Functions for matrix handling and creating special matrices. Here matrices
 * are represented as two-dimensional double arrays. 
 */

//###########MANIPULATING MATRICES#############################################
/**
 * Removes a row and a column from an m x n matrix. No memory is freed.  
 */
void shrinkMatrix(double **a, int m, int n, int row, int column);

/**
 * Extends a matrix with a row and a column. The m,n are the original sizes. 
 */
void extendMatrix(double **a, int m, int n, double *row, double *column);

/**
 * Fills a block of a matrix with a specified submatrix.
 */ 
void fillBlock(double **a, int ax, int ay, double **b, int bx, int by);

void constant_multiplication(double **a, int m, int n, double c);

//#########MATRIX CREATION#####################################################
/**
 * Allocates a matrix with a given size.
 */
double **safeMatrixAlloc(int m, int n);
void safeMatrixFree(double **a, int m, int n);

double **hypercycleAlloc(int n);
double **randomHypercycleAlloc(int n);
double **blockHypercycleAlloc(int n, int blocks);
double **uniformRandomAlloc(int n, double density);

/**
 * Prints a matrix to specified file. 
 */
void printMatrix(FILE *out, double **a, int m, int n);

/**
 * Reads a matrix from a specified file. 
 */
void readMatrix(FILE *in, double **a, int m, int n);

#endif //MATRICES_H
