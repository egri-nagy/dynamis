#ifndef ODE_SOLVER_H
#define  ODE_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "utils.h"

///the buffer size for reading inputfile
#define CHAR_BUFF_SIZE 256

///the maximum dimension of equations which the system can handle
#define MAX_DIMENSION 64

//######MEMBER FUNCTIONS FOR THE ODE DATA STRUCTURE -- OO style################################################################# 
///computes the values, parameters: number of equations, the t (time) value, the value of previous ys, the value of calculated ys
typedef void (*ode_functions)(int, double,  double[], double[], void *data);
///after a step there can be some modifications in the ode set
typedef void (*modify)(void*);
///after a step some information can be printed
typedef void (*dump)(FILE *out, double time, void *os);

///represents aan N dimensional set of ordinary differential equations
typedef struct{
  ///the dimenssion of the equation set i.e. the number of equations
  unsigned int N;
  ///the functions which calculate the values
  ode_functions functions;
  ///the functions for modification
  modify modification;
  ///information dumping
  dump data_dump;

  ///the current state and the next one
  double *state;
  double *nextstate;
  ///extra information for specific ODEs can be stored here
  void *data;
} ode_set;



//##############GLOBAL VARIABLES################################
///beginning of the interval to be integrated
double t0;
///end of the interval
double tn;
///number of steps within the interval
unsigned int number_of_steps;
///desired accuracy
double accuracy;



ode_set **load_ode_set(FILE *in);
void solve(ode_set **des);
void clean_up_ode_set(ode_set *des);
//shifts the second vector  to the first of the solution matrix (to be used if only  endvalues are  needed) 
void shift(ode_set *des);
ode_set *createOdeSet();
#endif  
//ODE_SOLVER_H
