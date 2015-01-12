#include  "ecology.h"

#include "matrices.h"
#include <time.h>
#include <stdlib.h>

///interaction matrix
double **im;
///number of species - the dimension of the system
int dim; 

///common data
double carrying_capacity;

void ecology_initialization(ode_set *des,FILE *in);
double ecology_score(int i, double y[]);
double ecology_total(double y[]);
void ecology(int N, double x, double y[], double result[],void *data);
void ecology_full_dump(FILE *out,double time, void *o_s);
void ecology_dump(FILE *out,double time, void *o_s);

void createEcologyODE(ode_set *des,FILE *in){
  printf("#Creating Ecology ODE \n");
  dim = des->N;
  des->functions = &ecology;
  ecology_initialization(des,in);
}

void ecology(int N, double x, double y[], double result[],void *data){
  int i;
  double growth_rate =pow( ((carrying_capacity - ecology_total(y)) / carrying_capacity), 3);

  for (i = 0; i < N; i++){
    double fitness = ecology_score(i,y);
    if (fitness > 0.0){
      result[i] = growth_rate * y[i] * (fitness);
    }
    else{
      result[i] = y[i] * fitness;   
 }
  }
}

double ecology_score(int  n, double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i]*im[n][i];
  }   
  return  sum;
}

double ecology_total(double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i];
  }
  return sum;
}

//time + total + state
void ecology_full_dump(FILE *out,double time, void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  fprintf(out, "%f %f ", time, ecology_total(os->state));
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  fprintf(out,"\n");
}

//dumps without time and total: only state
void ecology_dump(FILE *out,double time, void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  fprintf(out,"\n");
}

void ecology_initialization(ode_set *des,FILE *in){
  float tmp;
  time_t t;
  char buff[CHAR_BUFF_SIZE];


  while( fscanf(in, " %s ", buff) != EOF){
    if (!strcmp(buff, "CARRYING_CAPACITY")){
      //reading carrying capacity
      fscanf(in,"%f ",&tmp);
      carrying_capacity = tmp;
    }
    else if (!strcmp(buff, "INTERACTION_MATRIX")){
      //identifying the type of the interaction matrix
      fscanf(in, "%s ", buff);
      if (!strcmp(buff,"custom")){
	im = safeMatrixAlloc(dim,dim);
	readMatrix(in, im, dim,dim);
	time(&t);
	printf("#Init random with clock seed: %d\n", (int) t);
	srand(t);
      }
      else if (!strcmp(buff,"uniform_random")){
	int seed;
	fscanf(in, "%d", &seed);
	srand(seed);
	im = uniformRandomAlloc(dim,.7);
      }
    }
    else if (!strcmp(buff, "DUMP")){
      //identifying the dumping method
      fscanf(in, "%s ", buff);
      if (!strcmp(buff,"state")){
	des->data_dump = &ecology_dump;
      }
      else if (!strcmp(buff,"full")){
	des->data_dump = &ecology_full_dump;
      }
    }
  }
  

  
  time(&t);
  printf("#Init random with clock seed: %d\n", (int) t);
  srand(t);
  
  generalized_simplex(dim,carrying_capacity*0.4,des->state);
  //displaying the matrix 
  printMatrix(stdout, im, dim, dim); 

  fclose(in);
  fflush(stdout);
}
