#include "replicator.h"
#include "matrices.h"

#include <stdlib.h>
#include <time.h>

///interaction matrix
double **im;

///number of species (strategies, molecules...)  - the dimension of the system
int dim; 

//##########################LOCAL FUNCTION DEFINITIONS#############################################
void replicator_initialize(ode_set *des,FILE *in);
double replicator_score(int i, double y[]);
double replicator_total(double y[]);
void replicator_dump(FILE *out,double time, void *o_s);

//the constructor function
void createReplicatorODE(ode_set *des,FILE *in){
  dim = des->N;
  des->functions = &replicators;
  des->data_dump = &replicator_dump;
  replicator_initialize(des,in);
}


//##################################LOCAL FUNCTION IMPLEMENTATIONS#################################
void replicators(int N, double x, double y[], double result[],void *data){
  double tot = replicator_total(y);
  int i;
  for (i = 0; i < N; i++){
    result[i] = y[i]*(replicator_score(i,y) - tot);
  }
}

//the score (fitness, viability, payoff) of the nth species
double replicator_score(int  n, double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i]*im[n][i];
  }   
  return  sum;
}

//the total score of the system
double replicator_total(double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i] * replicator_score(i,y);
  }
  return sum;
}

//dumps the state of the system
void replicator_dump(FILE *out,double time, void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  fprintf(out, "%f ", time);
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  fprintf(out,"\n");
}


void replicator_initialize(ode_set *des,FILE *in){
  time_t t;
  char buff[CHAR_BUFF_SIZE];


  while( fscanf(in, " %s ", buff) != EOF){
    if (!strcmp(buff, "INTERACTION_MATRIX")){
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
      if (!strcmp(buff,"full")){
	des->data_dump = &replicator_dump;
      }
    }
  }

  time(&t);
  srand(t);

  fclose(in);
  fflush(stdout);
}
