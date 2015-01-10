#include "speciation.h"
#include "matrices.h"

#include <stdlib.h>
#include <time.h>

///interaction matrix
double **im;
///number of species - the dimension of the system - it changes all the time
int dim; 

//the result of mutate_vector
double mutated_row[MAX_DIMENSION], mutated_column[MAX_DIMENSION];

double extinction_threshold;
double speciation_threshold;

//######LOCAL FUNCTIONS##########################################################
void spec_load_parameters(ode_set *des,FILE *in);
double spec_score(int i, double y[]);
double spec_total(double y[]);
void spec_replicators(int N, double x, double y[], double result[], void *data);
void spec_modify(void *des);
void spec_dump(FILE *out,double time, void *des);
void mutate_vectors();
double chop(double x, double a, double b);
void spec_dump_dim(FILE *out, double time, void *o_s);



void createSpeciationODE(ode_set *des,FILE *in){
  dim = des->N;
  des->functions = &spec_replicators;
  des->data_dump = &spec_dump_dim;
  des->modification = &spec_modify;
  spec_load_parameters(des,in);
}

void spec_replicators(int N, double x, double y[], double result[], void *data){
  double tot = spec_total(y);
  int i;
  for (i = 0; i < N; i++){
    result[i] = y[i]*(spec_score(i,y) - tot);
  }
}

double spec_score(int  n, double y[]){
  int i;

  double sum = 0.0;
  //  printf("DEBUG: %d\n", n);
  for (i = 0; i < dim; i++){
    //printf("%f\n",im[n][i]);
    sum += (y[i]) * (im[n][i]);
  }   
  return  sum;
}

double spec_total(double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i] * spec_score(i,y);
  }
  return sum;
}

void spec_modify(void *o_s){
  int i,j;
  ode_set *os = (ode_set*) o_s;
  for (i = 0; i < dim; i++){
    if (os->state[i] < extinction_threshold){
      //printf("#shrinking\n");
      //printMatrix(stdout, im, dim,dim);
      shrinkMatrix(im, dim, dim, i,i);
      dim--;
      //printMatrix(stdout, im, dim,dim);
      for (j = i; j < dim; j++){
	os->state[j] = os->state[j+1];
      }
 
    }
  }// shrinking loop

  for (i = 0; i < dim; i++){
    if (os->state[i] > speciation_threshold){
      if( dim < MAX_DIMENSION-1){
	double tmp[MAX_DIMENSION];
	
	//printf("#Extending\n");
	for (j = 0; j < dim; j++){
	  tmp[j] = im[j][i];
	}
	mutate_vectors(im[i],tmp);
	extendMatrix(im, dim, dim, mutated_row, mutated_column);
	os->state[i] -= (2 * extinction_threshold);
	os->state[dim] = 2 * extinction_threshold;
	dim++;
	//printMatrix(stdout, im, dim,dim);
      }
      
    }// extension loop
  }
  os->N = dim;
  normalize(dim, os->state);
}

void spec_dump(FILE *out,double time, void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  fprintf(out, "%f ", time);
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  for (; i < MAX_DIMENSION; i++){
    fprintf(out,"0 "); 
  }
  fprintf(out,"\n");
}

void spec_dump_dim(FILE *out, double time, void *o_s){
  ode_set *os = (ode_set*) o_s;
  fprintf(out,"%d\n",os->N);
}

void spec_load_parameters(ode_set *des,FILE *in){
  int i,j;
  float tmp;
  int seed;
  im = (double**) calloc(MAX_DIMENSION, sizeof(double*));
  for (i = 0; i < MAX_DIMENSION; i++){
    im[i] = (double*) calloc(MAX_DIMENSION, sizeof(double));
  }
  fscanf(in, "%f " , &tmp);
  extinction_threshold = tmp;
  fscanf(in, "%f " , &tmp);
  speciation_threshold = tmp;
  fscanf(in, "%d", &seed);

  printf("#seed:%d\n",seed); 
  srand(seed);
 
  
  /*hypercycle
  im[0][4] = 1;
  im[1][0] = 1;
  im[2][1] = 1;  
  im[3][2] = 1;  
  im[4][3] = 1; 
  */
  //random matrix
  for(i = 0; i < 5; i++){
    for(j = 0; j < 5; j++){
      im[i][j] =  (1 - (2*(rand() / (double)RAND_MAX))); 
    }
  }
  
  printMatrix(stdout, im, dim,dim);
  fclose(in);
  fflush(stdout);
}

void mutate_vectors(double *row, double *column){
  int i;
  for (i = 0; i < dim; i++){
    mutated_row[i] = chop(row[i] + (.2*(rand() / (double)RAND_MAX)), -1.0,1.0);;
    mutated_column[i] =  chop(column[i] - (.5*(rand() / (double)RAND_MAX)), -1.0,1.0);
  }
  //  i = (int)(dim *(rand() / (double)RAND_MAX));

  mutated_row[dim] = 0;//0.1 * row[dim-1];//(1.0 - (2.0*(rand() / (double)RAND_MAX)));
}

double chop(double x, double a, double b){
  //sanity check
  if (a >= b){
    printf("Bad input forthe chop function: a >= b\n");
  }
  if (x > b){
    return b; 
  }
  else if(x < a){
    return a;  
  }
  else{
    return x; 
  }
}
