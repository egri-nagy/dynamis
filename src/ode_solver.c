#include "ode_solver.h"
#include "replicator.h"
#include "speciation.h"
#include "spatial_replicator.h"
#include "ecology.h"

//#####  LOCAL FUNCTION DEFINITIONS ######################################################
void runge_kutta_4th(double  t,  ode_set  *des, double h);
int isFixedPoint(ode_set *des);
//#####GLOBAL FUNCTIONS IMPLEMENTATION####################################################
ode_set **load_ode_set(FILE *in){
  char s[CHAR_BUFF_SIZE];
  char type[CHAR_BUFF_SIZE];
  int i;
  float  f;
  ode_set **os_array = (ode_set**) calloc(2, sizeof(ode_set**));
  ode_set *des = (ode_set*) calloc(1, sizeof(ode_set));
  ///loading the general values
  //type  
  fscanf(in,"%*s %s ", type);
  printf("#TYPE: %s\n", type);
  //type  
  fscanf(in,"%*s %s ", s);
  printf("#METHOD: %s\n", s);
  //dimension
  fscanf(in ,"%*s %d", &(des->N));
  printf("#DIMENSION: %d\n", des->N);
  //interval
  fscanf(in," %*s %f ", &f);
  t0 = f;
  fscanf(in, "%f ", &f);
  tn = f;
  printf("#INTERVAL STARTS: %f\n", t0);
  printf("#INTERVAL ENDS: %f\n", tn);
  //steps
  fscanf(in ,"%*s %d", &number_of_steps);
  printf("#STEPS: %d\n", number_of_steps);

  //memory allocation for the state variables
  des->state =  (double*) calloc(MAX_DIMENSION, sizeof(double));
  des->nextstate =  (double*) calloc(MAX_DIMENSION, sizeof(double));
  
  fscanf(in," %s ",s);
  if (!strcmp("INITIAL",s)){
    printf("#INITIAL VALUES: ");
    for  (i  = 0; i < des->N; i++){
      fscanf(in," %f ", &f);
    printf(" %f", f);
    
    des->state[i] = f;
    }
    printf("\n");
    fscanf(in, " %*s ");
  }

  ///loading the specific values from the rest of the file
  if (!strcmp(type,"replicator")){
    createReplicatorODE(des,in);
  }
  else if (!strcmp(type,"speciation")){
    createSpeciationODE(des, in);
  }
  else if (!strcmp(type,"ecology")){
    createEcologyODE(des, in);
  }
  else if (!strcmp(type,"spatial_replicator")){
    free(os_array);
    os_array = createSpatialReplicatorODE(des,in);
    return os_array;
  }
  os_array[0] =des;
  os_array[1] = NULL;
  return  os_array;
}

///solve the array of the ode sets, the array is closed with a nullpointer (it's assumed that they have the same interval)
void solve(ode_set **des){
  int i; 
  int non_coupled = 0;
  double h = (tn  -  t0)  /  (double) number_of_steps;
  ode_set **osp;
  non_coupled =!((int) des[1]);
  //the main loop
  for (i = 0; i < (number_of_steps - 1); i++){
    osp = des;
    while (*osp){
      runge_kutta_4th(t0 + (i * h), *osp, h);
      if(non_coupled){ //if it's not a coupled system
	if (isFixedPoint(*osp)){
	  printf("#Fixed-point reached!\n");
	  if ((*osp)->data_dump){(*osp)->data_dump(stdout,(t0+(i*h)), *osp);}   
	  return;
	}
      }
      shift(*osp);
      if ((*osp)->data_dump){(*osp)->data_dump(stdout,(t0+(i*h)), *osp);}   
      if ((*osp)->modification){(*osp)->modification(*osp);}
      osp++;
    }
  }
}

ode_set *createOdeSet(){
  return (ode_set*) calloc(1, sizeof(ode_set));
}



///midvectors for rungeKutta4th
double k1[MAX_DIMENSION];
double k2[MAX_DIMENSION];
double k3[MAX_DIMENSION];
double k4[MAX_DIMENSION];
double yhk1[MAX_DIMENSION];
double yhk2[MAX_DIMENSION];
double yhk3[MAX_DIMENSION];
double yhk4[MAX_DIMENSION];
double h1[MAX_DIMENSION];
double h2[MAX_DIMENSION];
double h3[MAX_DIMENSION];
double h4[MAX_DIMENSION];

//a step function: it calculates the next state with fourth order Runge-Kutta method
void runge_kutta_4th(double t, ode_set *des, double h){
  int i;

  des->functions(des->N,t, des->state, k1 , des->data);
  multiply_vector(des->N, h / (2.0), h1 , k1 );  
  add_vectors(des->N, yhk1 , des->state , h1 );
  
  des->functions(des->N, t + (h / 2.0), yhk1 , k2 ,des->data);
  multiply_vector(des->N, h / (2.0), h2 , k2 );  
  add_vectors(des->N, yhk2 , des->state , h2 );
  
  des->functions(des->N, t + (h / 2.0), yhk2 , k3 ,des->data);
  multiply_vector(des->N, h, h3 , k3 );  
  add_vectors(des->N, yhk3 , des->state , h3 );
  
  des->functions(des->N, t + h, yhk3 , k4 ,des->data); 
  
  for (i = 0; i  < des->N; i++){
    des->nextstate[i]  = des->state[i] + (h/6.0)*(k1 [i] + (2*k2 [i]) + (2*k3 [i]) + k4 [i]);
  }
}


void shift(ode_set *des){
  register int i;
  for (i = 0;  i < des->N;  i++){
    des->state[i] = des->nextstate[i];
  }
}

//##########IMPLEMENTATION OF LOCAL FUNCTIONS#############################

int isFixedPoint(ode_set *des){
  register int i;
  for (i = 0;  i < des->N;  i++){
    if (des->state[i] != des->nextstate[i]) return 0;
  }
  return 1;
}
