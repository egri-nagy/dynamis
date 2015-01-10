#include  "spatial_ecology.h"

#include "matrices.h"
#include "png_util.h"
#include <time.h>
#include <stdlib.h>

//TODO!!!!! make it common
typedef struct{
  int x;
  int y;
} location;


///interaction matrix
double **im;
///number of species - the dimension of the system
int dim; 

double diffusion_rate;
int xsize;
int ysize;
ode_set ***lattice;

//for png writing
int spe_serial = 0;
/*png_byte palette[39] = {195,0,0,
		     192,122,41,
		     138,239,229,
		     42,75,239,
		     235,239,7,
		     235,239,255,
		     0,0,7,
		     0,239,0,
		     235,0,7,
		     100,100,100,
		     20,20,130,
		     235,135,35,
		     33,168,23};

png_byte palette[30] = {123,144,118,
			120,146,132,
			101,146,121,
			65,146,101,
			3,146,68,
			162,143,140,
			162,120,112,
			162,93,81,
			162,75,60,
			162,50,30};

*/
png_byte spe_palette[30] = {253,67,25,
			252,98,21,
			252,155,10,
			252,204,10,
			251,243,10,
			0,202,251,
			23,168,252,
			23,132,250,
			23,84,252,
			23,27,252};

png_byte *spe_image = 0;
int width;
int height;
char spe_picname[CHAR_BUFF_SIZE];


///common data
double spe_carrying_capacity  = 2;

void spe_dump_pics(FILE *out, double time,void *o_s);
void spe_load_parameters(ode_set *des,FILE *in);
double spatial_ecology_score(int i, double y[]);
double spatial_ecology_total(double y[]);
void spatial_ecology(int N, double x, double y[], double result[],void *data);
void spatial_ecology_dump(FILE *out,double time, void *o_s);
void diffuse(void *o_s);
location *spe_createLocation(int x, int y);

ode_set **createSpatialEcologyODE(ode_set *des,FILE *in){
  ode_set **os_array;
  int i;
  int  j;
  int  k;
  
  dim = des->N;
  spe_load_parameters(des,in);
  os_array = (ode_set**) calloc((xsize-1) * (ysize-1), sizeof(ode_set*));
  k = 0; 
  for (i = 1; i < (xsize - 1);i++){
    for (j = 1; j < (ysize - 1);j++){
      os_array[k] = lattice[i][j];
      k++;
    }
  }
  return os_array;
}

void spatial_ecology(int N, double x, double y[], double result[],void *data){
  int i;
  for (i = 0; i < N; i++){
    result[i] = ((spe_carrying_capacity-y[i])/spe_carrying_capacity)*y[i]*spatial_ecology_score(i,y);
  }
}

double spatial_ecology_score(int  n, double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i]*im[n][i];
  }   
  return  sum;
}

double spatial_ecology_total(double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i];
  }
  return sum;
}

void diffuse(void *o_s){
  //diffusion
  int i,j,k;
  ode_set *os = (ode_set*) o_s;
  location *loc = (location*) os->data;  
  //the neighbour values are gathered in the next vector 
  for ( k = 0; k < dim; k++){
    os->nextstate[k] += ((rand() / (double) RAND_MAX) * diffusion_rate * (lattice[(loc->x)-1][loc->y]->state[k]
					     + lattice[(loc->x)+1][loc->y]->state[k]
					     + lattice[(loc->x)][(loc->y)-1]->state[k]
					     + lattice[(loc->x)][(loc->y)+1]->state[k]));
   os->nextstate[k] += ((rand() / (double) RAND_MAX) * diffusion_rate * (lattice[(loc->x)-1][(loc->y)-1]->state[k]
					     + lattice[(loc->x)-1][(loc->y)+1]->state[k]
					     + lattice[(loc->x)+1][(loc->y)-1]->state[k]
					     + lattice[(loc->x)+1][(loc->y)+1]->state[k]));
  }
  //the last cell fires the whole update process 
 if (loc->x == xsize-2 && loc->y == ysize-2){
    for (i = 1; i < xsize-1; i++){
      for (j = 1; j < ysize-1; j++){
	add_vectors(dim, lattice[i][j]->state, lattice[i][j]->state, lattice[i][j]->nextstate);
	normalize(dim, lattice[i][j]->state);
      }
    }
  }
}

void spatial_ecology_dump(FILE *out,double time, void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  fprintf(out, "%f %f ", time, spatial_ecology_total(os->state));
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  for (; i < dim; i++){
    fprintf(out,"0 "); 
  }
  fprintf(out,"\n");
}

void spe_dump_pics(FILE *out, double time,void *o_s){
  int i,j,k,imax;
  double max;
  ode_set *os = (ode_set*) o_s;
  ode_set *tmp;
  location *loc = (location*) os->data; 
  if (!spe_image){
     spe_image = (png_byte*) calloc((width * height), sizeof(png_byte));
  }

  if (loc->x == (xsize-2) && loc->y == (ysize-2)){

    for (i = 1; i < xsize-1; i++){
      for (j = 1; j < ysize-1; j++){

	tmp = lattice[i][j];
	max = -1;
	imax = 0;
	for (k=0; k < dim; k++){

	  if ((tmp->state[k]) > max){
	    max = tmp->state[k];
	    imax = k;
	  }//if
	}//k	  

	spe_image[((j-1)*width)+(i-1)] = imax;
      }
    }
    sprintf(spe_picname,"dynamis%6d.png", spe_serial++);
    for (i = 0; i < strlen(spe_picname);i++){
      if (spe_picname[i] == ' ') spe_picname[i] = '0';
    }
    write_indexed_png(spe_picname, spe_image, width, height, spe_palette, dim);
  }//if
}



void spe_load_parameters(ode_set *des,FILE *in){
  int i,j;
  float tmp;
  int seed;
  char puffer[CHAR_BUFF_SIZE];

  fscanf(in,"%*s %d ", &xsize);
  fscanf(in,"%*s %d ", &ysize);
  width = xsize -2;
  height = ysize -2;
  fscanf(in,"%*s %f ", &tmp);
  diffusion_rate = tmp;
  fscanf(in,"%*s %d ", &seed);
  srand(seed);
  fscanf(in,"%*s %s ", &(puffer[0]));
  
  if (!strcmp(puffer, "hypercycle")){
    printf("#Hypercycle!");
    im = hypercycleAlloc(dim); 
  }
  else if (!strcmp(puffer, "random_hypercycle")){
    im = randomHypercycleAlloc(dim); 
  }
  else if (!strcmp(puffer, "block_hypercycle")){
    im = blockHypercycleAlloc(dim, 5); 
  }
  else if (!strcmp(puffer, "hypercycle_blocks")){
    double **tm;
    fscanf(in,"%f ", &tmp);
    im = safeMatrixAlloc(dim,dim);//uniformRandomAlloc(dim,tmp);
    constant_multiplication(im, dim, dim, 0.1);
    tm = hypercycleAlloc(5);
    fillBlock(im, 0,0, tm, 5,5);
    safeMatrixFree(tm,5,5);
    tm = hypercycleAlloc(5);
    fillBlock(im, 5,5, tm, 5,5);
    safeMatrixFree(tm,5,5);    
    im[2][7] = 0.2;
    im[7][2] = 0.2;
  }
  else if (!strcmp(puffer, "uniform_random")){
    printf("#Uniform!");
    fscanf(in,"%f ", &tmp);
    im = uniformRandomAlloc(dim,tmp);
  }
  fflush(stdout);
  printMatrix(stdout, im, dim,dim);  fflush(stdout);
  
  lattice = (ode_set***) calloc(xsize, sizeof(ode_set**));

  for (i = 0; i < xsize; i++){
    lattice[i] = (ode_set**) calloc(ysize, sizeof(ode_set*));
  }

  for(i = 1; i < xsize-1; i++){
    for(j = 1; j < ysize-1; j++){
      ode_set *os = (lattice[i][j] = createOdeSet());
      os->N = dim;
      os->functions = &spatial_ecology;
      os->data_dump = &spatial_ecology_dump;
      os->modification = &diffuse;
      os->data = spe_createLocation(i,j);
      os->state = (double*) calloc(dim, sizeof(double));
      os->nextstate = (double*) calloc(dim, sizeof(double));
      simplex(dim, os->state);
    }
  }
  
  //periodic boundaries
  for(i = 1; i < xsize-1; i++){
    lattice[i][0] = lattice[i][ysize-2];
    lattice[i][ysize-1] = lattice[i][1];
  }
  for(j = 1; j < ysize-1; j++){
    lattice[0][j] = lattice[xsize-2][j];
    lattice[xsize-1][j] = lattice[1][j];
  }
  lattice[0][0] = lattice[xsize-2][ysize-2];
  lattice[xsize-1][ysize-1] = lattice[1][1];
  lattice[0][ysize-1] = lattice[xsize-2][1];
  lattice[xsize-1][0] = lattice[1][ysize-2];

  /*
    i  = (int) xsize / 2;
    j = (int) ysize / 2;
    simplex(dim, lattice[i][j]->state);    
  */
  fclose(in);

}


location *spe_createLocation(int x, int y){
  location *nl = (location*) calloc(1, sizeof(location));
  nl->x = x;
  nl->y = y;
  return nl;
}
