#include "spatial_replicator.h"
#include "matrices.h"

//for dumping pictures
#include "png_util.h"

#include <stdlib.h>

typedef struct{
  int x;
  int y;
} location;

///interaction matrix
double **im;

double diffusion_rate;
int xsize;
int ysize;
ode_set ***lattice;
char flag = 1;
int dim; 

//for png writing
int serial = 0;
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

kek-sarga
png_byte palette[30] = {253,67,25,
			252,98,21,
			252,155,10,
			252,204,10,
			251,243,10,
			0,202,251,
			23,168,252,
			23,132,250,
			23,84,252,
			23,27,252};

*/

png_byte palette[96] = {
  2,2,41,
  0,0,84,
  0,0,124,
  0,0,157,
  0,0,220,
  0,0,255,
  0,43,255,
  0,85,255,
  0,128,255,
  0,170,255,
  0,213,255,
  0,255,255,
  0,255,213,
  0,255,170,
  0,255,128,
  0,255,85,
  0,255,43,
  0,255,0,
  43,255,0,
  85,255,0,
  128,255,0,
  170,255,0,
  213,255,0,
  255,255,0,
  255,213,0,
  255,170,0,
  255,128,0,
  255,85,0,
  255,43,0,
  255,0,0,
  255,0,4,
  255,0,30
};
/*
  png_byte palette[57] = {
  255,255,255,
  220,255,218,
  195,255,188,
  169,255,154,
  144,255,122,
  118,255,87,
  94,255,54,
  32,222,3,
  24,204,0,
  25,172,2,
  22,154,2,
  16,116,1,
  12,84,1,
  3,20,0,
  49,50,49,
  89,90,89,
  123,124,123,
  165,167,165,
  207,209,207 };
*/
png_byte *image = 0;
png_byte *density_image = 0;
int density_species; //which species should be displayed
int width;
int height;
char picname[CHAR_BUFF_SIZE];

void spat_load_parameters(ode_set *des,FILE *in);
double spat_score(int i, double y[]);
double spat_total(double y[]);
void spat_replicators(int N, double x, double y[], double result[], void *data);
void spat_modify(void *des);
void spat_dump(FILE *in,double time, void *des);
void spat_dump_dominant(FILE *in,double time, void *des);
void spat_dump_density(FILE *in,double time, void *des);
void spat_rep_createlattice();
location *createLocation(int x, int y);
ode_set *createOdeSet();

ode_set **createSpatialReplicatorODE(ode_set *des,FILE *in){
  ode_set **os_array;
  int i;
  int  j;
  int  k;
  
  dim = des->N;
  spat_load_parameters(des,in);
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

void spat_replicators(int N, double x, double y[], double result[],void *data){
  double tot; 
  int i;
  tot = spat_total(y);
  for (i = 0; i < N; i++){
    result[i] = y[i]*(spat_score(i,y) - tot);
  }
}

double spat_score(int  n, double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += (y[i]) * (im[n][i]);
  }   
  return  sum;
}

double spat_total(double y[]){
  int i;
  double sum = 0.0;
  for (i = 0; i < dim; i++){
    sum += y[i] * spat_score(i,y);
  }
  return sum;
}

void spat_modify(void *o_s){
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
//local dynamics
void spat_dump(FILE *out, double time,void *o_s){
  int i;
  ode_set *os = (ode_set*) o_s;
  
  fprintf(out, "%f ", time);
  for (i = 0; i < dim; i++){
    fprintf(out,"%f ", os->state[i]);
  }
  fprintf(out, "\n");
  fflush(out);    
}

void spat_dump_dominant(FILE *out, double time,void *o_s){
  int i,j,k,imax;
  double max;
  ode_set *tmp;

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
      
      image[((j-1)*width)+(i-1)] = imax;
    }
  }
  sprintf(picname,"dominant%6d.png", serial++);
  for (i = 0; i < strlen(picname);i++){
    if (picname[i] == ' ') picname[i] = '0';
  }
  write_indexed_png(picname, image, width, height, palette, dim);
  
}

void spat_dump_density(FILE *out, double time,void *o_s){
  int i,j;
  double max,min;
  ode_set *tmp;
  double part;
  max = 0;
  min = 1;
  /*  
      for (i = 1; i < xsize-1; i++){
      for (j = 1; j < ysize-1; j++){
      tmp = lattice[i][j];
      if (tmp->state[density_species] < min){
      min = tmp->state[density_species];
      }
      if (tmp->state[density_species] > max){
      max = tmp->state[density_species];
      }
      }
      }
  */
  part = (max-min) / (double) 31;
  for (i = 1; i < xsize-1; i++){
    for (j = 1; j < ysize-1; j++){
      tmp = lattice[i][j];
      density_image[((j-1)*width)+(i-1)] = (png_byte) ((max - tmp->state[density_species]) / part);
    }
  }
  sprintf(picname,"density%6d.png", serial++);
  for (i = 0; i < strlen(picname);i++){
    if (picname[i] == ' ') picname[i] = '0';
  }
  write_indexed_png(picname, density_image, width, height, palette, 32);
  
}


void spat_load_parameters(ode_set *des,FILE *in){
  float tmp;
  int seed;
  char puffer[CHAR_BUFF_SIZE];

  while( fscanf(in, " %s ", puffer) != EOF){
    if (!strcmp(puffer, "INTERACTION_MATRIX")){
      fprintf(stdout,"#Type of interaction matrix:");
      //identifying the type of the interaction matrix
      fscanf(in, "%s ", puffer);
      if (!strcmp(puffer, "hypercycle")){
	printf("#Hypercycle!\n");
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
	//    im[2][7] = 0.2;
	//im[7][2] = 0.2;
      }
      else if (!strcmp(puffer, "uniform_random")){
	printf("#Uniform!\n");
	fscanf(in,"%f ", &tmp);
	im = uniformRandomAlloc(dim,tmp);
      }
      else if (!strcmp(puffer, "custom")){
	printf("#Custom!\n");
	im = safeMatrixAlloc(dim,dim);
	readMatrix(in, im, dim,dim);
      }
      
    }
    else if (!strcmp(puffer, "DUMP")){
      //identifying the dumping method
      fscanf(in, "%s ", puffer);
      if (!strcmp(puffer,"dominant")){
	int x,y;
	fscanf(in, "%d %d ", &x, &y);
	lattice[x][y]->data_dump = &spat_dump_dominant;
	image = (png_byte*) calloc((width * height), sizeof(png_byte));
      }
      else if (!strcmp(puffer,"density")){
	int x,y;
	fscanf(in, "%d %d %d", &x, &y, &density_species);
	lattice[x][y]->data_dump = &spat_dump_density;
	density_image = (png_byte*) calloc((width * height), sizeof(png_byte));
      }
      else if (!strcmp(puffer,"local")){
	int x,y;
	fscanf(in, "%d %d ", &x, &y);
	lattice[x][y]->data_dump = &spat_dump;
      }
    }
    else if (!strcmp(puffer, "LATTICE_SIZE")){
      fscanf(in,"%d %d ", &xsize, &ysize);
      width = xsize -2;
      height = ysize -2;
      spat_rep_createlattice();
    }
    else if (!strcmp(puffer, "DIFFUSION_RATE")){
      fscanf(in,"%f ", &tmp);
      diffusion_rate = tmp;

    }
    else if (!strcmp(puffer, "RANDOM_SEED")){
      fscanf(in,"%d ", &seed);
      srand(seed);
    }
  }

  fflush(stdout);
  printMatrix(stdout, im, dim,dim);  fflush(stdout);
  fclose(in);

}

void spat_rep_createlattice(){
  int i,j;
  lattice = (ode_set***) calloc(xsize, sizeof(ode_set**));

  for (i = 0; i < xsize; i++){
    lattice[i] = (ode_set**) calloc(ysize, sizeof(ode_set*));
  }

  for(i = 1; i < xsize-1; i++){
    for(j = 1; j < ysize-1; j++){
      ode_set *os = (lattice[i][j] = createOdeSet());
      os->N = dim;
      os->functions = &spat_replicators;
      os->data_dump = NULL;
      os->modification = &spat_modify;
      os->data = createLocation(i,j);
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


}

location *createLocation(int x, int y){
  location *nl = (location*) calloc(1, sizeof(location));
  nl->x = x;
  nl->y = y;
  return nl;
}


