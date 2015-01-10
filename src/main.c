#include "ode_solver.h"

int main(int argc, char **argv){
  FILE *in; 
  ode_set **des;

  //sanity check
  if(argc != 2){
    printf("Usage: ./ode inputfile\n");
    exit(1);
  }
  in = (FILE*) fopen(argv[1], "r");
  if (!in){printf("File could not be opened!\n"); exit(1);}

  //reading the inputfile
  des = load_ode_set(in);
  //solving the equations and dumping information 
  solve(des);
  return 0;
}
