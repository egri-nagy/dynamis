#ifndef REPLICATOR_H
#define REPLICATOR_H

#include  "ode_solver.h"

void createReplicatorODE(ode_set *des, FILE* in);

void replicators(int N, double x, double y[], double result[], void *data);

#endif
//REPLICATOR_H
