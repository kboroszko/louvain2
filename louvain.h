//
// Created by kajetan on 30.04.2020.
//

#include "graph-utils.h"

#ifndef LOUVAIN2_LOUVAIN_H
#define LOUVAIN2_LOUVAIN_H


float getKi(Graph *g, int vertice);

float getKiin(Graph *g, int vertice, int* cliques, int in );

float dQ(Graph*g, int vertice, int *clliques, int in, float sigma, float m);

void phaseOne(Graph *g, int *cliques);

#endif //LOUVAIN2_LOUVAIN_H
