//
// Created by kajetan on 30.04.2020.
//

#include "graph-utils.h"

#ifndef LOUVAIN2_LOUVAIN_H
#define LOUVAIN2_LOUVAIN_H

#ifndef DEBUG
#define DEBUG 0 //TODO
#endif

typedef struct {
    int vertice;
    int toClique;
    float gain;
} Move;

typedef struct {
    int clique;
    int vertice;
} Pair;

int cmpPair(const void *a, const void *b){
    Pair* ap = (Pair*) a;
    Pair* bp = (Pair*) b;
    if(ap->clique > bp->clique){
        return 1;
    } else if(ap->clique < bp->clique){
        return -1;
    } else if(ap->vertice > bp->vertice){
        return 1;
    }else if(ap->vertice < bp->vertice){
        return -1;
    }else{
        return 0;
    }
}


float getKi(Graph *g, int vertice);

float getKiin(Graph *g, int vertice, int* cliques, int in );

int bestClique(Graph *g, int vertice, int *cliques, float*sigmaTots, float m);

float dQ(Graph*g, int vertice, int *clliques, int in, float* sigma, float m);

int phaseOne(Graph *g, int *cliques, float minimum, float threshold);



void phaseTwo(Graph *g, int *cliques);

float modularity(Graph *g, int * cliques);

void updateOldCliques(Graph *g, int* cliques);

#endif //LOUVAIN2_LOUVAIN_H
