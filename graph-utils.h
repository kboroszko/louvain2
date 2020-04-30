//
// Created by kajetan on 30.04.2020.
//

#ifndef LOUVAIN2_GRAPH_UTILS_H
#define LOUVAIN2_GRAPH_UTILS_H

#include "mmio.h"

typedef struct{
    int from;
    int to;
    float value;
} Edge;

typedef struct {
    int size;
    int numEdges;
    Edge* edges;
    int * verticeLastEdgeExclusive;
} Graph;

void sortEdges(Graph *g);

int compareEdges(const void * a, const void * b);

Graph * initGraph(MData * data);

void destroyGraph(Graph* g);

void printGraph(Graph *g);

float hasEdge(Graph *g, int from, int to);

void addEdge(Graph *g, int index, int from, int to, float value);

float getKi(Graph *g, int vertice);

float getKiin(Graph *g, int vertice, int* cliques, int in );

#endif //LOUVAIN2_GRAPH_UTILS_H
