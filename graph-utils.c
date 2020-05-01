//
// Created by kajetan on 30.04.2020.
//

#include "graph-utils.h"


int compareEdges(const void * a, const void * b){
    Edge * edgeA = (Edge*) a;
    Edge * edgeB = (Edge*) b;
    if(edgeA->from > edgeB->from){
        return 1;
    } else if(edgeB->from > edgeA->from){
        return -1;
    } else {
        return 0;
    }
}

void addEdge(Graph *g, int index, int from, int to, float value){
    Edge e = {.from=from - 1, .to=to - 1, .value=value};
    g->edges[index] = e;
}

void sortEdges(Graph *g){
    qsort(g->edges, g->numEdges, sizeof(Edge), compareEdges);

    int counter = 0;
    int currVertice = 0;
    while(counter < g->numEdges && currVertice < g->size){
        Edge e = g->edges[counter];
        if(e.from > currVertice){
            g->verticeLastEdgeExclusive[currVertice] = counter;
            currVertice++;
        }
        counter++;
    }
    g->verticeLastEdgeExclusive[currVertice] = counter;
}

Graph * initGraph(MData * data){
    if(data->format.format == ARRAY){
        THROW("Array matrix type not supported!", 16);
    }
    int onDiagonal = 0;
    if(data->format.symmetry != SKEW){
        for(int i=0; i<data->size; i++){
            if(data->from[i] == data->to[i]){
                onDiagonal++;
            }
        }
    }
    int edges = (data->size - onDiagonal);
    if(data->format.symmetry == SYMMETRIC || data->format.symmetry == SKEW){
        edges = edges * 2;
    }
    edges += data->rows;

    Graph *g = malloc(sizeof(Graph));
    g->verticeLastEdgeExclusive = (int*) malloc(sizeof(int) * data->rows);
    g->edges = (Edge*) malloc(sizeof(Edge) * edges);
    g->numEdges = edges;
    g->size = data->rows;

    int counter=0;
    for(int i=0; i<data->size; i++){
        addEdge(g, counter++, data->from[i], data->to[i], data->value[i]);
        if(data->format.symmetry == SKEW || data->format.symmetry == SYMMETRIC){
            addEdge(g, counter++, data->to[i], data->from[i], data->value[i]);
        }
    }

    sortEdges(g);

//    int currVertice=0;
//    for(int i=0; i<g->numEdges; i++){
//        Edge e = g->edges[i];
//        printf("%d\t->\t%d\tw=%f", e.from, e.to, e.value);
//        if(i+1 == g->verticeLastEdgeExclusive[currVertice]){
//            printf("\tEND EDGES %d\n", currVertice);
//            currVertice++;
//        } else {
//            printf("\n");
//        }
//    }

    return g;
}


float hasEdge(Graph *g, int from, int to){
    for(int i=EDGES_IDX(g, from-1); i<EDGES_IDX(g, from); i++){
        if(g->edges[i].to == to){
            return g->edges[i].value;
        }
    }
    return -1.f;
}

void printGraph(Graph *g){
    for(int i=0; i<g->size; i++){
        float val = hasEdge(g, i, 0);
        printf("%1.0f", val >= 0 ? val : 0.f);
        for(int j=1; j < g->size; j++){
            val = hasEdge(g, i, j);
            printf(" %1.0f", val >= 0 ? val : 0.f);
        }
        printf("\n");
    }
}

void destroyGraph(Graph* g){
    free(g->edges);
    free(g->verticeLastEdgeExclusive);
    free(g);
}


