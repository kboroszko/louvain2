//
// Created by kajetan on 30.04.2020.
//

#include "graph-utils.h"


int compareEdges(const void * a, const void * b){
    Edge * edgeA = (Edge*) a;
    Edge * edgeB = (Edge*) b;
    if(edgeA->value == 0 || edgeB->value ==0){
        if(edgeA->value != 0){
            return -1;
        } else if(edgeB->value != 0){
            return 1;
        }
        return 0;
    }
    if(edgeA->from > edgeB->from){
        return 1;
    } else if(edgeB->from > edgeA->from){
        return -1;
    } else {
        if(edgeA->to > edgeB->to){
            return 1;
        } else {
            return -1;
        }
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
        if(e.value == 0){
            while(currVertice < g->size){
                g->verticeLastEdgeExclusive[currVertice] = counter;
                currVertice++;
            }
            counter++;
            break;
        }
        if(e.from > currVertice){
            while(currVertice < e.from){
                g->verticeLastEdgeExclusive[currVertice] = counter;
                currVertice++;
            }
        }
        counter++;
    }
    if(counter >= g->numEdges && currVertice < g->size){
        g->verticeLastEdgeExclusive[currVertice] = counter;
    } else {
        g->numEdges = counter;
    }
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
    int edges = (data->size);
    if(data->format.symmetry == SYMMETRIC || data->format.symmetry == SKEW){
        edges = edges * 2;
    }
    edges += data->rows;

    Graph *g = calloc(1, sizeof(Graph));
    g->verticeLastEdgeExclusive = (int*) calloc(data->rows, sizeof(int));
    g->edges = (Edge*) calloc(edges, sizeof(Edge));
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

    return g;
}

void printEdges(Graph *g){

    int currVertice=0;
    for(int i=0; i<g->numEdges; i++){
        Edge e = g->edges[i];
        printf("%d\t->\t%d\tw=%f", e.from, e.to, e.value);
        while(i+1 == g->verticeLastEdgeExclusive[currVertice]){
            printf("\tEND EDGES %d", currVertice);
            currVertice++;
        }
        printf("\n");
        if(currVertice >= g->size){
            break;
        }
    }
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


