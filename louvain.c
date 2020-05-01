//
// Created by kajetan on 30.04.2020.
//

#include "louvain.h"
#include "graph-utils.h"

float getKi(Graph *g, int vertice){
    float sum=0;
    for(int i=EDGES_IDX(g,vertice-1); i<EDGES_IDX(g,vertice); i++){
        sum+= g->edges[i].value;
    }
    return sum;
}

float getKiin(Graph *g, int vertice, int* cliques, int in ){
    float sum=0;
    for(int i=EDGES_IDX(g,vertice-1); i<EDGES_IDX(g,vertice); i++){
        Edge e = g->edges[i];
        if(cliques[e.to] == in)
            sum+= e.value;
    }
    return sum;
}

int bestClique(Graph *g, int vertice, int *cliques){
    float best = 0;
    int bestVert = -1;
    for(int i=EDGES_IDX(g,vertice-1); i<EDGES_IDX(g,vertice); i++){
        int to = g->edges[i].to;
        int in = cliques[to];
        float kiin = getKiin(g, vertice, cliques, in);
        if(kiin > best){
            best = kiin;
            bestVert = to;
        }
    }
    return bestVert;
}

float dQ(Graph*g, int vertice, int *cliques, int in, float sigma, float m){
    float kiin = getKiin(g, vertice, cliques, in);
    float ki = getKi(g, vertice);

    return kiin/m - (ki * sigma)/(2 * m * m);
}

void phaseOne(Graph *g, int *cliques){
    int changed = 1;
    while(changed != 0){
        changed = 0;

        for(int vert=0; vert < g->size; g++){

        }


    }
}


int main(){
    printf("hello world\n");

    MData * dat = readData("mycielskian4.mtx");
    printData(dat);

    Graph *g = initGraph(dat);
    destroyMData(dat);

    printGraph(g);

    destroyGraph(g);



    return 0;
}