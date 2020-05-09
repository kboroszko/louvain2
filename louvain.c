//
// Created by kajetan on 30.04.2020.
//

#include "louvain.h"
#include "graph-utils.h"
#include <stdbool.h>

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
        if(e.to != vertice && cliques[e.to] == in){
            sum+= e.value;
        }
    }
    return sum;
}

int bestClique(Graph *g, int vertice, int *cliques, float*sigmaTots, float m){
    float best = 0;
    int bestClique = -1;
    for(int i=EDGES_IDX(g,vertice-1); i<EDGES_IDX(g,vertice); i++){
        int to = g->edges[i].to;
        int in = cliques[to];
        if(in != cliques[vertice]){
            float deltaQ = dQ(g, vertice, cliques, in, sigmaTots[in], m);
            if(deltaQ > best){
                best = deltaQ;
                bestClique = in;
            } else if (deltaQ == best){
                if(bestClique > in){
                    bestClique = in;
                }
            }
        }
    }
    return bestClique;
}

float dQ(Graph*g, int vertice, int *cliques, int in, float sigma, float m){
    float kiin = getKiin(g, vertice, cliques, in);
    float ki = getKi(g, vertice);

    return kiin/m - (ki * sigma)/(2 * m * m);
}

int moveValid(int from, int to, int* cliqueSizes){
    if(from == to){
        return 0;
    }
    if(cliqueSizes[from] != 1){
        return 1;
    }
    if(cliqueSizes[to] > 1){
        return 1;
    }
    if(from > to){
        return 1;
    }
    return 0;
}

void printAll(Graph*g, int* cliques){
//    printf("%%graph:\n");
//    printGraph(g);
    for(int i=0; i<g->size; i++){
        printf("%d\n", cliques[i]);
    }
//    printf("%%end of graph\n");
}

void recalcSigmaTot(Graph*g, float* sigmaTot, int* cliques){
    for(int i=0; i < g->size; g++){
        sigmaTot[i] = 0;
    }
    for(int i = 0; i < g->numEdges; i++){
        Edge e = g->edges[i];
        int cl = cliques[e.to];
        sigmaTot[cl] += e.value;
    }
}

int phaseOne(Graph *g, int *cliques, float minimum){
    int changed = 1;
    int iters = 0;
    float* sigmaTot = (float*) malloc(sizeof(float) * g->size);
    int* cliqueSizes = (int*) calloc(sizeof(int), g->size);
    float m = 0;
    for(int i=0; i < g->size; i++){
        int c = cliques[i];
        cliqueSizes[c] +=1;
        sigmaTot[i] = getKi(g, i);
        m += sigmaTot[i];
    }
    m = m/2;

    printf("m=%f\n", m);
    printAll(g, cliques);
    while(changed != 0){

        printf("---------------------------- iter %d ------------------------------------------\n", iters);

        changed = 0;
        iters++;
        int* newCliques = (int*) malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);
        for(int vert=0; vert < g->size; vert++){
            int pretender = bestClique(g, vert, cliques, sigmaTot,m);
            if(pretender != -1){
                float deltaQ = dQ(g, vert, cliques, pretender, sigmaTot[pretender], m);
                if(deltaQ > minimum && moveValid(cliques[vert],pretender, cliqueSizes)){
                    printf("gonna move %2d from %2d to %2d   gain: %f \n", vert, cliques[vert], pretender, deltaQ );
                    changed = 1;
                    int oldClique = cliques[vert];
                    newCliques[vert] = pretender;
                    cliqueSizes[pretender] += 1;
                    cliqueSizes[oldClique] -= 1;
                }
            }
        }
        recalcSigmaTot(g, sigmaTot, cliques);
        memcpy(cliques, newCliques, g->size);
        //destroy newCliques
        free(newCliques);

        printAll(g, cliques);
    }

    free(sigmaTot);
    return iters;
}

float changeEdges(Graph *g, const int *cliques, const int *mins){
    for(int i=0; i<g->numEdges; i++){
        Edge *e = g->edges + i;
        int vertice = e->from;
        int cliq = cliques[vertice];
        int superVertice = mins[cliq];
        if(vertice != superVertice) {
            e->from = superVertice;
        }
        e->to = mins[cliques[e->to]];
    }
}


void phaseTwo(Graph *g, int *cliques){
    int *mins = (int*) malloc(g->size* sizeof(int)); //minimalny wierzchołek w klice
    for(int i=0; i<g->size; i++){
        mins[i] = -1;
    }
    for(int vertice=0; vertice < g->size; vertice++){
        int cl = cliques[vertice];
        if(mins[cl] == -1 || mins[cl] > vertice){
            mins[cl] = vertice;
        }
    }

    changeEdges(g, cliques, mins);

    sortEdges(g);

    float sum=0;
    int from=0;
    int to=0;
    Edge *lastEdge = g->edges;
    for(int i=0; i<g->numEdges; i++){
        Edge *e = g->edges + i;
        if(e->from != from || e->to != to){
            if(e->from > from){
                from = e->from;
            }
            to = e->to;
            if(sum > 0){
                lastEdge->value = sum;
            }
            sum = 0;
        }
        lastEdge = e;
        sum += e->value;
        e->value = 0;
    }

    sortEdges(g);
}

void moveClique(int size, int* cliques, int from, int to){
    for(int i=0; i<size; i++){
        if(cliques[i] == from){
            cliques[i] = to;
        }
    }
}

void updateOldCliques(int size, int* oldCliques, int*newCliques){
    printf("changing cliques\n");
    for(int i=0; i<size; i++){
        printf("%d -> %d\n", oldCliques[i], newCliques[i]);
    }

    for (int i = 0; i < size; ++i) {
        if(oldCliques[i] != newCliques[i]){
            moveClique(size, oldCliques, oldCliques[i], newCliques[i]);
        }
    }
}

void printCliques(int size, int*cliques){
    for (int i = 0; i < size; ++i) {
        printf("%d\n", cliques[i]);
    }
}

int main(){
    printf("hello world\n");

    MData * dat = readData("example.mtx");
    printData(dat);

    Graph *g = initGraph(dat);
    destroyMData(dat);



    int* cliques = (int*) malloc(sizeof(int) * g->size);
    for(int i=0; i<g->size; i++){
        cliques[i]=i;
    }

    for(int iter=0; iter<5; iter++){
        int* newCliques = (int*) malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);

        printf("========= PHASE 1 ==================\n");
        int iter = phaseOne(g, newCliques,0);
        if(iter == 1) {
            printf("converged!\n");
            break;
        }
//        printEdges(g);

        phaseTwo(g, newCliques);
        printf("========= PHASE 2 ==================\n");

//        printEdges(g);

        updateOldCliques(g->size, cliques, newCliques);

        free(newCliques);
    }
    printCliques(g->size, cliques);

    free(cliques);

    destroyGraph(g);



    return 0;
}