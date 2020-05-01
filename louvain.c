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
        if(cliques[e.to] == in){
            sum+= e.value;
        }
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
        } else if (kiin == best){
            if(bestVert > to){
                bestVert = to;
            }
        }
    }
    return bestVert;
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
//    printGraph(g);
    for(int i=0; i<g->size; i++){
        printf("%d\n", cliques[i]);
    }
}

void phaseOne(Graph *g, int *cliques, float minimum){
    int changed = 1;
    int iters = 0;
    float* sigmaTot = (float*) malloc(sizeof(float) * g->size);
    int* cliqueSizes = (int*) malloc(sizeof(int) * g->size);
    float m = 0;
    for(int i=0; i < g->size; i++){
        cliqueSizes[i] = 1;
        sigmaTot[i] = getKi(g, i);
        m += sigmaTot[i];
    }
    while(changed != 0){

        printf("iter - %d\n", iters);
        printf("--------------------------------------------------------------------------------\n");
        printAll(g, cliques);

        changed = 0;
        iters++;
        int* newCliques = (int*) malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);
        for(int vert=0; vert < g->size; vert++){
            int pretender = bestClique(g, vert, cliques);
            if(pretender != -1){
                float deltaQ = dQ(g, vert, cliques, pretender, sigmaTot[pretender], m);
                if(deltaQ > minimum && moveValid(cliques[vert],pretender, cliqueSizes)){
                    changed = 1;
                    newCliques[vert] = pretender;
                }
            }
        }
        //first update sigmaTot
        for(int i=0; i<g->size; i++){
            if(cliques[i] != newCliques[i]){
                printf("moving %2d from %2d to %2d\n", i, cliques[i], newCliques[i] );
                float ki = getKi(g, i);
                float kiin = getKiin(g, i, cliques, newCliques[i]);
                float kiinOld = getKiin(g, i, newCliques, cliques[i]);
                sigmaTot[newCliques[i]] += ki - kiin;
                sigmaTot[cliques[i]] -= ki - kiinOld;
                cliqueSizes[cliques[i]] -= 1;
                cliqueSizes[newCliques[i]] += 1;
            }
        }
        //update cliques
        memcpy(cliques, newCliques, sizeof(int) * g->size);
        //destroy newCliques
        free(newCliques);
    }

    free(sigmaTot);
}


int main(){
    printf("hello world\n");

    MData * dat = readData("mycielskian4.mtx");
    printData(dat);

    Graph *g = initGraph(dat);
    destroyMData(dat);

//    printGraph(g);
    int* cliques = (int*) malloc(sizeof(int) * g->size);
    for(int i=0; i<g->size; i++){
        cliques[i]=i;
    }
    phaseOne(g, cliques,0);

    destroyGraph(g);



    return 0;
}