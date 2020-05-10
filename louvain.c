//
// Created by kajetan on 30.04.2020.
//

#include <assert.h>
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
        if(in != bestClique && in != cliques[vertice]){
            float deltaQ = dQ(g, vertice, cliques, in, sigmaTots, m);
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

float selfLoop(Graph *g, int vertice){
    for(int k=EDGES_IDX(g,vertice-1); k<EDGES_IDX(g,vertice); k++){
        Edge e = g->edges[k];
        if(e.to == vertice){
            return e.value;
        }
    }
    return 0.0f;
}

int verticeHasEdges(Graph *g, int vertice){
    int has = EDGES_IDX(g, vertice-1) != EDGES_IDX(g, vertice);
    return has;
}

float modularity(Graph *g, int * cliques){
    float sum = 0;
    float m = 0;
    float * ac = (float*) calloc(g->size, sizeof(float));

    for(int i=0; i<g->size; i++){
        float ki = getKi(g, i);
        int clique = cliques[i];
        ac[clique] += ki;
        m += ki;
    }
    m = m/2.f;

    for(int i=0; i<g->size; i++){
        sum += ac[i] * ac[i];
    }

    sum = -sum/(2.f * m);

    for(int i=0; i < g->size; i++){
        float EiwCiBezi = getKiin(g, i, cliques, cliques[i]);
        sum += EiwCiBezi + selfLoop(g, i);
    }
    free(ac);
    return sum/(2.f*m);
}

int compareMoves( const void * a, const void * b){
    Move* ma = ((Move*)a);
    Move* mb = ((Move*)b);
    if(ma->gain > mb->gain) {
        return -1;
    } else if(mb->gain > ma->gain){
        return 1;
    } else {
        return 0;
    }
}

void applyBestMoves(int* cliques, Move* moves ,int nMoves, int nBest, int sort){
    if(nMoves == 0){
        return;
    }
    assert(nMoves >= nBest);
    if(nMoves != nBest && sort != 0){
        qsort(moves, nMoves, sizeof(Move), compareMoves);
    }
    for(int i=0; i < nBest; i++){
        Move m = moves[i];
        cliques[m.vertice] = m.toClique;
    }
}

float previewModularity(Graph * g, int*newCliques, Move* moves, int nMoves, int nBest, int sort){
    applyBestMoves(newCliques, moves, nMoves, nBest, sort);
    float newMod = modularity(g, newCliques);
    return newMod;
}

float dQ(Graph*g, int vertice, int *cliques, int in, float* sigmaTot, float m){
    float ki = getKi(g, vertice);
    float kiin = getKiin(g, vertice, cliques, in);
    float EiwCiBezi = getKiin(g, vertice, cliques, cliques[vertice]);
    float aciBezi= sigmaTot[cliques[vertice]] - ki;
    float acj = sigmaTot[in];
    float part1 = (kiin - EiwCiBezi)/m;
    float part2 = ki * (aciBezi - acj)/(2 * m * m);
    return  part1+part2;
}

int moveValid(int from, int to, int* cliqueSizes){
    if(from == to){
        return 0;
    }
    if(from > to || cliqueSizes[from] > 1 || cliqueSizes[to] > 1){
        return 1;
    }
    return 0;
}

void recalcSigmaTot(Graph*g, float* sigmaTot, int* cliques){
    for(int i=0; i < g->size; i++){
        sigmaTot[i] = 0;
    }
    for(int i = 0; i < g->size; i++){
        float ki = getKi(g, i);
        sigmaTot[cliques[i]] += ki;
    }
}

int calculateMovesToApply(int iters, int movesDone, int nMoves){
    int ret = movesDone;
    for(int i=0; i< iters; i++){
        ret = (ret + 1)/ 2;
    }
    return ret > 0 ? ret : 1;
}

int phaseOne(Graph *g, int *cliques, float minimum, float threshold){
    int changed = 1;
    int iters = 0;
    float* sigmaTot = (float*) malloc(sizeof(float) * g->size);
    int* cliqueSizes = (int*) calloc(g->size, sizeof(int));
    int nMoves = g->size;
    Move * moves = (Move*) calloc(nMoves, sizeof(Move));
    int movesDone = 0;
    float m = 0;
    for(int i=0; i < g->size; i++){
        float ki = getKi(g, i);
        if(ki > 0){
            int c = cliques[i];
            cliqueSizes[c] +=1;
            m += ki;
        }
    }
    m = m/2;
    recalcSigmaTot(g, sigmaTot, cliques);
    float mod = modularity(g, cliques);


    while(changed != 0 ){
        if(DEBUG){
            printf("---------------------------- small iter %d ------------------------------------------\n", iters);
        }
        changed = 0;
        iters++;
        movesDone = 0;
        for(int vert=0; vert < g->size; vert++){
            int pretender = bestClique(g, vert, cliques, sigmaTot,m);
            if(pretender != -1){
                float deltaQ = dQ(g, vert, cliques, pretender, sigmaTot, m);
                if(deltaQ > minimum && moveValid(cliques[vert],pretender, cliqueSizes)){
                    if(DEBUG) {
                        printf("%.8f > %.8f\n", deltaQ, minimum);
                        printf("gonna move %2d from %2d to %2d   gain: %f \n", vert, cliques[vert], pretender, deltaQ);
                    }
                    changed = 1;
                    int oldClique = cliques[vert];
                    Move * m = moves + movesDone;
                    m->vertice = vert;
                    m->gain = deltaQ;
                    m->toClique = pretender;
                    movesDone++;
                    cliqueSizes[pretender] += 1;
                    cliqueSizes[oldClique] -= 1;
                }
            }
        }
        if(DEBUG){
            int* newCliques = malloc(sizeof(int) * g->size);
            memcpy(newCliques, cliques, sizeof(int) * g->size);
            float newMod = previewModularity(g, newCliques, moves, movesDone, movesDone, 0);
            printf("modularity gain if all applied=%f\n", newMod - mod);
            free(newCliques);
        }
        int movesToApply = calculateMovesToApply(1, movesDone, nMoves);

        int* newCliques = malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);
        float newMod = previewModularity(g, newCliques, moves, movesDone, movesToApply, 1);

        if(DEBUG){
            printf("modularity gain if %d applied=%f\n",movesToApply, newMod - mod);
        }


        if(movesDone > 0){
            float bestdQ = moves[0].gain;
            int movesIter = 2;
            while((newMod - mod < threshold) && (movesToApply > 1 || bestdQ > threshold)){
                movesToApply = calculateMovesToApply(movesIter, movesDone, nMoves);
                memcpy(newCliques, cliques, sizeof(int) * g->size);
                newMod = previewModularity(g, newCliques, moves, movesDone, movesToApply, 0);
                movesIter++;
            }
            if (newMod - mod > threshold) {
                memcpy(cliques, newCliques, sizeof(int) * g->size);
                recalcSigmaTot(g, sigmaTot, cliques);
                mod = newMod;
//                printf("%f, \n", modularity(g, cliques));
            }
            if(movesToApply == 1 && bestdQ < threshold){
                changed = 0;
            }
        } else {
            changed = 0;
        }
        free(newCliques);
    }
    free(cliqueSizes);
    free(moves);
    free(sigmaTot);
    return iters;
}

/**
 * update all edges to go to and from superVertices aka cliques changed to vertices in phase2
 */
void changeEdges(Graph *g, const int *cliques){
    for(int i=0; i<g->numEdges; i++){
        Edge *e = g->edges + i;
        int vertice = e->from;
        int cliq = cliques[vertice];
        int superVertice = cliq;
        if(vertice != superVertice) {
            e->from = superVertice;
        }
        e->to = cliques[e->to];
    }
}

void changeCliqueToMin(Graph *g, int*cliques){
    int *mins = (int*) malloc(g->size* sizeof(int)); //minimalny wierzchołek w klice
    for(int i=0; i<g->size; i++){
        mins[i] = -1;
    }
    for(int vertice=0; vertice < g->size; vertice++){
        if(verticeHasEdges(g, vertice)){
            int cl = cliques[vertice];
            if(mins[cl] == -1 || mins[cl] > vertice){
                mins[cl] = vertice;
            }
        }
    }
    for(int i=0; i<g->size; i++){
        int cl = cliques[i];
        if(mins[cl] != -1){
            cliques[i] = mins[cl];
        }
    }
    free(mins);
}

void phaseTwo(Graph *g, int *cliques){
    changeCliqueToMin(g, cliques);
    changeEdges(g, cliques);
    sortEdges(g);

    float sum=0;
    int from=0;
    int to=0;
    //aggregate edges
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

/**
 * oldCliques represents vertice to clique mapping, this function updates it
 * when a super vertice changed place
 * @param size
 * @param oldCliques
 * @param newCliques
 */
void updateOldCliques(Graph *g, int* cliques){
    for(int i=0; i<g->size; i++){
        int index = cliques[i];
        while(index != cliques[index]){
            index = cliques[index];
        }
        if(cliques[i] != index){
            if(DEBUG){
            printf("CHANGED CLIQUE FOR %d FROM %d TO %d\n", i, cliques[i], index);
            }
            cliques[i] = index;
        }
    }
}

void printCliques(int size, int*cliques){
    printf("labs=[");
    for (int i = 0; i < size; ++i) {
//        printf("cliques[%d]=%d;\n", i, cliques[i]);
        printf("%d,", cliques[i]);
        if(i%500 == 0){
            printf("\n");
        }
    }
    printf("];\n");
}


void printUsage(char * name){
    printf("# Usage:\n");
    printf("# %s  [--verbose] <filename>\n", name);
    printf("#     --verbose   print out the links\n");
    printf("#     filename    name of file with MTX matrix\n");
}

int main(int argc, char **argv){
    char * fileName;
    int verbose = 0;
    if(argc < 2 || argc > 3){
        printf("wrong number of arguments!");
        printUsage(argv[0]);
        return 1;
    } else if(argc == 2){
        fileName = argv[1];
    } else {
        if(strcmp(argv[1], "--verbose") == 0){
            fileName = argv[2];
            verbose = 1;
        } else {
            printUsage(argv[0]);
            return 2;
        }
    }

    MData * dat = readData(fileName);

    Graph *g = initGraph(dat);
    destroyMData(dat);

    int* cliques = (int*) malloc(sizeof(int) * g->size);
    for(int i=0; i<g->size; i++){
        cliques[i]=i;
    }

    int bigLoopIteration = 0;
    float minimum = 0.1 / (2 + bigLoopIteration) - 0.02;

    float threshold = 0.00001f;

    // profiler at hangGlider_4 th=0.00001f

    float mod = modularity(g, cliques);
    printf("modularity:%f\n", mod);
    int iter = 10;
    while(iter > 1 || minimum > threshold/10.f){

//        printf("========= PHASE 1 ==================\n");
        minimum = 0.1 / (2 + bigLoopIteration) - 0.02;
        minimum = minimum < threshold/20.f ? threshold/20.f : minimum;
//        printf("min:%f\n", minimum);
        iter = phaseOne(g, cliques, minimum, threshold);

//        printCliques(g->size, cliques);

//        printf("========= PHASE 2 ==================\n");
        phaseTwo(g, cliques);
//        printEdges(g);
        updateOldCliques(g, cliques);
//        printf("modularity:%f\n", modularity(g, cliques));
//        printCliques(g->size, cliques);
        bigLoopIteration += 1;
    }
    printf("converged after %d iterations!\n", bigLoopIteration+1);
    if(verbose != 0){
        printCliques(g->size, cliques);
    }


    mod = modularity(g, cliques);
    printf("modularity:%f\n", mod);

    free(cliques);

    destroyGraph(g);



    return 0;
}