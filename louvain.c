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
    float * ac = (float*) calloc(sizeof(float), g->size);

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

void applyBestMoves(int* cliques, Move* moves ,int nMoves, int nBest){
    if(nMoves == 0){
        return;
    }
    assert(nMoves >= nBest);
    if(nMoves != nBest){
        qsort(moves, nMoves, sizeof(Move), compareMoves);
    }
    for(int i=0; i < nBest; i++){
        Move m = moves[i];
        cliques[m.vertice] = m.toClique;
    }
}

float previewModularity(Graph * g, int*newCliques, Move* moves, int nMoves, int nBest){
    applyBestMoves(newCliques, moves, nMoves, nBest);
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
    for(int i=1; i< iters; i++){
        ret = (ret + 1)/ 2;
    }
    return ret > 0 ? ret : 1;
}

int phaseOne(Graph *g, int *cliques, float minimum, float threshold){
    int changed = 1;
    int iters = 0;
    float* sigmaTot = (float*) malloc(sizeof(float) * g->size);
    int* cliqueSizes = (int*) calloc(sizeof(int), g->size);
    int nMoves = g->numEdges;
    Move * moves = (Move*) calloc(sizeof(Move) , nMoves);
    int movesDone = 0;
    printf("moves:%p\n", moves);
    float m = 0;
    for(int i=0; i < g->size; i++){
        int c = cliques[i];
        cliqueSizes[c] +=1;
        m += getKi(g, i);
    }
    m = m/2;

    recalcSigmaTot(g, sigmaTot, cliques);

//    printf("sigmatot:\n");
//    for(int i=0; i<g->size; i++){
//        printf("%f\n", sigmaTot[i]);
//    }

//
    float mod = modularity(g, cliques);
    printf("mod:%f\n", mod);
    float dq = dQ(g, 10, cliques, 10, sigmaTot, m);
    printf("dq=%f\n", dq);
    cliques[10] = 10;
    recalcSigmaTot(g, sigmaTot, cliques);

//    printf("sigmatot:\n");
//    for(int i=0; i<g->size; i++){
//        printf("%f\n", sigmaTot[i]);
//    }


    float dq2 = dQ(g, 10, cliques, 9, sigmaTot, m);
    printf("dq2=%f\n", dq2);

    printf("mod:%f\n", modularity(g, cliques));
    printf("delta:%f\n", modularity(g, cliques) - mod);





    while(changed != 0 ){
        printf("---------------------------- iter %d ------------------------------------------\n", iters);
        printf("%f, \n", modularity(g, cliques));
        changed = 0;
        iters++;
        movesDone = 0;
        for(int vert=0; vert < g->size; vert++){
            int pretender = bestClique(g, vert, cliques, sigmaTot,m);
            if(pretender != -1){
                float deltaQ = dQ(g, vert, cliques, pretender, sigmaTot, m);
                if(deltaQ > minimum && moveValid(cliques[vert],pretender, cliqueSizes)){
                    printf("gonna move %2d from %2d to %2d   gain: %f \n", vert, cliques[vert], pretender, deltaQ );
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
            float newMod = previewModularity(g, newCliques, moves, movesDone, movesDone);
            printf("modularity gain if all applied=%f\n", newMod - mod);
            free(newCliques);
        }
        int movesToApply = calculateMovesToApply(iters, movesDone, nMoves);

        int* newCliques = malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);
        float newMod = previewModularity(g, newCliques, moves, movesDone, movesToApply);

        if(DEBUG){
            printf("modularity gain if %d applied=%f\n",movesToApply, newMod - mod);
        }


        if(movesDone > 0 ){
            if (newMod - mod > threshold) {
                memcpy(cliques, newCliques, sizeof(int) * g->size);
                recalcSigmaTot(g, sigmaTot, cliques);
                mod = newMod;
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
 * update all edges to go to and from superVertices aka cliques changed to vertices after phase2
 */
void changeEdges(Graph *g, const int *cliques, const int *mins){
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
    free(mins);
}

/**
 * change all vertices in clique "from" to be in vertice "to"
 * @param size
 * @param cliques
 * @param from
 * @param to
 */
void moveClique(int size, int* cliques, int from, int to){
    for(int i=0; i<size; i++){
        if(cliques[i] == from){
            cliques[i] = to;
        }
    }
}

/**
 * oldCliques represents vertice to clique mapping, this function updates it
 * when a super vertice changed place
 * @param size
 * @param oldCliques
 * @param newCliques
 */
void updateOldCliques(int size, int* oldCliques, int*newCliques){
//    printf("changing cliques\n");
//    for(int i=0; i<size; i++){
//        printf("%d -> %d\n", oldCliques[i], newCliques[i]);
//    }

    for (int i = 0; i < size; ++i) {
        if(oldCliques[i] != newCliques[i]){
            moveClique(size, oldCliques, oldCliques[i], newCliques[i]);
        }
    }
}

void printCliques(int size, int*cliques){
    for (int i = 0; i < size; ++i) {
        printf("cliques[%d]=%d;\n", i, cliques[i]);
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

    cliques[0]=0;
    cliques[1]=1;
    cliques[2]=1;
    cliques[3]=0;
    cliques[4]=1;
    cliques[5]=0;
    cliques[6]=6;
    cliques[7]=0;
    cliques[8]=9;
    cliques[9]=9;
    cliques[10]=9;
    cliques[11]=10;
    cliques[12]=9;
    cliques[13]=10;
    cliques[14]=9;
    cliques[15]=8;


    float mod = modularity(g, cliques);
    printf("modularity:%f\n", mod);

    for(int iter=0; iter<5; iter++){
        int* newCliques = (int*) malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);

        printf("========= PHASE 1 ==================\n");
        int iter = phaseOne(g, newCliques,0.f, 0.0f);
        if(iter == 1) {
            printf("converged!\n");
            free(newCliques);
            break;
        }

        printf("========= PHASE 2 ==================\n");
        phaseTwo(g, newCliques);
        printf("modularity:%f\n", modularity(g, cliques));
        updateOldCliques(g->size, cliques, newCliques);
        free(newCliques);
    }
    printCliques(g->size, cliques);

    free(cliques);

    destroyGraph(g);



    return 0;
}