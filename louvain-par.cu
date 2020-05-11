//
// Created by kajetan on 30.04.2020.
//
extern "C" {
    #include <assert.h>
    #include "louvain.h"
    #include "graph-utils.h"
}

#include "errors.h"
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

//__device__ float atomicAdd(float* address, float val)
//{
//    unsigned int* address_as_ull =
//            (unsigned int*)address;
//    unsigned int old = *address_as_ull, assumed;
//
//    do {
//        assumed = old;
//        old = atomicCAS(address_as_ull, assumed,
//                        __float_as_int(val +
//                                               __int_as_float(assumed)));
//
//        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
//    } while (assumed != old);
//
//    return __int_as_float(old);
//}

__device__ int moveValid(int from, int to, int* cliqueSizes);

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
            float deltaQ =  0;//dQ(g, vertice, cliques, in, sigmaTots, m);
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

int compareMovesThrust( Move a, Move b){
    if(a.gain > b.gain) {
        return -1;
    } else if(b.gain > a.gain){
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


__device__ int moveValid(int from, int to, int* cliqueSizes){
    if(from == to){
        return 0;
    }
    if(from > to || cliqueSizes[from] > 1 || cliqueSizes[to] > 1){
        return 1;
    }
    return 0;
}



int calculateMovesToApply(int iters, int movesDone, int nMoves){
    int ret = movesDone;
    for(int i=0; i< iters; i++){
        ret = (ret + 1)/ 2;
    }
    return ret > 0 ? ret : 1;
}

__device__ float getKiDevice(int numEdges, Edge* edges){
    float sum = 0;
    for(int i=0; i<numEdges; i++){
        sum += edges[i].value;
    }
    return sum;
}

void copyGraphToDevice(Graph*g, Graph**deviceGraphPtr){

    Edge * edgesPtr ;
    int * vertPtr ;

    HANDLE_ERROR(cudaMalloc((void**) &edgesPtr, sizeof(Edge) * g->numEdges));
    HANDLE_ERROR(cudaMalloc((void**) &vertPtr, sizeof(int) * g->size));

//    printf("graph tables malloc succeded\n");


    HANDLE_ERROR(cudaMemcpy((void*) edgesPtr, (void*)g->edges, sizeof(Edge) * g->numEdges, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy((void*) vertPtr, (void*)g->verticeLastEdgeExclusive, sizeof(int) * g->size, cudaMemcpyHostToDevice));

//    printf("copying succeded\n");


    HANDLE_ERROR(cudaMalloc((void**)deviceGraphPtr, sizeof(Graph)));

    Graph gr = {.size=g->size, .numEdges=g->numEdges, .edges=edgesPtr, .verticeLastEdgeExclusive=vertPtr};

    HANDLE_ERROR(cudaMemcpy((void*)*deviceGraphPtr, (void*)&gr, sizeof(Graph), cudaMemcpyHostToDevice));

//    printf("graph init succeded\n");
}

__device__ float getKiinDevice(Graph *g, int vertice, int* cliques, int in ){
    float sum=0;
    for(int i=EDGES_IDX(g,vertice-1); i<EDGES_IDX(g,vertice); i++){
        Edge e = g->edges[i];
        if(e.to != vertice && cliques[e.to] == in){
            sum+= e.value;
        }
    }
    return sum;
}



__device__ float dQDevice(Graph*g, int vertice, int *cliques, int in, float* sigmaTot, float m, int numEdges, Edge* edges){

    float ki = getKiDevice(numEdges, edges);
    float kiin = getKiinDevice(g, vertice, cliques, in);
    float EiwCiBezi = getKiinDevice(g, vertice, cliques, cliques[vertice]);
    float aciBezi= sigmaTot[cliques[vertice]] - ki;
    float acj = sigmaTot[in];
    float part1 = (kiin - EiwCiBezi)/m;
    float part2 = ki * (aciBezi - acj)/(2 * m * m);
    return  part1+part2;
}




void copyArrayToDevice(int * arr, int size, int** deviceArray){
    HANDLE_ERROR(cudaMalloc((void**) deviceArray, sizeof(int) * size));
    HANDLE_ERROR(cudaMemcpy((void*) *deviceArray, (void*)arr, sizeof(int) * size, cudaMemcpyHostToDevice));
}

void copyFloatArrayToDevice(float * arr, int size, float** deviceArray){
    HANDLE_ERROR(cudaMalloc((void**) deviceArray, sizeof(float) * size));
    HANDLE_ERROR(cudaMemcpy((void*) *deviceArray, (void*)arr, sizeof(float) * size, cudaMemcpyHostToDevice));
}

__global__ void recalcSigmaTotPar(Graph*g, float* sigmaTot, int* cliques) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int vertice = tid;
    if(tid < g->size){ // there is a chance that a whole block except one thread will be stuck but oh, well
        int clique = cliques[vertice];

        int edgesStart =  EDGES_IDX(g, vertice - 1);
        int edgesEnd =  EDGES_IDX(g, vertice);
        Edge * edgesPtr = g->edges + edgesStart;
        int numEdges = edgesEnd - edgesStart;
        float ki = getKiDevice(numEdges, edgesPtr);
        atomicAdd(sigmaTot + clique, ki);
    }
}

__global__ void calculateCliqueSizes(Graph*g, int* cliques, int * cliqueSizes) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid < g->size){
        int vertice = tid;
        int clique = cliques[vertice];
        atomicAdd(cliqueSizes + clique, 1);
    }
}

__global__ void calcNeighbours(Graph *g, int *sizes){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int vertice = tid;
    if(vertice < g->size){
        int edgesStart =  EDGES_IDX(g, vertice - 1);
        int edgesEnd =  EDGES_IDX(g, vertice);
        int numEdges = edgesEnd - edgesStart;
        sizes[vertice] = numEdges;
    }
}



__global__ void calculateMoves(Graph *g, int* cliques, int*cliqueSizes,
        Move* moves, float m, float* sigmaTot,
        float minimum, int * nMoves){
    extern __shared__ float bestOutcomes[];
    int vertice = blockIdx.x;
    int edgesStart =  EDGES_IDX(g, vertice - 1);
    int edgesEnd =  EDGES_IDX(g, vertice);
    Edge * edgesPtr = g->edges + edgesStart;
    int numEdges = edgesEnd - edgesStart;
    int tid = threadIdx.x;
    bestOutcomes[tid + blockDim.x] = -1.f;
    bestOutcomes[tid] = 0;
    int cliqueFrom = cliques[vertice];
    if(tid < numEdges){
        Edge e = g->edges[edgesStart + tid];
        if(cliqueFrom != cliques[e.to]){
            int pretender = cliques[e.to];
            float deltaQ = dQDevice(g, vertice, cliques, pretender, sigmaTot, m, numEdges, edgesPtr);

            if(deltaQ > minimum && moveValid(cliqueFrom, pretender, cliqueSizes)){
                bestOutcomes[tid + blockDim.x] = __int2float_rn(pretender);
                bestOutcomes[tid] = deltaQ;
            }
        }
    }
    //reduce within a block
    for (int stride=1;stride<blockDim.x;stride*=2)
    {
        __syncthreads();
        if (tid%(2*stride)==0){
            if(tid+stride < blockDim.x && bestOutcomes[tid] < bestOutcomes[tid+stride]){
                bestOutcomes[tid] = bestOutcomes[tid+stride];
                bestOutcomes[tid + blockDim.x] = bestOutcomes[tid + blockDim.x + stride];
            }
        }
    }
    if (tid==0 && bestOutcomes[0] > 0){
        int toClique = __float2int_rn(bestOutcomes[blockDim.x]);
        float gain =  bestOutcomes[0];
        int myMove = atomicAdd(nMoves, 1) - 1;
        Move m = {.vertice=vertice, .toClique = toClique, .gain=gain};
        moves[myMove] = m;

    }
}



void destroyDeviceGraph(Graph * deviceGraph){
    Graph g;
    HANDLE_ERROR(cudaMemcpy((void*) &g, (void*)deviceGraph, sizeof(Graph), cudaMemcpyDeviceToHost));
    cudaFree(g.edges);
    cudaFree(g.verticeLastEdgeExclusive);
    cudaFree(deviceGraph);
}
















int phaseOne(Graph *g, int *cliques, float minimum, float threshold){
    int changed = 1;
    int iters = 0;

    Graph * deviceGraph;
    copyGraphToDevice(g, &deviceGraph);

    float * deviceSigmaTot;
    HANDLE_ERROR(cudaMalloc((void**) &deviceSigmaTot, sizeof(float) * g->size));
    thrust::device_ptr<float> deviceSigmaTot_ptr(deviceSigmaTot);
    thrust::fill(deviceSigmaTot_ptr, deviceSigmaTot_ptr + g->size, (float) 0);


    int * deviceCliques;
    copyArrayToDevice(cliques, g->size, &deviceCliques);

    recalcSigmaTotPar<<<(g->size + 255)/256, 256>>>(deviceGraph, deviceSigmaTot, deviceCliques);

    int * deviceCliqueSizes;
    HANDLE_ERROR(cudaMalloc((void**) &deviceCliqueSizes, sizeof(int) * g->size));
    thrust::device_ptr<int> deviceCliqueSizes_ptr(deviceCliqueSizes);
    thrust::fill(deviceCliqueSizes_ptr, deviceCliqueSizes_ptr + g->size, (int) 0);

    int nMoves = g->size;

    int movesDone = 0;
    int * movesDoneDevice;
    HANDLE_ERROR(cudaMalloc((void**) &movesDoneDevice, sizeof(int)));
    HANDLE_ERROR(cudaMemcpy((void*) movesDoneDevice, (void*)&movesDone, sizeof(int), cudaMemcpyHostToDevice));


    int * deviceSizes;
    HANDLE_ERROR(cudaMalloc((void**) &deviceSizes, sizeof(int) * g->size));
    thrust::device_ptr<int> deviceSizes_ptr(deviceSizes);
    thrust::fill(deviceSizes_ptr, deviceSizes_ptr + g->size, (float) 0);


    calcNeighbours<<<(g->size + 255)/256, 256>>>(deviceGraph, deviceSizes);


    int maxNeighbours = thrust::reduce(deviceSizes_ptr, deviceSizes_ptr + g->size, (int) 0, thrust::maximum<int>());


    float m = thrust::reduce(deviceSigmaTot_ptr, deviceSigmaTot_ptr + g->size, (float) 0, thrust::plus<float>());


    m = m/2;
    if(DEBUG){
        printf("calculated:\n");
        printf("m=%f\n, maxN=%d", m, maxNeighbours);
    }


//    if(minimum < 1){
//        printf("exiting\n");
//        exit(10);
//
//    }

    float mod = modularity(g, cliques);

    if(DEBUG){
        printf("modularity: %f\n", mod);
    }


    while(changed != 0 ){

        Move empty = {.vertice=0,.toClique=0,.gain=0};


        Move * deviceMoves;
        HANDLE_ERROR(cudaMalloc((void**) &deviceMoves, sizeof(Move) * nMoves));
        thrust::device_ptr<Move> deviceMoves_ptr(deviceMoves);
        thrust::fill(deviceMoves_ptr, deviceMoves_ptr + nMoves, empty);

        HANDLE_ERROR(cudaMemcpy(deviceCliques, cliques, sizeof(int) * g->size, cudaMemcpyHostToDevice));

        calculateCliqueSizes<<<(g->size + 255)/256, 256>>>(deviceGraph, deviceCliques, deviceCliqueSizes);

        thrust::fill(deviceSigmaTot_ptr, deviceSigmaTot_ptr + g->size, (float) 0);
        recalcSigmaTotPar<<<(g->size + 255)/256, 256>>>(deviceGraph, deviceSigmaTot, deviceCliques);

        if(DEBUG){
            printf("---------------------------- small iter %d ------------------------------------------\n", iters);
        }
        changed = 0;
        iters++;
        movesDone = 0;
        HANDLE_ERROR(cudaMemcpy((void*) movesDoneDevice, (void*)&movesDone, sizeof(int), cudaMemcpyHostToDevice));

        //todo change that to max threads 256 :)
        calculateMoves<<<g->size, maxNeighbours, maxNeighbours * 2 * sizeof(float)>>>(deviceGraph, deviceCliques, deviceCliqueSizes, deviceMoves, m,deviceSigmaTot, minimum, movesDoneDevice);

        HANDLE_ERROR(cudaMemcpy((void*)&movesDone, (void*) movesDoneDevice, sizeof(int), cudaMemcpyDeviceToHost));
        movesDone = movesDone > 0 ?  movesDone -1 : 0;

        if(DEBUG){
            printf("calculated %d moves\n", movesDone);
        }


        if(movesDone > 0){
            changed = 1;
        }
        //sort moves //TODO
//        thrust::stable_sort(deviceMoves.begin(),deviceMoves.end(), compareMovesThrust);

        if(DEBUG){
            printf("moves sorted\n");
        }

        Move * moves = (Move*) malloc(nMoves * sizeof(Move));

        //wydobyć moves
        HANDLE_ERROR(cudaMemcpy( moves, deviceMoves, sizeof(Move)*nMoves, cudaMemcpyDeviceToHost));

        if(DEBUG){
            printf("moves moved to host\n");
        }

        int movesToApply = calculateMovesToApply(1, movesDone, nMoves);

        int* newCliques = (int*) malloc(sizeof(int) * g->size);
        memcpy(newCliques, cliques, sizeof(int) * g->size);
        float newMod = previewModularity(g, newCliques, moves, movesDone, movesToApply, 1);

        if(DEBUG){
            printf("moves:\n");

            for(int i=0; i<movesDone; i++){
                Move m = moves[i];
                printf("move %d from %d to %d \tgain:%f\n", m.vertice, cliques[m.vertice], m.toClique, m.gain);
            }


            printf("modularity gain if %d applied=%f\n",movesToApply, newMod - mod);
        }


        if(movesDone > 0){
            Move bestMove = moves[0];
            float bestdQ = bestMove.gain;
            int movesIter = 2;
            while(bestdQ > 0 && (newMod - mod < threshold) && (movesToApply > 1 || bestdQ > threshold)){
                movesToApply = calculateMovesToApply(movesIter, movesDone, nMoves);
                memcpy(newCliques, cliques, sizeof(int) * g->size);
                newMod = previewModularity(g, newCliques, moves, movesDone, movesToApply, 0);

                if(DEBUG){
                    printf("modularity gain if %d applied=%f\n",movesToApply, newMod - mod);
                }

                movesIter++;
            }
            if (newMod - mod > threshold) {
                memcpy(cliques, newCliques, sizeof(int) * g->size);
                mod = newMod;
//                printf("modularity: %f\n", modularity(g, cliques));
            }
            if(movesToApply == 1 && bestdQ < threshold){
                changed = 0;
            }
        } else {
            changed = 0;
        }
        cudaFree(deviceMoves);
        free(moves);
        free(newCliques);

    }

    cudaFree(deviceSigmaTot);
    cudaFree(deviceCliques);
    cudaFree(deviceCliqueSizes);
    cudaFree(movesDoneDevice);
    cudaFree(deviceSizes);

    destroyDeviceGraph(deviceGraph);

    //distroy all ptrs and graph
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
    printf("# %s  -f mtx-matrix-file -g min-gain [-v]\n", name);
    printf("#     mtx-matrix-file   matrix in mtx format representing undirected weighted graph\n");
    printf("#     min-gain    minimal modularity gain to move a node between communities\n");
    printf("#     -v    verbose mode, printing communities\n");
}


int main(int argc, char **argv){
    char * fileName;
    int verbose = 0;
    float min_gain = 0;
    if(argc != 5 && argc != 6){
        printf("wrong number of arguments!\n");
        printUsage(argv[0]);
        return 1;
    } else if(strcmp(argv[1], "-f") && strcmp(argv[3], "-g")){
        if(argc == 6){
            if(strcmp(argv[5], "-v") != 0){
                printf("what is that gibberish?!\n");
                printUsage(argv[0]);
                return 1;
            }
            verbose = 1;
        }
        fileName = argv[2];
        min_gain = strtof(argv[4], NULL);
    } else {
        printf("what is that gibberish?!\n");
        printUsage(argv[0]);
        return 1;
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

    float threshold = min_gain;

    // profiler at hangGlider_4 th=0.00001f

    float mod = modularity(g, cliques);
//    printf("modularity:%f\n", mod);

    cudaEvent_t start, stop;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start, 0));


    int iter = 10;
    while(iter > 1 || minimum > threshold/10.f){

        if(DEBUG){
            printf("========= PHASE 1 ==================\n");
        }
        minimum = 0.1 / (2 + bigLoopIteration) - 0.02;
        minimum = minimum < threshold/20.f ? threshold/20.f : minimum;
//        printf("min:%f\n", minimum);
        iter = phaseOne(g, cliques, minimum, threshold);
//        printCliques(g->size, cliques);
        if(DEBUG){
            printf("========= PHASE 2 ==================\n");
        }
        if(iter > 1){
            phaseTwo(g, cliques);
        }
//        printEdges(g);
        updateOldCliques(g, cliques);
//        printf("modularity:%f\n", modularity(g, cliques));
//        printCliques(g->size, cliques);
        bigLoopIteration += 1;
    }
    if(DEBUG){
        printf("converged after %d iterations!\n", bigLoopIteration+1);
    }
//    if(verbose != 0){
//        printCliques(g->size, cliques);
//    }


    HANDLE_ERROR(cudaEventRecord(stop, 0));
    HANDLE_ERROR(cudaEventSynchronize(stop));

    float elapsedTime;
    HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));

    printf("%f\n", modularity(g, cliques));

    printf("%f %f\n", elapsedTime, elapsedTime);

    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));

    free(cliques);

    destroyGraph(g);



    return 0;
}