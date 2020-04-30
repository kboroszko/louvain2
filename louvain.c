//
// Created by kajetan on 30.04.2020.
//

#include "louvain.h"

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