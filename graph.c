//
// Created by Erik Nikulski on 16.05.22.
//
#include "graph.h"

#include "utility.h"
#include "voronoi.h"

#include "stdlib.h"


int nonNegCount(int* l, unsigned int n) {
    int res = 0;
    for (int i = 0; i < n; ++i) {
        if (l[i] >= 0) ++res;
    }
    return res;
}


Graph* toGraph(Cell** cells, int n) {
    Cell* cell;
    gNode* node;
    Graph* graph = malloc(sizeof(Graph));
    graph->count = n;
    graph->nodes = malloc(sizeof(Node) * n);

    int c;

    for (int i = 0; i < n; ++i) {
        c = 0;
        cell = cells[i];
        node = &graph->nodes[i];

        node->id = cell->id;
        node->neighbourCount = nonNegCount(cell->neighbors, cell->faceCount);

        node->neighbors = malloc(sizeof(int) * cell->faceCount);

        for (int j = 0; j < cell->faceCount; ++j) {
            if (cell->neighbors[j] >= 0) {
                node->neighbors[c] = cell->neighbors[j];
                ++c;
            }
        }
    }

    return graph;
}

int* growSeeds(Graph* graph, int* seeds, int n) {
    int* awail = malloc(sizeof(int) * graph->count);
    int count = (int)graph->count - n;
    int* assign = malloc(sizeof(int) * graph->count);
    int* neigh = malloc(sizeof(int) * (int)graph->count);

    for (int i = 0; i < graph->count; ++i) {
        awail[i] = 1;
    }

    for (int i = 0; i < n; ++i) {
        awail[seeds[i]] = 0;
        assign[seeds[i]] = seeds[i];
    }

    while (count > 0) {
        int c = 0;
        int randSeed = seeds[getRand(0, n - 1)];

        for (int i = 0; i < (int)graph->count; ++i) {
            if (assign[i] != randSeed) continue;
            for (int j = 0; j < (int)graph->nodes[i].neighbourCount; ++j) {
                if (awail[graph->nodes[i].neighbors[j]]) neigh[c++] = graph->nodes[i].neighbors[j];
            }
        }

        if (c == 0) continue;

        int randNeigh = getRandFrom(neigh, c);
        assign[randNeigh] = randSeed;
        awail[randNeigh] = 0;
        --count;
    }
    free(neigh);
    return assign;
}
