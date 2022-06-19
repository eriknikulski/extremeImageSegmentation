//
// Created by Erik Nikulski on 16.05.22.
//

#include "voronoi.h"

#ifndef EXTREMEIMAGESEGMENTATION_GRAPH_H
#define EXTREMEIMAGESEGMENTATION_GRAPH_H

typedef struct gNode {
    unsigned int id;
    int* neighbors;
    unsigned int neighbourCount;
} gNode;

typedef struct Graph {
    gNode* nodes;
    unsigned int count;
} Graph;

int nonNegCount(int* l, unsigned int n);

Graph* toGraph(Cell** cells, int n);

int* growSeeds(Graph* graph, int* seeds, int n);

#endif //EXTREMEIMAGESEGMENTATION_GRAPH_H
