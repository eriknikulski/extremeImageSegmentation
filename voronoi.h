//
// Created by Erik Nikulski on 19.04.22.
//

#include "vec.h"

#ifndef EXTREMEIMAGESEGMENTATION_VORONOI_H
#define EXTREMEIMAGESEGMENTATION_VORONOI_H

typedef Vec Node;

typedef struct Cell {
    Node* nodes;
    unsigned int count;
} Cell;

Vec* getRandParticles(int n);

void createParticleFile(int n, char* fname);

void createCells(int n, char* fname);

Cell** readCells(int n, char* fname);

Cell** getCells(int n, char* fname);

void discretizeCell(Cell* cell, int size);

void discretizeCells(Cell** cells, int n, int size);

double voronoiDist(Vec* v, Cell** cells, int nCells);

#endif //EXTREMEIMAGESEGMENTATION_VORONOI_H
