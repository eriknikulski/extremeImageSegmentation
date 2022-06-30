//
// Created by Erik Nikulski on 19.04.22.
//

#include "vec.h"

#ifndef EXTREMEIMAGESEGMENTATION_VORONOI_H
#define EXTREMEIMAGESEGMENTATION_VORONOI_H

typedef Vec Node;

typedef struct FaceCalc {
    Node* n1;
    Node* n2;
    Node* n3;

    Vec* u;
    Vec* v;
    Vec* n;

    double nn;
} FaceCalc;

typedef struct Face {
    Node* nodes;
    unsigned int count;
    Node* normalVec;
    FaceCalc* faceCalcs;
    int nFaceCalcs;
} Face;

typedef struct Cell {
    int id;
    Node* nodes;
    Face* faces;
    int* neighbors;
    unsigned int nodeCount;
    unsigned int faceCount;
    unsigned int neighborCount;
} Cell;

Vec* getRandParticles(int n);

void createParticleFile(int n, char* fname);

void createCells(int n, char* fname);

void setParticleId(Cell* cell, char* s);

void setNodeCountStr(Cell* c, char* s);

void setVerticesStr(Cell* c, char* s);

int getNodeCountStr(char* s);

void setFaceCountStr(Cell* c, char* s);

void setFaceNodesStr(Cell* c, Face* f, char* s);

void setFacesStr(Cell* c, char* s);

void setNormalVecFace(Cell* c, char* s);

void setNeighboringCells(Cell* cell, char* s);

Cell** readCells(int n, char* fname);

Cell** getCells(int n, char* fname);

int isNodeInsideFace(Face* f, Node* n);

void removeDupNodes(Cell* c);

void removeDupFaces(Cell* c);

Cell** actualizeAssignment(Cell** cs, int* seeds, int* assign, int nOld, int nNew);

Cell** mergeCells(Cell** cells, int nOld, int nNew);

void discretizeCell(Cell* cell, int size);

void discretizeCells(Cell** cells, int n, int size);

double getDistLineSeg(Vec* v1, Vec* v2, Vec* v);

double getDistFace(Face* f, Vec* v);

double voronoiDist(Vec* v, Cell** cells, int nCells);

int equalFaces(Face* f1, Face* f2);

void printCells(Cell** cells, int count);

#endif //EXTREMEIMAGESEGMENTATION_VORONOI_H
