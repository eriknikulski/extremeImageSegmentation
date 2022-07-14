//
// Created by Erik Nikulski on 10.04.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_VEC_H
#define EXTREMEIMAGESEGMENTATION_VEC_H

typedef struct Vec {
    double x;
    double y;
    double z;
} Vec;

void printVec(Vec* v);

void printVecs(Vec* v, int dim);

void addVecs(Vec* v1, Vec* v2, Vec* res);

void subtractVec(Vec* v1, Vec* v2, Vec* res);

Vec* getSubVec(Vec* v1, Vec* v2);

Vec* addVecLists(Vec* v1, Vec* v2, int dim);

void multVecs(double t, Vec* vec, Vec* res);

Vec* listMultVec(double* t, Vec* vec, int dim);

Vec* listMultVecs(double* t, Vec* vecs, int dim);

Vec* getRandVec();

Vec* getVecLogNearCube(int base);

Vec* getRandVecOnCube();

double getDistSq(Vec* v1, Vec* v2);

double getDist(Vec* v1, Vec* v2);

double getLength(Vec* v);

Vec* calcCrossProduct(Vec* v1, Vec* v2, Vec* res);

Vec* getCrossProduct(Vec* v1, Vec* v2);

double getDotProd(Vec* v1, Vec* v2);

int negEqualVecs(Vec* v1, Vec* v2);

int equalVecs(Vec* v1, Vec* v2);

int isZero(double d);

int isAlmostZero(double d);

void discretizeVec(Vec* v, int size);

void discretizeVecs(Vec* vecs, int dim, int size);

int isOutsideUnitCube(Vec* vec);

#endif //EXTREMEIMAGESEGMENTATION_VEC_H
