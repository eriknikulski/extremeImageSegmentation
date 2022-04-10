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

Vec* addVecLists(Vec* v1, Vec* v2, int dim);

void multVecs(double t, Vec* vec, Vec* res);

Vec* listMultVec(double* t, Vec* vec, int dim);

Vec* listMultVecs(double* t, Vec* vecs, int dim);

#endif //EXTREMEIMAGESEGMENTATION_VEC_H
