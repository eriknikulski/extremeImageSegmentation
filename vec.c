//
// Created by Erik Nikulski on 10.04.22.
//

#include <memory.h>
#include "vec.h"

#include "stdlib.h"
#include "stdio.h"

void printVec(Vec* v) {
    printf("Vec { .x = %lf , .y = %lf, .z = %lf }\n", v->x, v->y, v->z);
}

void printVecs(Vec* v, int dim) {
    for (int i = 0; i < dim; ++i) {
        printf("%d: ", i);
        printVec(&v[i]);
    }
}

void addVecs(Vec* v1, Vec* v2, Vec* res) {
    res->x = v1->x + v2->x;
    res->y = v1->y + v2->y;
    res->z = v1->z + v2->z;
}

void* malloc0(size_t s) {
    void* p = malloc(s);
    memset(p, '\0', s);
    return p;
}

Vec* addVecLists(Vec* l1, Vec* l2, int dim) {
    Vec* res = malloc0(sizeof(Vec) * dim);
    for (int i = 0; i < dim; ++i) {
        addVecs(&l1[i], &l2[i], &res[i]);
    }
    return res;
}

void multVecs(double t, Vec* vec, Vec* res) {
    res->x = t * vec->x;
    res->y = t * vec->y;
    res->z = t * vec->z;
}

Vec* listMultVec(double* t, Vec* vec, int dim) {
    Vec* res = malloc0(sizeof(Vec) * dim);
    for (int i = 0; i < dim; ++i) {
        multVecs(t[i], vec, &res[i]);
    }
    return res;
}

Vec* listMultVecs(double* t, Vec* vecs, int dim) {
    Vec* res = malloc0(sizeof(Vec) * dim);
    for (int i = 0; i < dim; ++i) {
        multVecs(t[i], &vecs[i], &res[i]);
    }
    return res;
}
