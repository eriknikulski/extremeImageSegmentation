//
// Created by Erik Nikulski on 10.04.22.
//

#include "vec.h"

#include "params.h"
#include "utility.h"

#include "float.h"
#include "math.h"
#include "memory.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

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

void subtractVec(Vec* v1, Vec* v2, Vec* res) {
    res->x = v1->x - v2->x;
    res->y = v1->y - v2->y;
    res->z = v1->z - v2->z;
}

Vec* getSubVec(Vec* v1, Vec* v2) {
    Vec* res = malloc(sizeof(Vec));
    subtractVec(v1, v2, res);
    return res;
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

Vec* getRandVec() {
    Vec* vec = malloc(sizeof(Vec));

    vec->x = getRandD();
    vec->y = getRandD();
    vec->z = getRandD();

    return vec;
}

Vec* getVecLogNearCube(int base) {
    Vec* vec = malloc(sizeof(Vec));

    vec->x = 1.0 + log(rand() + 1) / log(base);
    vec->y = 1.0 + log(rand() + 1) / log(base);
    vec->z = 1.0 + log(rand() + 1) / log(base);

    return vec;
}

Vec* getRandVecOnCube() {
    Vec* vec = malloc(sizeof(Vec));
    int fix =  rand() % 6;

    vec->x = getRandD();
    vec->y = getRandD();
    vec->z = getRandD();

    if (fix == 0) {
        vec->x = 0;
    } else if(fix == 1) {
        vec->x = 1;
    } else if(fix == 2) {
        vec->y = 0;
    } else if(fix == 3) {
        vec->y = 1;
    } else if(fix == 4) {
        vec->z = 0;
    } else if(fix == 5) {
        vec->z = 1;
    }
    return vec;
}

double getDistSq(Vec* v1, Vec* v2) {
    return pow(v2->x - v1->x, 2) + pow(v2->y - v1->y, 2) + pow(v2->z - v1->z, 2);
}

double getDist(Vec* v1, Vec* v2) {
    return sqrt(pow(v2->x - v1->x, 2) + pow(v2->y - v1->y, 2) + pow(v2->z - v1->z, 2));
}

double getLength(Vec* v) {
    return sqrt(pow(v->x, 2) + pow(v->y, 2) + pow(v->z, 2));
}

Vec* calcCrossProduct(Vec* v1, Vec* v2, Vec* res) {
    res->x = v1->y * v2->z - v1->z * v2->y;
    res->y = v1->z * v2->x - v1->x * v2->z;
    res->z = v1->x * v2->y - v1->y * v2->x;
    return res;
}

Vec* getCrossProduct(Vec* v1, Vec* v2) {
    Vec* cross = malloc(sizeof(Vec));
    cross->x = v1->y * v2->z - v1->z * v2->y;
    cross->y = v1->z * v2->x - v1->x * v2->z;
    cross->z = v1->x * v2->y - v1->y * v2->x;
    return cross;
}

double getDotProd(Vec* v1, Vec* v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

int negEqualVecs(Vec* v1, Vec* v2) {
    if (!v1 && !v2) return 1;
    if (!v1 || !v2) return 0;
    if (fabs(v1->x + v2->x) < DBL_EPSILON && fabs(v1->y + v2->y) < DBL_EPSILON && fabs(v1->z + v2->z) < DBL_EPSILON) return 1;
    return 0;
}

int equalVecs(Vec* v1, Vec* v2) {
    if (!v1 && !v2) return 1;
    if (!v1 || !v2) return 0;
    if (fabs(v1->x - v2->x) < DBL_EPSILON && fabs(v1->y - v2->y) < DBL_EPSILON && fabs(v1->z - v2->z) < DBL_EPSILON) return 1;
    return 0;
}

int isZero(double d) {
    return fabs(d) < DBL_EPSILON;
}

int isAlmostZero(double d) {
    return fabs(d) < IDC_SMALL;
}

void discretizeVec(Vec* v, int size) {
    v->x = round(v->x * (double)(size - 1));
    v->y = round(v->y * (double)(size - 1));
    v->z = round(v->z * (double)(size - 1));
}

void discretizeVecs(Vec* vecs, int count, int size) {
    for (int i = 0; i < count; ++i) {
        discretizeVec(&vecs[i], size);
    }
}

int isOutsideUnitCube(Vec* vec) {
    if (vec->x < 0 || vec->x > 1 || vec->y < 0 || vec->y > 1 || vec->z < 0 || vec->z > 1)
        return 1;
    return 0;
}
