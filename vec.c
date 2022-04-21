//
// Created by Erik Nikulski on 10.04.22.
//

#include "vec.h"

#include "utility.h"

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

    vec->x = getRand();
    vec->y = getRand();
    vec->z = getRand();

    return vec;
}

Vec* getRandVecOnCube() {
    Vec* vec = malloc(sizeof(Vec));
    int fix =  rand() % 6;

    vec->x = getRand();
    vec->y = getRand();
    vec->z = getRand();

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

double getDist(Vec* v1, Vec* v2) {
    return sqrt(pow(v2->x - v1->x, 2) + pow(v2->y - v1->y, 2) + pow(v2->z - v1->z, 2));
}

double getLength(Vec* v) {
    return sqrt(pow(v->x, 2) + pow(v->y, 2) + pow(v->z, 2));
}

double getDistLine(Vec* lv1, Vec* lv2, Vec* v) {
    Vec dir = {.x = lv2->x - lv1->x,
               .y = lv2->y - lv1->y,
               .z = lv2->z - lv1->z};
    Vec tmp = {.x = lv1->x - v->x,
               .y = lv1->y - v->y,
               .z = lv1->z - v->z};
    Vec cross = {.x = tmp.y * dir.z - tmp.z * dir.y,
                 .y = tmp.x * dir.z - tmp.z * dir.x,
                 .z = tmp.x * dir.y - tmp.y * dir.x};
    return getLength(&cross) / getLength(&dir);
}
