//
// Created by Erik Nikulski on 07.04.22.
//

#include "spline.h"

#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "vec.h"


double* listDiv(double* l, double d, int dim) {
    double* res = malloc(sizeof(double) * dim);
    for (int i = 0; i < dim; ++i) {
        res[i] = l[i] / d;
    }
    return res;
}

double* listSubR(double min, double* l, int dim) {
    double* res = malloc(sizeof(double) * dim);
    for (int i = 0; i < dim; ++i) {
        res[i] = min - l[i];
    }
    return res;
}

double* listSubL(double* l, double min, int dim) {
    double* res = malloc(sizeof(double) * dim);
    for (int i = 0; i < dim; ++i) {
        res[i] = l[i] - min;
    }
    return res;
}

double* linespace(double start, double stop, int num) {
    assert(start < stop);
    assert(num > 0);

    double range = stop - start;
    double delta = range / (num -1);
    double* numbers = malloc(sizeof(double) * num);

    numbers[0] = start;
    for (int i = 1; i < num; ++i) {
        numbers[i] = numbers[i - 1] + delta;
    }
    return numbers;
}

double tj(double ti, Vec* pi, Vec* pj, double alpha) {
    return pow(pow(pj->x - pi->x, 2) + pow(pj->y - pi->y, 2) + pow(pj->z - pi->z, 2), alpha/2) + ti;
}

Vec* calcLineSeq(double t0, double t1, double* t, Vec* p0, Vec* p1, int dim) {
    return addVecLists(listMultVec(listDiv(listSubR(t1, t, dim), t1 - t0, dim), p0, dim),
                       listMultVec(listDiv(listSubL(t, t0, dim), t1 - t0, dim), p1, dim),
                       dim);
}

Vec* calcLineSeqVecList(double t0, double t1, double* t, Vec* vecs1, Vec* vecs2, int dim) {
    return addVecLists(listMultVecs(listDiv(listSubR(t1, t, dim), t1 - t0, dim), vecs1, dim),
                       listMultVecs(listDiv(listSubL(t, t0, dim), t1 - t0, dim), vecs2, dim),
                       dim);
}

Vec* getCatmullRomSpline(Vec* p0, Vec* p1, Vec* p2, Vec* p3, int nPoints) {
    double alpha = 0.5;
    int dim = nPoints;

    double t0 = 0;
    double t1 = tj(t0, p0, p1, alpha);
    double t2 = tj(t1, p1, p2, alpha);
    double t3 = tj(t2, p2, p3, alpha);

    double* t = linespace(t1, t2, dim);

    Vec* A1 = calcLineSeq(t0, t1, t, p0, p1, dim);
    Vec* A2 = calcLineSeq(t1, t2, t, p1, p2, dim);
    Vec* A3 = calcLineSeq(t2, t3, t, p2, p3, dim);

    Vec* B1 = calcLineSeqVecList(t0, t2, t, A1, A2, dim);
    Vec* B2 = calcLineSeqVecList(t1, t3, t, A2, A3, dim);

    Vec* C = calcLineSeqVecList(t1, t2, t, B1, B2, dim);
    return C;
}
