//
// Created by Erik Nikulski on 07.04.22.
//

#include "spline.h"

#include "vec.h"

#include "assert.h"
#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"


void printSplines(Vec** splines, int n, int dim) {
    for (int i = 0; i < n; ++i) {
        printf("Spline %d:\n", i);
        printVecs(splines[i], dim);
        printf("\n");
    }
}

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

Vec* getCatmullRomSpline(Vec* p0, Vec* p1, Vec* p2, Vec* p3, double alpha, int nPoints) {
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

    free(A1);
    free(A2);
    free(A3);
    free(B1);
    free(B2);
    return C;
}

double getMinDist(Vec* s1, Vec* s2, int dim) {
    double minDist = DBL_MAX;
    double dist;
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            dist = getDist(&s1[i], &s2[j]);
            minDist = dist < minDist ? dist : minDist;
        }
    }
    return minDist;
}

static int isClose(Vec** splines, Vec* spline, int splinesN, double minDist, int dim) {
    for (int i = 0; i < splinesN; ++i) {
        if (getMinDist(splines[i], spline, dim) < minDist) return 1;
    }
    return 0;
}

Vec** getNSplines(SplineParams* splineParams) {
    assert(splineParams->minDist < 1 && splineParams->minDist > 0);
    Vec** splines = malloc(sizeof(Vec) * splineParams->nPoints * splineParams->nSplines);
    Vec* spline;
    int i = 0;
    int base = 1000;

    Vec* c0;
    Vec* c1;
    Vec* c2;
    Vec* c3;

    while (i < splineParams->nSplines) {
        c0 = getVecLogNearCube(base);
        c1 = getRandVecOnCube();
        c2 = getRandVecOnCube();
        c3 = getVecLogNearCube(base);

        spline = getCatmullRomSpline(c0, c1, c2, c3, splineParams->alpha, splineParams->nPoints);
        if (isClose(splines, spline, i, splineParams->minDist, splineParams->nPoints) ||
                isOutsideUnitCube(&spline[splineParams->nPoints / 4]) ||
                isOutsideUnitCube(&spline[splineParams->nPoints / 2]) ||
                isOutsideUnitCube(&spline[splineParams->nPoints * 3 / 4])) {
            free(spline);
            continue;
        }
        printVec(c0);
        printVec(c1);
        printVec(c2);
        printVec(c3);
        printf("\n\n");
        splines[i] = spline;
        ++i;
    }
    return splines;
}

void discretizeSplines(Vec** splines, int n, int dim, int size) {
    for (int i = 0; i < n; ++i) {
        discretizeVecs(splines[i], dim, size);
    }
}

double splineDist(Vec* point, Vec* spline, int dim) {
    double sDist = DBL_MAX;
    double cDist;
    for (int i = 0; i < dim; ++i) {
        cDist = getDist(point, &spline[i]);
        if (cDist < sDist) sDist = cDist;
        if (isZero(sDist)) return 0.0;
    }
    return sDist;
}

double splinesDist(Vec* point, Vec** splines, SplineParams* splineParams, int* spline) {
    double sDist = DBL_MAX;
    double cDist;
    *spline = 0;
    for (int i = 0; i < splineParams->nSplines; ++i) {
        cDist = splineDist(point, splines[i], splineParams->nPoints);
        if (cDist < sDist) {
            sDist = cDist;
            *spline = i;
        }
        if (isZero(sDist)) return 0.0;
    }
    return sDist;
}
