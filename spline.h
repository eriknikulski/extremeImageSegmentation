//
// Created by Erik Nikulski on 07.04.22.
//

#include "params.h"
#include "vec.h"

#ifndef EXTREMEIMAGESEGMENTATION_SPLINE_H
#define EXTREMEIMAGESEGMENTATION_SPLINE_H


void printSplines(Vec** splines, int n, int dim);

double* listDiv(double* l, double d, int dim);

double* listSubR(double min, double* l, int dim);

double* listSubL(double* l, double min, int dim);

double* linespace(double start, double stop, int num);

double tj(double ti, Vec* pi, Vec* pj, double alpha);

Vec* calcLineSeq(double t0, double t1, double* t, Vec* p0, Vec* p1, int dim);

Vec* getCatmullRomSpline(Vec* p0, Vec* p1, Vec* p2, Vec* p3, double alpha, int nPoints);

double getMinDist(Vec* s1, Vec* s2, int dim);

Vec** getNSplines(SplineParams* splineParams);

void discretizeSpline(Vec* spline, int dim, int size);

void discretizeSplines(Vec** splines, int n, int dim, int size);

double splineDist(Vec* point, Vec* spline, int dim);

double splinesDist(Vec* point, Vec** splines, SplineParams* splineParams, int* spline);

Vec* getSplineSeeds(Vec** splines, SplineParams* splineParams);

#endif //EXTREMEIMAGESEGMENTATION_SPLINE_H
