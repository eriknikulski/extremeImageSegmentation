//
// Created by Erik Nikulski on 07.04.22.
//

#include "vec.h"

#ifndef EXTREMEIMAGESEGMENTATION_SPLINE_H
#define EXTREMEIMAGESEGMENTATION_SPLINE_H


double* listDiv(double* l, double d, int dim);

double* listSubR(double min, double* l, int dim);

double* listSubL(double* l, double min, int dim);

double* linespace(double start, double stop, int num);

double tj(double ti, Vec* pi, Vec* pj, double alpha);

Vec* calcLineSeq(double t0, double t1, double* t, Vec* p0, Vec* p1, int dim);

Vec* getCatmullRomSpline(Vec* p0, Vec* p1, Vec* p2, Vec* p3, int nPoints);

#endif //EXTREMEIMAGESEGMENTATION_SPLINE_H
