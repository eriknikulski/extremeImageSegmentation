//
// Created by Erik Nikulski on 11.04.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_IMAGE_H
#define EXTREMEIMAGESEGMENTATION_IMAGE_H

#include "vec.h"
#include "voronoi.h"

#include "stdint.h"

typedef struct Pixel {
    Vec* v;
    uint8_t value;
    double dist;
    struct Seed* grouping;
    double delta;
    int inSSL;
    Vec* particle;
} Pixel;

typedef struct Bitmap {
    Pixel* pixels;
    int size;
    int size2;
} Bitmap;

typedef struct ImageParams {
    int imageSize;
    double theta_0;
    double theta_1;
    double sigma_b;
    double mu_b;
    double sigma_c;
    double mu_c;
} ImageParams;

typedef struct VoronoiParams {
    int nInitialCells;
    int nCells;
    int srgPrecision;
} VoronoiParams;

typedef struct SplineParams {
    double alpha;
    double minDist;
    int nPoints;
    int nSplines;
} SplineParams;

Pixel* getPixel(Bitmap* bitmap, int x, int y, int z);

Pixel** getNeighbors(Bitmap* bitmap, Pixel* p, int* count);

void writeBitmap(Bitmap* bitmap, char* fname);

void createSplineImage(Vec** splines, int nSplines, int dim, ImageParams* imageParams, char* fname);

void createVoronoiImage(Cell** cells, int nCells, ImageParams* imageParams, char* fname);

#endif //EXTREMEIMAGESEGMENTATION_IMAGE_H
