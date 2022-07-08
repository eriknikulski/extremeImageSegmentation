//
// Created by Erik Nikulski on 11.04.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_IMAGE_H
#define EXTREMEIMAGESEGMENTATION_IMAGE_H

#include "params.h"
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

Pixel* getPixel(Bitmap* bitmap, int x, int y, int z);

Pixel** getNeighbors(Bitmap* bitmap, Pixel* p, int* count);

void writeBitmap(Bitmap* bitmap, char* fname);

void setValuesBitmap(Bitmap* bitmap, ImageParams* imageParams);

Bitmap* createSplineImage(Vec** splines, SplineParams* splineParams, ImageParams* imageParams);

Bitmap* createVoronoiImage(Cell** cells, VoronoiParams* voronoiParams, ImageParams* imageParams);

double getRandsIndex(Bitmap* orig, Bitmap* srg);

double getVariationOfInformation(Bitmap* orig, Bitmap* srg, Vec* particles, int n);


#endif //EXTREMEIMAGESEGMENTATION_IMAGE_H
