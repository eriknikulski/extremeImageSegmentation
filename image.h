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
    int particleId;
    Face* closestFace;
    int cellId;
} Pixel;

typedef struct Bitmap {
    Pixel* pixels;
    int size;
    int size2;
    int reverse;
} Bitmap;

Pixel* getPixel(Bitmap* bitmap, int x, int y, int z);

Pixel** getNeighbors(Bitmap* bitmap, Pixel* p, int* count);

Vec** getSeeds(Bitmap* bitmap, uint8_t value, int* count);

Vec** getSeedsWithBlockRad(Bitmap *bitmap, int value, int *count, double radius);

Vec** addBorderSeed(Vec** seeds, int *count);

void setNoiseValuesBitmap(Bitmap* bitmap, ImageParams* imageParams);

void setValuesBitmap(Bitmap* bitmap, ImageParams* imageParams);

void freeBitmap(Bitmap *bitmap);

Bitmap* copyBitmap(Bitmap *orig);

Bitmap* initializeBitmap(ImageParams* imageParams);

void writeBitmap(Bitmap* bitmap, char* fname);

Bitmap* calcVoronoiDist(Bitmap* bitmap, Cell** cells, int nCells);

Bitmap* removeDistMissingFaces(Bitmap* bitmap, Cell** cells, int nCells);

Bitmap* calcMissingDist(Bitmap* bitmap, Cell** cells, int nCells);

Bitmap* createSplineImage(Vec** splines, SplineParams* splineParams, ImageParams* imageParams);

Bitmap* createVoronoiImage(Cell** cells, VoronoiParams* voronoiParams, ImageParams* imageParams);

void printBitmapPixelDist(Bitmap* bitmap);

void setParticleIds(Bitmap* bitmap, int n);

double getRandsIndex(Bitmap* orig, Bitmap* srg, int nParticles, int nSeeds);

double getVariationOfInformation(Bitmap* orig, Bitmap* srg, int nParticles, int nSeeds, double *falseJoins, double *falseCuts);


#endif //EXTREMEIMAGESEGMENTATION_IMAGE_H
