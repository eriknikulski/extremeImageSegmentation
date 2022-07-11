//
// Created by Erik Nikulski on 01.07.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_SRG_H
#define EXTREMEIMAGESEGMENTATION_SRG_H

#include "vec.h"
#include "image.h"

typedef struct SSL {
    struct SSL* next;
    struct SSL* nextHigher;
    Pixel* p;
    double delta;
} SSL;

typedef struct Seed {
    int id;
    double sum;
    int count;
    Pixel* p;
} Seed;

Bitmap* srg(Bitmap* bitmap, Vec** vSeeds, int nSeeds, int precision, char* imagePath, ImageParams* imageParams);

void setGroupings(Bitmap* bitmap);

#endif //EXTREMEIMAGESEGMENTATION_SRG_H
