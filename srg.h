//
// Created by Erik Nikulski on 01.07.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_SRG_H
#define EXTREMEIMAGESEGMENTATION_SRG_H

#include "vec.h"
#include "image.h"

typedef struct SSL {
    struct SSL* next;
    Pixel* p;
    double delta;
} SSL;

typedef struct Seed {
    int id;
    double sum;
    int count;
    Pixel* p;
} Seed;

Bitmap* srg(char* path, Vec* vSeeds, int nSeeds);

void setGroupings(Bitmap* bitmap);

#endif //EXTREMEIMAGESEGMENTATION_SRG_H
