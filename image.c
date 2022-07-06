//
// Created by Erik Nikulski on 11.04.22.
//

#include "image.h"

#include "utility.h"
#include "spline.h"
#include "voronoi.h"
#include "png.h"

#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


/* Given "bitmap", this returns the pixel of bitmap at the point
   ("x", "y", "z"). */

Pixel* getPixel(Bitmap* bitmap, int x, int y, int z) {
    return bitmap->pixels + bitmap->size2 * z + bitmap->size * y + x;
}

static double getPixelValue(double d, ImageParams* imageParams) {
    double p = (double)1.0 / (double)((double)1 + exp(-imageParams->theta_0 - imageParams->theta_1 * d));

    // TODO: N1 and N2 in task are diff for splines and voronoi. Does it matter??
    if (getRandD() > p) {
        // select from N2
        return normalCDF(getRandD(), imageParams->sigma_b, imageParams->mu_b);
    } else {
        // select from N1
        return normalCDF(getRandD(), imageParams->sigma_c, imageParams->mu_c);
    }
}

static uint8_t to8Bit(double value) {
    if (value < 0) value = 0;
    if (value > 1) value = 1;

    return (uint8_t) (value * 255);
}

void writeBitmap(Bitmap* bitmap, char* fname) {
    size_t len;
    char* activeName;

    for (int z = 0; z < bitmap->size; ++z) {
        len = snprintf(NULL, 0, "%simg_%d.png", fname, z) + 1;
        activeName = malloc(len);
        snprintf(activeName, len, "%simg_%d.png", fname, z);

        if (save_bitmap_slice_to_file(bitmap, z, activeName)) {
            fprintf(stderr, "Error writing file %d.\n", z);
        }
        free(activeName);
    }

}

Bitmap* calcSplineDistance(Vec** splines, int nSplines, int dim, ImageParams* imageParams) {
    Vec p;
    double dist;
    double value;

    Bitmap* bitmap = malloc(sizeof(Bitmap));
    bitmap->size = imageParams->imageSize;
    bitmap->size2 = imageParams->imageSize * imageParams->imageSize;
    bitmap->pixels = malloc(sizeof(Pixel) * imageParams->imageSize * imageParams->imageSize * imageParams->imageSize);

    for (int z = 0; z < imageParams->imageSize; ++z) {
        for (int y = 0; y < imageParams->imageSize; ++y) {
            for (int x = 0; x < imageParams->imageSize; ++x) {
                p.x = (double) x;
                p.y = (double) y;
                p.z = (double) z;

                dist = splinesDist(&p, splines, nSplines, dim);
                value = getPixelValue(dist, imageParams);
                Pixel* pixel = getPixel(bitmap, x, y, z);
                pixel->value = to8Bit(value);
            }
        }
        printf("%d\n", z);
    }
    return bitmap;
}

void createSplineImage(Vec** splines, int nSplines, int dim, ImageParams* imageParams, char* fname) {
    discretizeSplines(splines, nSplines, dim, imageParams->imageSize);
    Bitmap* bitmap = calcSplineDistance(splines, nSplines, dim, imageParams);
    writeBitmap(bitmap, fname);
}

Pixel** getNeighbors(Bitmap* bitmap, Pixel* p, int* count) {
    int iter = 0;
    *count = 26;
    Pixel** neighbors = malloc(sizeof(Pixel*) * (*count));

    int x_min = -1;
    int x_max = 1;
    int y_min = -1;
    int y_max = 1;
    int z_min = -1;
    int z_max = 1;

    if (p->v->x == 0)
        x_min = 0;
    if (p->v->x == bitmap->size - 1)
        x_max = 0;
    if (p->v->y == 0)
        y_min = 0;
    if (p->v->y == bitmap->size - 1)
        y_max = 0;
    if (p->v->z == 0)
        z_min = 0;
    if (p->v->z == bitmap->size - 1)
        z_max = 0;

    for (int z = z_min; z <= z_max; ++z) {
        for (int y = y_min; y <= y_max; ++y) {
            for (int x = x_min; x <= x_max; ++x) {
                if (z == 0 && y == 0 && x == 0) continue;
                neighbors[iter] = getPixel(bitmap, (int)p->v->x + x, (int)p->v->y + y, (int)p->v->z + z);
                ++iter;
            }
        }
    }
    *count = iter;
    neighbors = realloc(neighbors, sizeof(Pixel*) * (*count));
    return neighbors;
}

Vec* getClosestParticle(Vec* v, Cell** cells, int nCells) {
    double sDist = DBL_MAX;
    double cDist;
    Vec* closest = NULL;
    for (int i = 0; i < nCells; ++i) {
        cDist = getDist(cells[i]->particle, v);
        if (cDist < sDist) {
            sDist = cDist;
            closest = cells[i]->particle;
        }
    }
    return closest;
}

Bitmap* initializeBitmap(Cell** cells, int nCells, ImageParams* imageParams) {
    Pixel* pixel;
    Bitmap* bitmap = malloc(sizeof(Bitmap));
    bitmap->size = imageParams->imageSize;
    bitmap->size2 = imageParams->imageSize * imageParams->imageSize;
    bitmap->pixels = malloc(sizeof(Pixel) * imageParams->imageSize * imageParams->imageSize * imageParams->imageSize);

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                pixel->v = malloc(sizeof(Vec));
                pixel->v->x = x;
                pixel->v->y = y;
                pixel->v->z = z;
                pixel->value = UINT8_MAX;
                pixel->dist = DBL_MAX;
                pixel->particle = getClosestParticle(pixel->v, cells, nCells);
            }
        }
    }

    return bitmap;
}

Bitmap* calcVoronoiDist(Bitmap* bitmap, Cell** cells, int nCells) {
    Pixel* pixel;
    double cDist;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);

                for (int i = 0; i < nCells; ++i) {
                    if (equalVecs(cells[i]->particle, pixel->particle)) {
                        for (int j = 0; j < cells[i]->faceCount; ++j) {
                            cDist = getDistPlane(pixel->v, &cells[i]->faces[j]);
                            if (cDist < pixel->dist)
                                pixel->dist = cDist;
                        }
                    }
                }
            }
        }
    }
    return bitmap;
}

void setValuesBitmap(Bitmap* bitmap, ImageParams* imageParams) {
    Pixel* pixel;
    double value;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                value = getPixelValue(pixel->dist, imageParams);
                pixel->value = to8Bit(value);
            }
        }
    }
}

Bitmap* setVoronoiValues(Cell** cells, int nCells, ImageParams* imageParams) {
    Bitmap* bitmap = initializeBitmap(cells, nCells, imageParams);
    printf("        Calculating distances\n");
    calcVoronoiDist(bitmap, cells, nCells);
    printf("        Setting bitmap values\n");
    setValuesBitmap(bitmap, imageParams);
    return bitmap;
}

void createVoronoiImage(Cell** cells, int nCells, ImageParams* imageParams, char* fname) {
    discretizeCells(cells, nCells, imageParams->imageSize);
    printf("    Calculating image values\n");
    Bitmap* bitmap = setVoronoiValues(cells, nCells, imageParams);
    printf("    Writing images\n");
    writeBitmap(bitmap, fname);
}
