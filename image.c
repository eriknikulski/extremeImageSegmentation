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

Vec** getSeeds(Bitmap* bitmap, uint8_t value, int* count) {
    Pixel* pixel;
    Vec** seeds;
    *count = 0;
    int i = 0;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                if (pixel->value <= value) ++(*count);
            }
        }
    }

    seeds = malloc(sizeof(Vec*) * *count);

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                if (pixel->value <= value) {
                    seeds[i] = pixel->v;
                    ++i;
                }
            }
        }
    }

    return seeds;
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

void setValuesBitmap(Bitmap* bitmap, ImageParams* imageParams) {
    Pixel* pixel;
    double value;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                value = getPixelValue(pixel->dist, imageParams);
                if (bitmap->reverse) value = 1 - value;
                pixel->value = to8Bit(value);
            }
        }
    }
}

Bitmap* initializeBitmap(ImageParams* imageParams) {
    Pixel* pixel;
    Bitmap* bitmap = malloc(sizeof(Bitmap));
    bitmap->reverse = 0;
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
                pixel->particle = NULL;
            }
        }
    }

    return bitmap;
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

Bitmap* calcSplineDistance(Bitmap* bitmap, Vec** splines, SplineParams* splineParams, ImageParams* imageParams) {
    Vec p;
    double dist;
    int spline;

    for (int z = 0; z < imageParams->imageSize; ++z) {
        for (int y = 0; y < imageParams->imageSize; ++y) {
            for (int x = 0; x < imageParams->imageSize; ++x) {
                p.x = (double) x;
                p.y = (double) y;
                p.z = (double) z;

                Pixel* pixel = getPixel(bitmap, x, y, z);
                dist = splinesDist(&p, splines, splineParams, &spline);
                pixel->dist = dist;
                if (dist <= 2)
                    pixel->particle = &splines[spline][splineParams->nPoints / 2];
            }
        }
    }
    return bitmap;
}

Bitmap* setSplineValues(Vec** splines, SplineParams* splineParams, ImageParams* imageParams) {
    Bitmap* bitmap = initializeBitmap(imageParams);
    bitmap->reverse = 1;
    printf("        Calculating distances\n");
    bitmap = calcSplineDistance(bitmap, splines, splineParams, imageParams);
    printf("        Setting bitmap values\n");
    setValuesBitmap(bitmap, imageParams);
    return bitmap;
}

Bitmap* createSplineImage(Vec** splines, SplineParams* splineParams, ImageParams* imageParams) {
    discretizeSplines(splines, splineParams->nSplines, splineParams->nPoints, imageParams->imageSize);
    printf("    Calculating image values\n");
    Bitmap* bitmap = setSplineValues(splines, splineParams, imageParams);
    printf("    Writing images\n");
    writeBitmap(bitmap, splineParams->imagePath);
    return bitmap;
}

Bitmap* calcVoronoiDist(Bitmap* bitmap, Cell** cells, int nCells) {
    Pixel* pixel;
    double cDist;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                pixel->particle = getClosestParticle(pixel->v, cells, nCells);

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

Bitmap* setVoronoiValues(Cell** cells, int nCells, ImageParams* imageParams) {
    Bitmap* bitmap = initializeBitmap(imageParams);
    printf("        Calculating distances\n");
    calcVoronoiDist(bitmap, cells, nCells);
    printf("        Setting bitmap values\n");
    setValuesBitmap(bitmap, imageParams);
    return bitmap;
}

Bitmap* createVoronoiImage(Cell** cells, VoronoiParams* voronoiParams, ImageParams* imageParams) {
    discretizeCells(cells, voronoiParams->nCells, imageParams->imageSize);
    printf("    Calculating image values\n");
    Bitmap* bitmap = setVoronoiValues(cells, voronoiParams->nCells, imageParams);
    printf("    Writing images\n");
    writeBitmap(bitmap, voronoiParams->imagePath);
    return bitmap;
}


/* *********************************************************************************************************************
 * *********************************************************************************************************************
 * ***************************************************   MEASURES   ****************************************************
 * *********************************************************************************************************************
 * *********************************************************************************************************************
 */

double getRandsIndex(Bitmap* orig, Bitmap* srg) {
    Pixel* pSRG;
    Pixel* pOrig;
    int counter = 0;
    int n = orig->size * orig->size * orig->size;

    for (int z = 0; z < orig->size; ++z) {
        for (int y = 0; y < orig->size; ++y) {
            for (int x = 0; x < orig->size; ++x) {
                pSRG = getPixel(srg, x, y, z);
                pOrig = getPixel(orig, x, y, z);
                if (pSRG->particle == pOrig->particle) ++counter;
            }
        }
    }
    return (double)counter / (double)n;
}

double getVariationOfInformation(Bitmap* orig, Bitmap* srg, Vec** particles, int n) {
    int *p = malloc(sizeof(int) * (n + 1));
    int *q = malloc(sizeof(int) * (n + 1));
    int *r = malloc(sizeof(int) * (n + 1));
    memset(p, 0, sizeof(int) * (n + 1));
    memset(q, 0, sizeof(int) * (n + 1));
    memset(r, 0, sizeof(int) * (n + 1));
    int size = orig->size * orig->size * orig->size;
    double res = 0;
    Pixel* pSRG;
    Pixel* pOrig;

    for (int z = 0; z < orig->size; ++z) {
        for (int y = 0; y < orig->size; ++y) {
            for (int x = 0; x < orig->size; ++x) {
                pSRG = getPixel(srg, x, y, z);
                pOrig = getPixel(orig, x, y, z);

                if (!pOrig->particle) ++p[n];
                if (!pOrig->particle && !pSRG->particle) ++r[n];
                if (!pSRG->particle) ++q[n];

                for (int i = 0; i < n; ++i) {
                    if (equalVecs(pSRG->particle, particles[i])) ++q[i];
                    if (equalVecs(pOrig->particle, particles[i])) ++p[i];
                    if (equalVecs(pOrig->particle, particles[i]) && equalVecs(pSRG->particle, particles[i]))
                        ++r[i];
                }
            }
        }
    }

    for (int i = 0; i < n + 1; ++i) {
        if (r[i] == 0 || p[i] == 0 || q[i] == 0) continue;
        res -= (double)r[i] / size * (log((double)r[i] / p[i]) + log((double)r[i] / q[i]));
    }
    
    return res;
}
