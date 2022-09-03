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
                if (pixel->value >= value) ++(*count);
            }
        }
    }

    seeds = malloc(sizeof(Vec*) * *count);

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                if (pixel->value >= value) {
                    seeds[i] = pixel->v;
                    ++i;
                }
            }
        }
    }

    return seeds;
}

Vec** getSeedsWithBlockRad(Bitmap *bitmap, int value, int *count, double radius) {
    Vec** seeds = getSeeds(bitmap, value, count);

    for (int i = 0; i < *count; ++i) {
        for (int j = 0; j < *count; ++j) {
            if (getDist(seeds[i], seeds[j]) < radius) {
                if (j != *count - 1) seeds[j] = seeds[(*count) - 1];
                --(*count);
            }
        }
    }

    return seeds;
}

Vec** addBorderSeed(Vec** seeds, int *count) {
    ++(*count);
    seeds = realloc(seeds, sizeof(Vec*) * *count);
    Vec* vec = malloc(sizeof(Vec));
    vec->x = 0;
    vec->y = 0;
    vec->z = 0;
    seeds[(*count) - 1] = vec;
    return seeds;
}

Vec* getClosestParticle(Vec* v, Cell** cells, int nCells) {
    double sDist = DBL_MAX;
    double cDist;
    Vec* closest = NULL;
    for (int i = 0; i < nCells; ++i) {
        cDist = getDistSq(cells[i]->particle, v);
        if (cDist < sDist) {
            sDist = cDist;
            closest = cells[i]->particle;
        }
    }
    return closest;
}

double getLongestDist(Bitmap* bitmap) {
    double dist = 0;
    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                Pixel* pixel = getPixel(bitmap, x, y, z);
                if (pixel->dist < DBL_MAX && pixel->dist > dist) dist = pixel->dist;
            }
        }
    }
    return dist;
}

void setNoiseValuesBitmap(Bitmap* bitmap, ImageParams* imageParams) {
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

void setValuesBitmap(Bitmap* bitmap, ImageParams* imageParams) {
    Pixel* pixel;
    double value;
    int longest = (int)getLongestDist(bitmap) + 1;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                value = pixel->dist / longest * imageParams->distScalingFactor;
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

void createAgif(char *path) {
    char* command = malloc(sizeof(char) * (24 + 2 * strlen(path) + 1));
    sprintf(command, "convert %s/*.png %s/anim.gif\n", path, path);

    system(command);
    free(command);
}

void writeBitmap(Bitmap* bitmap, char* fname) {
    size_t len;
    char* activeName;

    for (int z = 0; z < bitmap->size; ++z) {
        len = snprintf(NULL, 0, "%simg_%03d.png", fname, z) + 1;
        activeName = malloc(len);
        snprintf(activeName, len, "%simg_%03d.png", fname, z);

        if (save_bitmap_slice_to_file(bitmap, z, activeName)) {
            fprintf(stderr, "Error writing file %d.\n", z);
        }
        free(activeName);
    }

    createAgif(fname);

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
                            if (cDist < pixel->dist) {
                                pixel->dist = cDist;
                                pixel->closestFace = getFaceCopy(&cells[i]->faces[j]);
                                pixel->cellId = cells[i]->id;
                            }
                        }
                    }
                }
            }
        }
    }
    return bitmap;
}

Bitmap* removeDistMissingFaces(Bitmap* bitmap, Cell** cells, int nCells) {
    Pixel* pixel;
    int found;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                found = 0;

                for (int i = 0; i < nCells; ++i) {
                    for (int j = 0; j < cells[i]->faceCount; ++j) {
                        if (equalFaces(pixel->closestFace, &cells[i]->faces[j])) {
                            found = 1;
                            break;
                        }
                    }
                    if (found) break;
                }

                if (!found) {
                    pixel->dist = DBL_MAX;
                    pixel->closestFace = NULL;
                }
            }
        }
    }
    return bitmap;
}

Bitmap* calcMissingDist(Bitmap* bitmap, Cell** cells, int nCells) {
    Pixel* pixel;
    int i;
    int found;

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);

                if (pixel->dist == DBL_MAX) {
                    found = 0;

                    for (i = 0; i < nCells; ++i) {
                        for (int j = 0; j < cells[i]->nIds; ++j) {
                            if (cells[i]->ids[j] == pixel->cellId) {
                                found = 1;
                                break;
                            }
                        }
                        if (found) break;
                    }
                    pixel->dist = getDistCell(pixel->v, cells[i]);
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

void printBitmapPixelDist(Bitmap* bitmap) {
    int *counts = malloc(sizeof(int) * 256);

    for (int i = 0; i < 256; ++i) {
        counts[i] = 0;
    }

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                Pixel *pixel = getPixel(bitmap, x, y, z);
                ++counts[pixel->value];
            }
        }
    }

    printf("Bitmap pixel value counts: \n");
    for (int i = 0; i < 256; ++i) {
        printf("    %d :  %d\n", i, counts[i]);
    }
}


/* *********************************************************************************************************************
 * *********************************************************************************************************************
 * ***************************************************   MEASURES   ****************************************************
 * *********************************************************************************************************************
 * *********************************************************************************************************************
 */

void setParticleIds(Bitmap* bitmap, int n) {
    Vec** mapping = malloc(sizeof(Vec*) * n);

    for (int i = 0; i < n; ++i) {
        mapping[i] = NULL;
    }

    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                Pixel* pixel = getPixel(bitmap, x, y, z);

                for (int i = 0; i < n; ++i) {
                    if (mapping[i] == NULL) {
                        mapping[i] = pixel->particle;
                        pixel->particleId = i;
                        break;
                    }
                    if (mapping[i] == pixel->particle) {
                        pixel->particleId = i;
                        break;
                    }
                }
            }
        }
    }
    free(mapping);
}

double getRandsIndex(Bitmap* orig, Bitmap* srg, int nParticles, int nSeeds) {
    Pixel *pSRG, *pOrig;

    size_t n = orig->size * orig->size * orig->size;
    n = n * (n - 1) / 2;

    size_t *nMap = malloc(sizeof(size_t) * nParticles * nSeeds);
    memset(nMap, 0, sizeof(size_t) * nParticles * nSeeds);

    size_t *truthSum = malloc(sizeof(size_t) * nParticles);
    memset(truthSum, 0, sizeof(size_t) * nParticles);

    size_t *predSum = malloc(sizeof(size_t) * nSeeds);
    memset(predSum, 0, sizeof(size_t) * nSeeds);

    size_t falseJoins = 0, falseCuts = 0, trueJoins = 0, trueCuts = 0;

    for (int z = 0; z < orig->size; ++z) {
        for (int y = 0; y < orig->size; ++y) {
            for (int x = 0; x < orig->size; ++x) {

                pOrig = getPixel(orig, x, y, z);
                pSRG = getPixel(srg, x, y, z);

                ++nMap[pOrig->particleId + nParticles * pSRG->particleId];
                ++truthSum[pOrig->particleId];
                ++predSum[pSRG->particleId];
            }
        }
    }

    for (int i = 0; i < nSeeds; ++i) {
        falseJoins += predSum[i] * predSum[i];
    }

    for (int i = 0; i < nParticles; ++i) {
        falseCuts += truthSum[i] * truthSum[i];
    }

    for (int i = 0; i < nParticles * nSeeds; ++i) {
        size_t n_ij = nMap[i];
        if (n_ij == 0) continue;

        trueJoins += n_ij * (n_ij - 1) / 2;
        falseCuts -= n_ij * n_ij;
        falseJoins -= n_ij * n_ij;
    }

    falseJoins /= 2;
    falseCuts /= 2;

    trueCuts = n - (trueJoins + falseJoins) - falseCuts;

    free(nMap);
    free(truthSum);
    free(predSum);

    return (double)(trueJoins + trueCuts) / (double)n;
}

double getVariationOfInformation(Bitmap* orig, Bitmap* srg, int nParticles, int nSeeds) {

    size_t n = orig->size * orig->size * orig->size;

    double *nMap = malloc(sizeof(double) * nParticles * nSeeds);
    memset(nMap, 0, sizeof(double) * nParticles * nSeeds);

    double *truthSum = malloc(sizeof(double) * nParticles);
    memset(truthSum, 0, sizeof(double) * nParticles);

    double *predSum = malloc(sizeof(double) * nSeeds);
    memset(predSum, 0, sizeof(double) * nSeeds);

    double H0 = 0, H1 = 0, I = 0;
    int l1, l2;


    for (int z = 0; z < orig->size; ++z) {
        for (int y = 0; y < orig->size; ++y) {
            for (int x = 0; x < orig->size; ++x) {
                Pixel *pSRG = getPixel(srg, x, y, z);
                Pixel *pOrig = getPixel(orig, x, y, z);

                ++nMap[pOrig->particleId + nParticles * pSRG->particleId];
                ++truthSum[pOrig->particleId];
                ++predSum[pSRG->particleId];
            }
        }
    }

    for (int i = 0; i < nParticles; ++i) {
        truthSum[i] /= n;
    }

    for (int i = 0; i < nSeeds; ++i) {
        predSum[i] /= n;
    }

    for (int i = 0; i < nParticles * nSeeds; ++i) {
        nMap[i] /= n;
    }


    for (int i = 0; i < nParticles; ++i) {
        if (isZero(truthSum[i])) continue;
        H0 -= truthSum[i] * log2(truthSum[i]);
    }

    for (int i = 0; i < nSeeds; ++i) {
        if (isZero(predSum[i])) continue;
        H1 -= predSum[i] * log2(predSum[i]);
    }

    double div;

    for (int i = 0; i < nParticles * nSeeds; ++i) {
        if (isZero(nMap[i])) continue;

        l1 = i % nParticles;
        l2 = i / nParticles;

        div = nMap[i]  / (truthSum[l1] * predSum[l2]);

        if (isZero(div)) continue;

        I += nMap[i] * log2(div);
    }

    free(nMap);
    free(truthSum);
    free(predSum);

    return H0 + H1 - 2.0 * I;
}
