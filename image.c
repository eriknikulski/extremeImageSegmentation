//
// Created by Erik Nikulski on 11.04.22.
//

#include "image.h"

#include "utility.h"
#include "spline.h"
#include "voronoi.h"
#include "png.h"


#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

static double getPixelValue(double d, ImageParams imageParams) {
    double p = (double)1.0 / (double)((double)1 + exp(-imageParams.theta_0 - imageParams.theta_1 * d));

    // TODO: N1 and N2 in task are diff for splines and voronoi. Does it matter??
    if (getRandD() > p) {
        // select from N2
        return normalCDF(getRandD(), imageParams.sigma_b, imageParams.mu_b);
    } else {
        // select from N1
        return normalCDF(getRandD(), imageParams.sigma_c, imageParams.mu_c);
    }
}

static uint8_t to8Bit(double value) {
    if (value < 0) value = 0;
    if (value > 1) value = 1;

    return (uint8_t) (value * 255);
}

static int writeImage(bitmap_t* img, char* fname, int x) {
    int status = 0;
    size_t len = snprintf(NULL, 0, "%simg_%d.png", fname, x) + 1;
    char* activeName = malloc(len);
    snprintf(activeName, len, "%simg_%d.png", fname, x);

    if (save_png_to_file(img, activeName)) {
        fprintf(stderr, "Error writing file %d.\n", x);
        status = -1;
    }

    free(activeName);
    return status;
}

void createSplineImage(Vec** splines, int nSplines, int dim, int imageSize, ImageParams imageParams, char* fname) {
    discretizeSplines(splines, nSplines, dim, imageSize);
    Vec p;
    double dist;
    double value;

    bitmap_t img;

    for (int x = 0; x < imageSize; ++x) {
        img.width = imageSize;
        img.height = imageSize;

        img.pixels = calloc(img.width * img.height, sizeof(pixel_t));
        for (int y = 0; y < imageSize; ++y) {
            for (int z = 0; z < imageSize; ++z) {
                p.x = (double) x;
                p.y = (double) y;
                p.z = (double) z;

                dist = splinesDist(&p, splines, nSplines, dim);
                value = getPixelValue(dist, imageParams);
                pixel_t* pixel = pixel_at(&img, y, z);
                pixel->value = to8Bit(value);
            }
        }
        writeImage(&img, fname, x);
        free(img.pixels);
        printf("%d\n", x);
    }
}

void createVoronoiImage(Cell** cells, int nCells, int imageSize, ImageParams imageParams, char* fname) {
    Vec p;
    Vec lastP;
    double dist;
    double lastDist;
    double secondLastDist = -1;
    double value;

    bitmap_t img;

    discretizeCells(cells, nCells, imageSize - 1);

    printCells(cells, nCells);

    for (int x = 0; x < imageSize; ++x) {
        img.width = imageSize;
        img.height = imageSize;

        img.pixels = calloc(img.width * img.height, sizeof(pixel_t));
        for (int y = 0; y < imageSize; ++y) {
            for (int z = 0; z < imageSize; ++z) {
                p.x = (double) x;
                p.y = (double) y;
                p.z = (double) z;

                // TODO: something is still wrong here
                if (secondLastDist >= 0 &&
                    lastDist + getDist(&p, &lastP) <= secondLastDist - getDist(&p, &lastP)) {
                    dist = lastDist + getDist(&p, &lastP);
                    secondLastDist = secondLastDist - getDist(&p, &lastP);
                } else {
                    dist = voronoiDist(&p, cells, nCells, &secondLastDist);
                }
                lastDist = dist;

                value = getPixelValue(dist, imageParams);
                pixel_t* pixel = pixel_at(&img, y, z);
                pixel->value = to8Bit(value);
                lastP = p;
            }
        }
        writeImage(&img, fname, x);
        free(img.pixels);
        printf("%d\n", x);
    }
}
