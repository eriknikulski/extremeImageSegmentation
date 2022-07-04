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


/* Given "bitmap", this returns the pixel of bitmap at the point
   ("x", "y", "z"). */

Pixel* getPixel(Bitmap* bitmap, int x, int y, int z) {
    return bitmap->pixels + bitmap->size2 * z + bitmap->size * y + x;
}

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

Bitmap* calcSplineDistance(Vec** splines, int nSplines, int dim, ImageParams imageParams) {
    Vec p;
    double dist;
    double value;

    Bitmap* bitmap = malloc(sizeof(Bitmap));
    bitmap->size = imageParams.imageSize;
    bitmap->size2 = imageParams.imageSize * imageParams.imageSize;
    bitmap->pixels = malloc(sizeof(Pixel) * imageParams.imageSize * imageParams.imageSize * imageParams.imageSize);

    for (int z = 0; z < imageParams.imageSize; ++z) {
        for (int y = 0; y < imageParams.imageSize; ++y) {
            for (int x = 0; x < imageParams.imageSize; ++x) {
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

void createSplineImage(Vec** splines, int nSplines, int dim, ImageParams imageParams, char* fname) {
    discretizeSplines(splines, nSplines, dim, imageParams.imageSize);
    Bitmap* bitmap = calcSplineDistance(splines, nSplines, dim, imageParams);
    writeBitmap(bitmap, fname);
}

Bitmap* calcVoronoiDistances(Cell** cells, int nCells, ImageParams imageParams) {
    Vec p;
    double dist;
    double value;

    Bitmap* bitmap = malloc(sizeof(Bitmap));
    bitmap->size = imageParams.imageSize;
    bitmap->size2 = imageParams.imageSize * imageParams.imageSize;
    bitmap->pixels = malloc(sizeof(Pixel) * imageParams.imageSize * imageParams.imageSize * imageParams.imageSize);

    for (int z = 0; z < imageParams.imageSize; ++z) {
        for (int y = 0; y < imageParams.imageSize; ++y) {
            for (int x = 0; x < imageParams.imageSize; ++x) {
                p.x = (double) x;
                p.y = (double) y;
                p.z = (double) z;

                dist = voronoiDist(&p, cells, nCells);
                value = getPixelValue(dist, imageParams);
                Pixel* pixel = getPixel(bitmap, x, y, z);
                pixel->value = to8Bit(value);
            }
        }
        printf("%d\n", z);
    }
    return bitmap;
}

void createVoronoiImage(Cell** cells, int nCells, ImageParams imageParams, char* fname) {
    discretizeCells(cells, nCells, imageParams.imageSize);
    Bitmap* bitmap = calcVoronoiDistances(cells, nCells, imageParams);
    writeBitmap(bitmap, fname);
}
