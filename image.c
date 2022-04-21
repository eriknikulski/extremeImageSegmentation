//
// Created by Erik Nikulski on 11.04.22.
//

#include "image.h"

#include "utility.h"
#include "spline.h"
#include "voronoi.h"

#include "math.h"
#include "stdio.h"

static double getPixelValue(double d, double theta_0, double theta_1, double sigma_b, double mu_b, double sigma_c, double mu_c) {
    double p = (double)1.0 / (double)((double)1 + exp(-theta_0 - theta_1 * d));
    if (getRand() < p) {
        // select from N2
        return normDist(sigma_b, mu_b);
    } else {
        // select from N1
        return normDist(sigma_c, mu_c);
    }
}

static unsigned char to8Bit(double value) {
    if (value < 0) value = 0;
    if (value > 1) value = 1;

    return (unsigned char) (value * 255);
}

static void writeImage(int size, void* image, char* fname) {
    FILE* fp;
    fp = fopen(fname, "wb");
    fwrite(&image, size * size * size, sizeof(unsigned char), fp);
    fclose(fp);
}

void createSplineImage(Vec** splines, int nSplines, int dim, int imageSize,
                 double theta_0, double theta_1, double sigma_b, double mu_b, double sigma_c, double mu_c,
                 char* fname) {
    discretizeSplines(splines, nSplines, dim, imageSize);
    Vec p;
    double dist;
    double value;

    unsigned char image[imageSize][imageSize][imageSize];

    for (int i = 0; i < imageSize; ++i) {
        for (int j = 0; j < imageSize; ++j) {
            for (int k = 0; k < imageSize; ++k) {
                p.x = (double) i;
                p.y = (double) j;
                p.z = (double) k;

                dist = splinesDist(&p, splines, nSplines, dim);
                value = getPixelValue(dist, theta_0, theta_1, sigma_b, mu_b, sigma_c, mu_c);
                image[i][j][k] = to8Bit(value);
            }
        }
    }

    writeImage(imageSize, &image, fname);
    printf("Created spline image!\n");
}

void createVoronoiImage(Cell** cells, int nCells, int imageSize,
                        double theta_0, double theta_1, double sigma_b, double mu_b, double sigma_c, double mu_c,
                        char* fname) {
    Vec p;
    double dist;
    double value;

    unsigned char image[imageSize][imageSize][imageSize];

    discretizeCells(cells, nCells, imageSize);

    for (int i = 0; i < imageSize; ++i) {
        for (int j = 0; j < imageSize; ++j) {
            for (int k = 0; k < imageSize; ++k) {
                p.x = (double) i;
                p.y = (double) j;
                p.z = (double) k;

                dist = voronoiDist(&p, cells, nCells);
                value = getPixelValue(dist, theta_0, theta_1, sigma_b, mu_b, sigma_c, mu_c);
                image[i][j][k] = to8Bit(value);
            }
        }
    }

    writeImage(imageSize, &image, fname);
    printf("Created voronoi image!\n");
}
