//
// Created by Erik Nikulski on 11.04.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_IMAGE_H
#define EXTREMEIMAGESEGMENTATION_IMAGE_H

#include "vec.h"
#include "voronoi.h"

void createSplineImage(Vec** splines, int nSplines, int dim, int imageSize,
                 double theta_0, double theta_1, double sigma_b, double mu_b, double sigma_c, double mu_c,
                 char* fname);

void createVoronoiImage(Cell** cells, int nCells, int imageSize,
                        double theta_0, double theta_1, double sigma_b, double mu_b, double sigma_c, double mu_c,
                        char* fname);

#endif //EXTREMEIMAGESEGMENTATION_IMAGE_H
