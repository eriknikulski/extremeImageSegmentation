//
// Created by Erik Nikulski on 07.07.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_PARAMS_H
#define EXTREMEIMAGESEGMENTATION_PARAMS_H

typedef struct ImageParams {
    int imageSize;
    double theta_0;
    double theta_1;
    double sigma_b;
    double mu_b;
    double sigma_c;
    double mu_c;
} ImageParams;

typedef struct VoronoiParams {
    int nInitialCells;
    int nCells;
    int srgPrecision;
    char* imagePath;
    char* srgImagePath;
} VoronoiParams;

typedef struct SplineParams {
    double alpha;
    double minDist;
    int nPoints;
    int nSplines;
    int srgPrecision;
    char* imagePath;
    char* srgImagePath;
} SplineParams;

#endif //EXTREMEIMAGESEGMENTATION_PARAMS_H