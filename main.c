#include "image.h"
#include "params.h"
#include "spline.h"
#include "voronoi.h"
#include "srg.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "time.h"


void spline(SplineParams* splineParams, ImageParams* imageParams) {
    Vec** s = getNSplines(splineParams);
    printf("Created splines\n");
    Bitmap* bitmapOrig = createSplineImage(s, splineParams, imageParams);
    printf("Created spline images!\n");

    int nSeeds;
    Vec** seeds = getSeeds(bitmapOrig, splineParams->seedThreshold, &nSeeds);
    printf("nSeeds: %d\n", nSeeds);

    printf("Applying srg\n");
    Bitmap* bitmapSRG = srg(bitmapOrig, seeds, nSeeds, splineParams->srgPrecision,
                            splineParams->imagePath, imageParams);
    printf("Writing bitmap\n");
    writeBitmap(bitmapSRG, splineParams->srgImagePath);

    double randsIndex = getRandsIndex(bitmapOrig, bitmapSRG);
    printf("Rands Index: %lf\n", randsIndex);
    double variationOfInformation = getVariationOfInformation(bitmapOrig, bitmapSRG, seeds, nSeeds);
    printf("Variation Of Information: %lf\n", variationOfInformation);
}

void voronoi(VoronoiParams* voronoiParams, ImageParams* imageParams) {
    Vec* particles = malloc(sizeof(Vec) * voronoiParams->nInitialCells);
    printf("Creating cells\n");
    Cell** cells = getCells(voronoiParams->nInitialCells, "test", &particles);

//    Vec* particles = malloc(sizeof(Vec) * voronoiParams->nCells);
//    printf("Creating cells\n");
//    Cell** cells = getCells(voronoiParams->nCells, "test", &particles);

//    printf("Particles:\n");
//    printVecs(particles, voronoiParams->nInitialCells);
//    printVecs(particles, voronoiParams->nCells);

    discretizeCells(cells, voronoiParams->nInitialCells, imageParams->imageSize);
    Bitmap* bitmapOrig = initializeBitmap(imageParams);
    calcVoronoiDist(bitmapOrig, cells, voronoiParams->nInitialCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePathPre);

    cells = mergeCells(cells, voronoiParams);
    printf("Merged cells\n\n");
    printCells(cells, voronoiParams->nCells);

    removeDistMissingFaces(bitmapOrig, cells, voronoiParams->nCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePathRem);

    calcMissingDist(bitmapOrig, cells, voronoiParams->nCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePath);

//    printf("Creating voronoi images\n");
//    Bitmap* bitmapOrig = createVoronoiImage(cells, voronoiParams, imageParams);

//    int nSeeds;
//    Vec** seeds = getSeeds(bitmapOrig, voronoiParams->seedThreshold, &nSeeds);
//    printf("nSeeds: %d\n", nSeeds);
//
//    printf("Applying srg\n");
//    Bitmap* bitmapSRG = srg(bitmapOrig, seeds, nSeeds, voronoiParams->srgPrecision,
//                            voronoiParams->imagePath, imageParams);
//    printf("Writing bitmap\n");
//    writeBitmap(bitmapSRG, voronoiParams->srgImagePath);
//
//    double randsIndex = getRandsIndex(bitmapOrig, bitmapSRG);
//    printf("Rands Index: %lf\n", randsIndex);
//    double variationOfInformation = getVariationOfInformation(bitmapOrig, bitmapSRG, seeds, nSeeds);
//    printf("Variation Of Information: %lf\n", variationOfInformation);
}

int main() {
    srand(time(0));

    // distinct
    ImageParams imageParams = {
            .imageSize = 128,
            .theta_0 = -3.0,
            .theta_1 = 1,
            .sigma_b = 5,
            .mu_b = 3,
            .sigma_c = 5,
            .mu_c = -3,
    };

    VoronoiParams voronoiParams = {
            .nInitialCells = 4,
            .nCells = 2,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/"),
            .imagePathPre = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/pre/"),
            .imagePathRem = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/rem/"),
            .srgImagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/srg/"),
            .srgPrecision = 100,
            .seedThreshold = 54,
    };

    SplineParams splineParams = {
            .alpha = 0.5,
            .minDist = 0.1,
            .nPoints = 128,
            .nSplines = 5,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/"),
            .srgImagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/srg/"),
            .srgPrecision = 100,
            .seedThreshold = 80,
    };

    clock_t begin;
    clock_t end;
    double time_spent;


    begin = clock();

    spline(&splineParams, &imageParams);
    voronoi(&voronoiParams, &imageParams);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME: %lf\n", time_spent);

    // TODO: cgal

    return 0;
}
