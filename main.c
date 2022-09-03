#include "image.h"
#include "params.h"
#include "spline.h"
#include "voronoi.h"
#include "srg.h"

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "time.h"


void spline(SplineParams* splineParams, ImageParams* imageParams) {
    imageParams->distScalingFactor = splineParams->imageDistScalingFactor;
    Vec** s = getNSplines(splineParams);
    printf("Created splines\n");
    Bitmap* bitmapOrig = createSplineImage(s, splineParams, imageParams);
    printf("Created spline images!\n");

    setNoiseValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, splineParams->imagePath);

    printBitmapPixelDist(bitmapOrig);

    Vec *particles = getParticles(s, bitmapOrig, splineParams, imageParams);
    calcSeedValueMetrics(bitmapOrig, imageParams, splineParams->statsPath,
                         splineParams->srgImagePath,splineParams->imagePath,
                         splineParams->nSplines, splineParams->blockingRadius,
                         splineParams->srgPrecision, particles);
}

void voronoi(VoronoiParams* voronoiParams, ImageParams* imageParams) {
    imageParams->distScalingFactor = voronoiParams->imageDistScalingFactor;

    Vec* particles = malloc(sizeof(Vec) * voronoiParams->nInitialCells);
    printf("Creating cells\n");
    Cell** cells = getCells(voronoiParams->nInitialCells, "test", &particles);

    discretizeCells(cells, voronoiParams->nInitialCells, imageParams->imageSize);
    Bitmap* bitmapOrig = initializeBitmap(imageParams);
    calcVoronoiDist(bitmapOrig, cells, voronoiParams->nInitialCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePathPre);

    printf("Cells\n\n");
    printCells(cells, voronoiParams->nInitialCells);
    printf("\n\n\n");

    cells = mergeCells(cells, voronoiParams);
    printf("Merged cells\n\n");
    printCells(cells, voronoiParams->nCells);

    removeDistMissingFaces(bitmapOrig, cells, voronoiParams->nCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePathRem);

    calcMissingDist(bitmapOrig, cells, voronoiParams->nCells);
    setValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePathMerged);

    setNoiseValuesBitmap(bitmapOrig, imageParams);
    writeBitmap(bitmapOrig, voronoiParams->imagePath);

    printBitmapPixelDist(bitmapOrig);

    calcSeedValueMetrics(bitmapOrig, imageParams, voronoiParams->statsPath,
                         voronoiParams->srgImagePath,voronoiParams->imagePath,
                         voronoiParams->nCells, voronoiParams->blockingRadius,
                         voronoiParams->srgPrecision, particles);
}

int main() {
    srand(time(0));

    // distinct
    ImageParams imageParams = {
            .imageSize = 128,
            .theta_0 = -2.0,
            .theta_1 = 0.5,
            .sigma_b = (double) 255 / 255,
            .mu_b = (double) 434 / 255,
            .sigma_c = sqrt((double) 70 / 255),
            .mu_c = (double) 25.5 / 255,
            .distScalingFactor = 1.0,
    };

    VoronoiParams voronoiParams = {
            .nInitialCells = 16,
            .nCells = 4,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/fin/"),
            .imagePathPre = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/pre/"),
            .imagePathRem = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/rem/"),
            .imagePathMerged = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/merged/"),
            .srgImagePathBase = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/srg/"),
            .srgImagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/srg/"),
            .statsPath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/voronoi_stats.csv"),
            .srgPrecision = 100,
            .seedThreshold = 54,
            .blockingRadius = 40,
            .imageDistScalingFactor = 1.0,
    };

    SplineParams splineParams = {
            .alpha = 0.5,
            .minDist = 0.1,
            .nPoints = 128,
            .nSplines = 5,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/"),
            .srgImagePathBase = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/srg/"),
            .srgImagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/srg/"),
            .statsPath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/spline_stats.csv"),
            .srgPrecision = 100,
            .seedThreshold = 80,
            .blockingRadius = 40,
            .imageDistScalingFactor = 5.0,
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

    return 0;
}
