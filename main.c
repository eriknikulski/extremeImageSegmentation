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
//    Vec* seeds = getSeeds(bitmapOrig, splineParams->seedThreshold, &nSeeds);
//    printf("nSeeds: %d\n", nSeeds);

//    nSeeds = splineParams->nSplines;
//    Vec* seeds = malloc(sizeof(Vec) * nSeeds);
//    for (int i = 0; i < splineParams->nSplines; ++i) {
//        seeds[i] = s[i][splineParams->nPoints / 2];
//    }

    nSeeds = splineParams->nSplines + 1;
    Vec* seeds = malloc(sizeof(Vec) * nSeeds);
    for (int i = 0; i < splineParams->nSplines; ++i) {
        seeds[i] = s[i][splineParams->nPoints / 2];
    }
    for (int i = 0; i < imageParams->imageSize; ++i) {
        Pixel* p = getPixel(bitmapOrig, i, 0, 0);
        if (p->particle == NULL) {
            seeds[nSeeds - 1] = *p->v;
            break;
        }
    }

    printf("Applying srg\n");
    Bitmap* bitmapSRG = initializeBitmap(imageParams);
    bitmapSRG->pixels = memcpy(bitmapSRG->pixels, bitmapOrig->pixels, sizeof(Pixel) * imageParams->imageSize * imageParams->imageSize * imageParams->imageSize);
    bitmapSRG = srg(bitmapSRG, seeds, nSeeds, splineParams->srgPrecision,
                    splineParams->imagePath, imageParams);
    printf("Writing bitmap\n");
    writeBitmap(bitmapSRG, splineParams->srgImagePath);

    printf("Calculating Metrics: \n");
    setParticleIds(bitmapOrig, splineParams->nSplines);
    setParticleIds(bitmapSRG, nSeeds);

    double randsIndex = getRandsIndex(bitmapOrig, bitmapSRG, splineParams->nSplines, nSeeds);
    printf("Rands Index: %lf\n", randsIndex);
    double variationOfInformation = getVariationOfInformation(bitmapOrig, bitmapSRG, splineParams->nSplines, nSeeds);
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

//    setNoiseValuesBitmap(bitmapOrig, imageParams);
//    writeBitmap(bitmapOrig, voronoiParams->imagePath);

//    printf("Creating voronoi images\n");
//    Bitmap* bitmapOrig = createVoronoiImage(cells, voronoiParams, imageParams);

//    int nSeeds;
//    Vec** seeds = getSeeds(bitmapOrig, voronoiParams->seedThreshold, &nSeeds);
//    printf("nSeeds: %d\n", nSeeds);

//    int nSeeds = voronoiParams->nInitialCells;

    int nSeeds = voronoiParams->nInitialCells + 1;
    particles = realloc(particles, sizeof(Vec) * (voronoiParams->nInitialCells + 1));
    particles[voronoiParams->nInitialCells].x = 0.0;
    particles[voronoiParams->nInitialCells].y = 0.0;
    particles[voronoiParams->nInitialCells].z = 0.0;

    printf("Applying srg\n");
    Bitmap* bitmapSRG = initializeBitmap(imageParams);
    bitmapSRG->pixels = memcpy(bitmapSRG->pixels, bitmapOrig->pixels, sizeof(Pixel) * imageParams->imageSize * imageParams->imageSize * imageParams->imageSize);
    bitmapSRG = srg(bitmapSRG, particles, nSeeds, voronoiParams->srgPrecision,
                            voronoiParams->imagePath, imageParams);
    printf("Writing bitmap\n");
    writeBitmap(bitmapSRG, voronoiParams->srgImagePath);

    printf("Calculating Metrics: \n");
    setParticleIds(bitmapOrig, voronoiParams->nCells);
    setParticleIds(bitmapSRG, nSeeds);

    double randsIndex = getRandsIndex(bitmapOrig, bitmapSRG, voronoiParams->nCells, nSeeds);
    printf("Rands Index: %lf\n", randsIndex);

    double variationOfInformation = getVariationOfInformation(bitmapOrig, bitmapSRG, voronoiParams->nCells, nSeeds);
    printf("Variation Of Information: %lf\n", variationOfInformation);
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
    };

    VoronoiParams voronoiParams = {
            .nInitialCells = 4,
            .nCells = 4,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/fin/"),
            .imagePathPre = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/pre/"),
            .imagePathRem = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/rem/"),
            .imagePathMerged = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/merged/"),
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
