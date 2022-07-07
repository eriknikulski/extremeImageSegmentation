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
    createSplineImage(s, splineParams, imageParams);
    printf("Created spline images!\n");
}

void voronoi(VoronoiParams* voronoiParams, ImageParams* imageParams) {
//    Vec* particles = malloc(sizeof(Vec) * voronoiParams->nInitialCells);
//    Cell** cells = getCells(voronoiParams->nInitialCells, "test", &particles);

    Vec* particles = malloc(sizeof(Vec) * voronoiParams->nCells);
    printf("Creating cells\n");
    Cell** cells = getCells(voronoiParams->nCells, "test", &particles);

//    printf("Particles:\n");
//    printVecs(particles, voronoiParams->nInitialCells);
//    printVecs(particles, voronoiParams->nCells);

    // TODO: when cells are merged, seeds need to be merged (choose one)
    // TODO: when merging: some distances need to be recalculated (careful)
//    cells = mergeCells(cells, voronoiParams);
//    printf("Merged cells\n\n");
//    printCells(cells, voronoiParams->nCells);

    printf("Creating voronoi images\n");
    createVoronoiImage(cells, voronoiParams, imageParams);

    printf("Particles:\n");
    printVecs(particles, voronoiParams->nCells);
//
    printf("Applying srg\n");
    Bitmap* bitmap = srg(particles, voronoiParams, imageParams);
    printf("Writing bitmap\n");
    writeBitmap(bitmap, voronoiParams->srgImagePath);
}

int main() {
    srand(time(0));

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
            .nInitialCells = 16,
            .nCells = 4,
            .srgPrecision = 100,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/"),
            .srgImagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/srg/"),
    };

    SplineParams splineParams = {
            .alpha = 0.5,
            .minDist = 0.1,
            .nPoints = 100,
            .nSplines = 5,
            .imagePath = strdup("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/"),
    };

    clock_t begin;
    clock_t end;
    double time_spent;


    begin = clock();

//    spline(&splineParams, &imageParams);

    voronoi(&voronoiParams, &imageParams);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME: %lf\n", time_spent);

    // TODO: use initial points for cell creation as  seeds; and then one for cell walls
    // TODO: for image creation: dont calc dist from pixel to cells; calc cell
    // TODO: cgal

    return 0;
}
