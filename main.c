#include "image.h"
#include "spline.h"
#include "voronoi.h"
#include "srg.h"

#include "stdlib.h"
#include "stdio.h"
#include "time.h"


void spline(SplineParams splineParams, ImageParams* imageParams) {
    Vec** s = getNSplines(splineParams.nSplines, splineParams.alpha, splineParams.minDist, splineParams.nPoints);
    printf("Created splines\n");
    createSplineImage(s, splineParams.nSplines, splineParams.nPoints, imageParams,
                      "/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/");
    printf("Created spline images!\n");
}

void voronoi(VoronoiParams voronoiParams, ImageParams* imageParams) {
//    Vec* particles = malloc(sizeof(Vec) * voronoiParams.nInitialCells);
//    Cell** cells = getCells(voronoiParams.nInitialCells, "test", &particles);

    Vec* particles = malloc(sizeof(Vec) * voronoiParams.nCells);
    Cell** cells = getCells(voronoiParams.nCells, "test", &particles);
    printf("Created cells\n");

//    printf("Particles:\n");
//    printVecs(particles, voronoiParams.nInitialCells);
//    printVecs(particles, voronoiParams.nCells);

    // TODO: when cells are merged, seeds need to be merged (choose one)
//    cells = mergeCells(cells, voronoiParams.nInitialCells, voronoiParams.nCells);
//    printf("Merged cells\n\n");
//    printCells(cells, voronoiParams.nCells);

    createVoronoiImage(cells, voronoiParams.nCells, imageParams,
                       "/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/");
    printf("Created voronoi images!\n");

//    discretizeVecs(particles, voronoiParams.nCells, imageParams.imageSize);
//    printf("Particles:\n");
//    printVecs(particles, voronoiParams.nCells);
//
//    Bitmap* bitmap = srg("/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/",
//                         particles, voronoiParams.nCells);
//    writeBitmap(bitmap, "/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/srg/");
}

int main() {
    srand(time(0));

    ImageParams imageParams = {
            .imageSize = 100,
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
    };

    SplineParams splineParams = {
            .alpha = 0.5,
            .minDist = 0.1,
            .nPoints = 100,
            .nSplines = 5,
    };

    clock_t begin;
    clock_t end;
    double time_spent;


    begin = clock();

//    spline(splineParams, &imageParams);

    voronoi(voronoiParams, &imageParams);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME: %lf\n", time_spent);

    // TODO: use initial points for cell creation as  seeds; and then one for cell walls
    // TODO: for image creation: dont calc dist from pixel to cells; calc cell
    // TODO: cgal

    return 0;
}
