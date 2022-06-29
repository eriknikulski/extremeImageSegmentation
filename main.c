#include "image.h"
#include "spline.h"
#include "voronoi.h"

#include "stdlib.h"
#include "stdio.h"
#include "time.h"

int main() {
    srand(time(0));

    int imageSize = 100;

    double alpha = 0.5;
    double minDist = 0;
    int nPoints = 100;

    int nSplines = 5;

    ImageParams imageParams = {
            .imageSize = 100,
            .theta_0 = -3.0,
            .theta_1 = 1,
            .sigma_b = 5,
            .mu_b = 3,
            .sigma_c = 5,
            .mu_c = -3,
    };

    int mCells = 4;
    int nCells = 4 * mCells;

    clock_t begin;
    clock_t end;
    double time_spent;


    begin = clock();

    Vec** s = getNSplines(nSplines, alpha, minDist, nPoints);
    printf("Created splines\n");
    createSplineImage(s, nSplines, nPoints, imageSize, imageParams,
                      "/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/spline/");
    printf("Created spline images!\n");

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME splines: %lf", time_spent);


    begin = clock();

    Cell** cells = getCells(nCells, "test");
    printf("Created cells\n");

    cells = mergeCells(cells, nCells, mCells);
    printf("Merged cells\n\n");

    printCells(cells, mCells);

    createVoronoiImage(cells, mCells, imageSize, imageParams,
                       "/Users/eriknikulski/CLionProjects/extremeImageSegmentation/images/voronoi/");
    printf("Created voronoi images!\n");

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME voronoi: %lf", time_spent);


    return 0;
}
