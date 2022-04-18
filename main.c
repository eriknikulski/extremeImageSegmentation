#include "spline.h"

#include "image.h"

#include "stdlib.h"
#include "stdio.h"
#include "time.h"

int main() {
    srand(time(0));

    int imageSize = 100;

    double alpha = 0.5;
    double minDist = 0;
    int nPoints = 100;

    int nSplines = 1;

    double theta_0 = 5.0;
    double theta_1 = -0.1;
    double sigma_b = 0.3;
    double mu_b = 0.3;
    double sigma_c = 0.1;
    double mu_c = 0.9;


    clock_t begin = clock();

    Vec** s = getNSplines(nSplines, alpha, minDist, nPoints);
    createImage(s, nSplines, nPoints, imageSize,
                theta_0, theta_1, sigma_b, mu_b, sigma_c, mu_c, "image");

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("TIME: %lf", time_spent);
    return 0;
}
