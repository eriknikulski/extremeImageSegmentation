//
// Created by Erik Nikulski on 11.04.22.
//

#include "math.h"
#include "stdlib.h"

double getRand() {
    return (double)rand() / (double)(RAND_MAX);
}

double normDist(double sigma, double mu) {
    double randomNum_normal = sqrt(-2 * log((rand() + 1.0) / (RAND_MAX + 1.0)))
            * cos(2 * M_PI * (rand() + 1.0) / (RAND_MAX + 1.0));
    return mu + sigma * randomNum_normal;
}
