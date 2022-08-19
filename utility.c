//
// Created by Erik Nikulski on 11.04.22.
//

#include "assert.h"
#include "utility.h"

#include "math.h"
#include "stdlib.h"

double getRandD() {
    return (double)rand() / (double)(RAND_MAX);
}

int getRand(int lower, int upper) {
    return (rand() % (upper - lower + 1)) + lower;
}

int* getRandPair(int lower, int upper) {
    int r1 = (rand() % (upper - lower + 1)) + lower;
    int r2 = r1;
    while (r2 == r1) {
        r2 = (rand() % (upper - lower + 1)) + lower;
    }

    int* res = malloc(sizeof(void*) * 2);
    res[0] = r1;
    res[1] = r2;
    return res;
}

int* getNRandExc(int n, int lower, int upper) {
    assert(n > 0);
    assert(lower < upper);
    assert(n <= (upper - lower + 1));

    int* res = malloc(sizeof(int) * n);
    int a[upper - lower + 1];
    int tmp;
    int j;

    // fill array
    for (int i = 0; i < upper - lower + 1; ++i) {
        a[i] = lower + i;
    }

    // shuffle
    for (int i = 0; i < n; ++i) {
        j = getRand(i, upper - lower);
        tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }

    // fill results
    for (int i = 0; i < n; ++i) {
        res[i] = a[i];
    }

    return res;
}

int getRandFrom(int* l, int n) {
    assert(n > 0);
    return l[getRand(0, n - 1)];
}


double normDist(double sigma, double mu) {
    double randomNum_normal = sqrt(-2 * log((rand() + 1.0) / (RAND_MAX + 1.0)))
            * cos(2 * M_PI * (rand() + 1.0) / (RAND_MAX + 1.0));
    return mu + sigma * randomNum_normal;
}

double normalCDF(double value, double sigma, double mu) {
//    return 0.5 * erfc((value - mu) / (sigma * M_SQRT2));
    return 0.5 * erfc((mu - value) / (sigma * M_SQRT2));
}
