#include "spline.h"

int main() {
    int nPoints = 100;
    Vec p0 = {.x = 0,
              .y = 1.5,
              .z = 1};
    Vec p1 = {.x = 1,
              .y = 2,
              .z = 1};
    Vec p2 = {.x = 3,
              .y = 2,
              .z = 1};
    Vec p3 = {.x = 4,
              .y = 2,
              .z = 1.5};

    Vec* result = getCatmullRomSpline(&p0, &p1, &p2, &p3, nPoints);
    printVecs(result, nPoints);

    return 0;
}
