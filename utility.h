//
// Created by Erik Nikulski on 11.04.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_UNTIL_H
#define EXTREMEIMAGESEGMENTATION_UNTIL_H

double getRandD();

int getRand(int lower, int upper);

int* getRandPair(int lower, int upper);

int* getNRandExc(int n, int lower, int upper);

int getRandFrom(int* l, int n);

double normDist(double sigma, double mu);

double normalCDF(double value, double sigma, double mu);

#endif //EXTREMEIMAGESEGMENTATION_UNTIL_H
