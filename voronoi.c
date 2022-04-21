//
// Created by Erik Nikulski on 19.04.22.
//

#include "voronoi.h"

#include "utility.h"

#include "assert.h"
#include "float.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

Vec* getRandParticles(int n) {
    Vec* particles = malloc(sizeof(Vec) * n);

    for (int i = 0; i < n; ++i) {
        particles[i].x = getRand();
        particles[i].y = getRand();
        particles[i].z = getRand();
    }

    return particles;
}

void createParticleFile(int n, char* fname) {
    Vec* particles = getRandParticles(n);

    FILE* fp;
    fp = fopen(fname, "w");
    if (fp == NULL) exit(EXIT_FAILURE);

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%d %lf %lf %lf\n", i, particles[i].x, particles[i].y, particles[i].z);
    }

    fclose(fp);
    free(particles);
}

void createCells(int n, char* fname) {
    char* command = malloc(sizeof(char) * (30 + strlen(fname) + 1));
    strcpy(command, "voro++ -c \"%w %P\" 0 1 0 1 0 1 ");
    strcat(command, fname);
    strcat(command, "\0");

    createParticleFile(n, fname);
    system(command);
    free(command);
}

Cell** readCells(int n, char* fname) {
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    FILE* fp;
    fp = fopen(fname, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

    char* spaceDelimiter = " ";
    char* commaDelimiter = ",";
    int isNodeCount;
    int nodeCount = -1;

    Cell** cells = malloc(sizeof(void*) * n);
    int iCells = 0;
    int iNodes = 0;

    char* elem;

    while ((read = getline(&line, &len, fp)) != -1) {
        isNodeCount = 1;
        line[strlen(line) - 1] = 0;

        for (char* split = strtok_r(line, spaceDelimiter, &line); split != NULL; split = strtok_r(line, spaceDelimiter, &line)) {
            if (isNodeCount) {
                nodeCount = (int)strtol(split, NULL, 10);
                isNodeCount = 0;

                cells[iCells] = malloc(sizeof(Cell));
                cells[iCells]->count = nodeCount;
                cells[iCells]->nodes = malloc(sizeof(Node) * nodeCount);
                continue;
            }

            split++;
            split[strlen(split) - 1] = 0;

            elem = strtok_r(split, commaDelimiter, &split);
            cells[iCells]->nodes[iNodes].x = strtod(elem, NULL);
            elem = strtok_r(split, commaDelimiter, &split);
            cells[iCells]->nodes[iNodes].y = strtod(elem, NULL);
            elem = strtok_r(split, commaDelimiter, &split);
            cells[iCells]->nodes[iNodes].z = strtod(elem, NULL);

            iNodes++;
        }
        iCells++;
        iNodes = 0;
    }

    free(line);
    return cells;
}

Cell** getCells(int n, char* fname) {
    char* fout = malloc(sizeof(char) * (strlen(fname) + 4 + 1));
    strcpy(fout, fname);
    strcat(fout, ".vol\0");

    createCells(n, fname);
    Cell** cells = readCells(n, fout);
    free(fout);
    return cells;
}

void discretizeCell(Cell* cell, int size) {
    assert(size % 10 == 0);
    for (int i = 0; i < cell->count; ++i) {
        cell->nodes[i].x = round(cell->nodes[i].x * (double)size);
        cell->nodes[i].y = round(cell->nodes[i].y * (double)size);
        cell->nodes[i].z = round(cell->nodes[i].z * (double)size);
    }
}

void discretizeCells(Cell** cells, int n, int size) {
    assert(size % 10 == 0);
    for (int i = 0; i < n; ++i) {
        discretizeCell(cells[i], size);
    }
}

double voronoiDist(Vec* v, Cell** cells, int nCells) {
    double sDist = DBL_MAX;
    double cDist;
    for (int i = 0; i < nCells; ++i) {
        for (int j = 0; j < cells[i]->count - 1; ++j) {
            cDist = getDistLine(&cells[i]->nodes[j], &cells[i]->nodes[j + 1], v);
            if (cDist < sDist) sDist = cDist;
            if (sDist == 0.0) return 0.0;
        }
    }
    return sDist;
}
