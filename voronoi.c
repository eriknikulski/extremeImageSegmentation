//
// Created by Erik Nikulski on 19.04.22.
//

#include "voronoi.h"

#include "params.h"
#include "graph.h"
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
        particles[i].x = (double)(rand() + 0.05*RAND_MAX) / (double)(1.1*RAND_MAX);
        particles[i].y = (double)(rand() + 0.05*RAND_MAX) / (double)(1.1*RAND_MAX);
        particles[i].z = (double)(rand() + 0.05*RAND_MAX) / (double)(1.1*RAND_MAX);
    }

    return particles;
}

void createParticleFile(Vec* particles, int n, char* fname) {
    FILE* fp;
    fp = fopen(fname, "w");
    if (fp == NULL) exit(EXIT_FAILURE);

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%d %lf %lf %lf\n", i, particles[i].x, particles[i].y, particles[i].z);
    }

    fclose(fp);
}

void createCells(Vec* particles, int n, char* fname) {
    char* command = malloc(sizeof(char) * (45 + strlen(fname) + 1));
    strcpy(command, "voro++ -c \"%i:%w:%P:%s:%t:%l:%n\" 0 1 0 1 0 1 ");
    strcat(command, fname);
    strcat(command, "\0");

    createParticleFile(particles, n, fname);
    system(command);
    free(command);
}

void setParticleId(Cell* cell, char* s) {
    cell->id = (int)strtol(s, NULL, 10);
}

void setNodeCountStr(Cell* c, char* s) {
    c->nodeCount = (int)strtol(s, NULL, 10);
}

void setVerticesStr(Cell* c, char* s) {
    char* spaceDelimiter = " ";
    char* commaDelimiter = ",";
    char* elem;
    int iNodes = 0;

    char* split = NULL;
    char* saveptr1;
    char* saveptr2;

    for (split = strtok_r(s, spaceDelimiter, &saveptr1); split != NULL; split = strtok_r(NULL, spaceDelimiter, &saveptr1)) {
        // Remove brackets from vertices
        ++split;
        split[strlen(split) - 1] = 0;

        elem = strtok_r(split, commaDelimiter, &saveptr2);
        c->nodes[iNodes].x = strtod(elem, NULL);
        elem = strtok_r(NULL, commaDelimiter, &saveptr2);
        c->nodes[iNodes].y = strtod(elem, NULL);
        elem = strtok_r(NULL, commaDelimiter, &saveptr2);
        c->nodes[iNodes].z = strtod(elem, NULL);

        ++iNodes;
    }
}

int getNodeCountStr(char* s) {
    char* commaDelimiter = ",";
    char* saveptr;
    int i = 0;
    for (char* split = strtok_r(s, commaDelimiter, &saveptr); split != NULL; split = strtok_r(NULL, commaDelimiter, &saveptr)) {
        ++i;
    }
    return i;
}

void setFaceCountStr(Cell* c, char* s) {
    c->faceCount = (int)strtol(s, NULL, 10);
}

void setFaceNodesStr(Cell* c, Face* f, char* s) {
    char* commaDelimiter = ",";
    char* saveptr;
    int i = 0;

    for (char* split = strtok_r(s, commaDelimiter, &saveptr); split != NULL; split = strtok_r(NULL, commaDelimiter, &saveptr)) {
        f->nodes[i] = c->nodes[(int)strtol(split, NULL, 10)];
        ++i;
    }
}

void setFacesStr(Cell* c, char* s) {
    char* spaceDelimiter = " ";
    char* saveptr;
    int i = 0;

    for (char* split = strtok_r(s, spaceDelimiter, &saveptr); split != NULL; split = strtok_r(NULL, spaceDelimiter, &saveptr)) {
        // Remove brackets from vertices
        ++split;
        split[strlen(split) - 1] = 0;

        char* tmp = malloc(sizeof(char) * (strlen(split) + 1));
        memset(tmp, 0, sizeof(char) * (strlen(split) + 1));
        strcpy(tmp, split);
        c->faces[i].count = getNodeCountStr(tmp);
        free(tmp);

        c->faces[i].nodes = malloc(sizeof(Node) * c->faces[i].count);
        setFaceNodesStr(c, &c->faces[i], split);

        c->faces[i].faceCalcs = NULL;
        c->faces[i].nFaceCalcs = 0;

        ++i;
    }
}

void setNormalVecFace(Cell* c, char* s) {
    char* spaceDelimiter = " ";
    char* commaDelimiter = ",";
    int i = 0;
    char* elem;

    char* saveptr1;
    char* saveptr2;

    for (char* split = strtok_r(s, spaceDelimiter, &saveptr1); split != NULL; split = strtok_r(NULL, spaceDelimiter, &saveptr1)) {

        c->faces[i].normalVec = malloc(sizeof(Node));

        // Remove brackets from vertices
        ++split;
        split[strlen(split) - 1] = 0;

        elem = strtok_r(split, commaDelimiter, &saveptr2);
        c->faces[i].normalVec->x = strtod(elem, NULL);
        elem = strtok_r(NULL, commaDelimiter, &saveptr2);
        c->faces[i].normalVec->y = strtod(elem, NULL);
        elem = strtok_r(NULL, commaDelimiter, &saveptr2);
        c->faces[i].normalVec->z = strtod(elem, NULL);
        ++i;
    }
}

void setNeighboringCells(Cell* cell, char* s) {
    char* spaceDelimiter = " ";
    char* saveptr;
    int i = 0;

    cell->neighbors = malloc(sizeof(int) * cell->faceCount);
    cell->neighborCount = cell->faceCount;

    for (char* split = strtok_r(s, spaceDelimiter, &saveptr); split != NULL; split = strtok_r(NULL, spaceDelimiter, &saveptr)) {
        cell->neighbors[i] = (int)strtol(split, NULL, 10);
        ++i;
    }
}

Cell** readCells(int n, char* fname) {
    char* line = NULL;
    char* linetmp = NULL;
    size_t len = 0;

    FILE* fp;
    fp = fopen(fname, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

    Cell** cells = malloc(sizeof(void*) * n);

    char* colonDelimiter = ":";
    char* split;
    char* saveptr;

    for (int i = 0; getline(&linetmp, &len, fp) != -1; ++i) {
        line = realloc(line, sizeof(char) * strlen(linetmp) + 1);
        strcpy(line, linetmp);
        line[strlen(line) - 1] = 0;     // strip new line

        cells[i] = malloc(sizeof(Cell));

        split = strtok_r(line, colonDelimiter, &saveptr);
        setParticleId(cells[i], split);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setNodeCountStr(cells[i], split);
        cells[i]->nodes = malloc(sizeof(Node) * cells[i]->nodeCount);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setVerticesStr(cells[i], split);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setFaceCountStr(cells[i], split);
        cells[i]->faces = malloc(sizeof(Face) * cells[i]->faceCount);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setFacesStr(cells[i], split);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setNormalVecFace(cells[i], split);
        split = strtok_r(NULL, colonDelimiter, &saveptr);
        setNeighboringCells(cells[i], split);
    }

    free(line);
    return cells;
}

Cell** getCells(int n, char* fname, Vec** particles) {
    char* fout = malloc(sizeof(char) * (strlen(fname) + 5 + 1));
    strcpy(fout, fname);
    strcat(fout, ".vol\0");

    *particles = getRandParticles(n);
    createCells(*particles, n, fname);
    Cell** cells = readCells(n, fout);
    free(fout);

    for (int i = 0; i < n; ++i) {
        cells[i]->particle = &((*particles)[cells[i]->id]);
        cells[i]->nIds = 0;
    }

    return cells;
}

void removeDupNodes(Cell* c) {
    unsigned int count = c->nodeCount;
    for (int i = 0; i < count; ++i) {
        for (int j = i + 1; j < count; ++j) {
            if (equalVecs(&c->nodes[i], &c->nodes[j])) {
                if (j == count - 1) {
                    --count;
                    break;
                }
                c->nodes[j] = c->nodes[count - 1];
                --count;
                --j;
            }
        }
    }

    c->nodes = realloc(c->nodes, sizeof(Node) * count);
    c->nodeCount = count;
}

void removeDupFaces(Cell* c) {
    unsigned int count = c->faceCount;
    for (int i = 0; i < count - 1; ++i) {
        for (int j = i + 1; j < count; ++j) {
            if (equalFaces(&c->faces[i], &c->faces[j])) {
                if (i == count - 2 && j == count - 1) {
                    count -= 2;
                    break;
                }
                if (j != count - 1) {
                    c->faces[j] = c->faces[count - 1];
                }
                c->faces[i] = c->faces[count - 2];

                count -= 2;
                j = i + 1;
            }
        }
    }

    c->faces = realloc(c->faces, sizeof(Face) * count);
    c->faceCount = count;
}

Cell** actualizeAssignment(Cell** cs, int* seeds, int* assign, int nOld, int nNew) {
    Cell** cells = malloc(sizeof(Cell*) * nNew);

    for (int i = 0; i < nNew; ++i) {
        unsigned int nodeCount = 0;
        unsigned int faceCount = 0;
        int cNode = 0;
        int cFace = 0;
        int cCells = 0;
        int counter = 0;

        cells[i] = malloc(sizeof(Cell));
        Cell* cell = cells[i];

        int seed = seeds[i];
        cell->id = seed;
        cell->particle = NULL;

        for (int j = 0; j < nOld; ++j) {
            if (assign[j] != seed) continue;
            nodeCount += cs[j]->nodeCount;
            faceCount += cs[j]->faceCount;
            ++cCells;
        }

        cell->nIds = cCells;
        cell->ids = malloc(sizeof(int) * cell->nIds);

        cell->faces = malloc(sizeof(Face) * faceCount);
        cell->nodes = malloc(sizeof(Node) * nodeCount);
        cell->faceCount = faceCount;
        cell->nodeCount = nodeCount;

        for (int j = 0; j < nOld; ++j) {
            if (assign[j] != seed) continue;
            for (int k = 0; k < cs[j]->nodeCount; ++k) {
                cell->nodes[cNode] = cs[j]->nodes[k];
                ++cNode;
            }
            for (int k = 0; k < cs[j]->faceCount; ++k) {
                cell->faces[cFace] = cs[j]->faces[k];
                ++cFace;
            }
            cell->ids[counter] = cs[j]->id;
            ++counter;
        }

        removeDupNodes(cell);   // TODO: could be removed -> performance
        removeDupFaces(cell);
    }
    return cells;
}

Cell** mergeCells(Cell** cells, VoronoiParams* voronoiParams) {
    int nOld = voronoiParams->nInitialCells;
    int nNew = voronoiParams->nCells;
    assert(nOld >= nNew);
    if (nOld == nNew) return cells;

    Graph* graph = toGraph(cells, nOld);
    int* seeds = getNRandExc(nNew, 0, (int)graph->count - 1);
    int* assign = growSeeds(graph, seeds, nNew);

    return actualizeAssignment(cells, seeds, assign, nOld, nNew);
}

void discretizeFace(Face* face, int size) {
    discretizeVecs(face->nodes, (int)face->count, size);
}

void discretizeCell(Cell* cell, int size) {
    if (cell->particle)
        discretizeVec(cell->particle, size);
    discretizeVecs(cell->nodes, (int)cell->nodeCount, size);
    for (int i = 0; i < cell->faceCount; ++i) {
        discretizeFace(&cell->faces[i], size);
    }
}

void discretizeCells(Cell** cells, int n, int size) {
    for (int i = 0; i < n; ++i) {
        discretizeCell(cells[i], size);
    }
}

double getDistPlane(Vec* v, Face* f) {
    Vec vProj;
    Vec tmp;
    double d;

    subtractVec(v, &f->nodes[0], &tmp);
    d = getDotProd(f->normalVec, &tmp);

    multVecs(d, f->normalVec, &vProj);
    subtractVec(v, &vProj, &vProj);

    return getDist(&vProj, v);
}

double getDistLineSegSq(Vec* v1, Vec* v2, Vec* v) {
    // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    double lineDist = getDistSq(v1, v2);
    if (isAlmostZero(lineDist)) return getDistSq(v, v1);

    Vec tmpV;
    Vec tmp;
    double t;

    subtractVec(v, v1, &tmpV);
    subtractVec(v2, v1, &tmp);
    t = getDotProd(&tmpV, &tmp) / lineDist;

    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    multVecs(t, &tmp, &tmp);
    addVecs(v1, &tmp, &tmp);
    return getDistSq(v, &tmp);
}

double getDistFacePlane(Vec* v, Face* f) {
    // https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle/544947
    double alpha;
    double beta;
    double gamma;
    Vec w;
    Vec tmp;
    Vec res;

    FaceCalc* faceCalc;

    double cDist;
    double sDist = DBL_MAX;

    for (int i = 0; i < f->count - 2; ++i) {
        Vec* n1 = &f->nodes[0];
        Vec* n2 = &f->nodes[i + 1];
        Vec* n3 = &f->nodes[i + 2];

        if (f->nFaceCalcs < i + 1) {
            f->faceCalcs = realloc(f->faceCalcs, sizeof(FaceCalc) * (f->nFaceCalcs + 1));
            f->nFaceCalcs++;
            faceCalc = &f->faceCalcs[f->nFaceCalcs - 1];

            faceCalc->u = getSubVec(n2, n1);
            faceCalc->v = getSubVec(n3, n1);
            faceCalc->n = getCrossProduct(faceCalc->u, faceCalc->v);
            faceCalc->nn = getDotProd(faceCalc->n, faceCalc->n);
        }
        faceCalc = &f->faceCalcs[i];

        subtractVec(v, n1, &w);

        calcCrossProduct(faceCalc->u, &w, &tmp);
        gamma = getDotProd(&tmp, faceCalc->n) / faceCalc->nn;

        calcCrossProduct(&w, faceCalc->v, &tmp);
        beta = getDotProd(&tmp, faceCalc->n) / faceCalc->nn;

        alpha = 1 - gamma - beta;

        if (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1) {
            res.x = alpha * n1->x + beta * n2->x + gamma * n3->x;
            res.y = alpha * n1->y + beta * n2->y + gamma * n3->y;
            res.z = alpha * n1->z + beta * n2->z + gamma * n3->z;

            cDist = getDist(&res, v);
            if (cDist < sDist) sDist = cDist;
            if (isAlmostZero(sDist)) return 0;
        }
    }

    return sDist;
}

// TODO: improve performance
double getDistFaceLineSegs(Vec* v, Face* f) {
    double sDist = DBL_MAX;
    double cDist;

    // TODO: if nodes in order -> iterate over i, i+1
    for (int i = 0; i < f->count - 1; ++i) {
        for (int j = i + 1; j < f->count; ++j) {
            cDist = getDistLineSegSq(&f->nodes[i], &f->nodes[j], v);
            if (cDist < sDist) sDist = cDist;
            if (isAlmostZero(sDist)) return 0;
        }
    }

    return sDist == DBL_MAX ? DBL_MAX : sqrt(sDist);
}

double getDistFace(Vec* v, Face* f) {
    double plane = getDistFacePlane(v, f);
    double lineSeg = getDistFaceLineSegs(v, f);
    return plane < lineSeg ? plane : lineSeg;
}

double getDistCell(Vec* v, Cell* cell) {
    double sDist = DBL_MAX;
    double cDist;

    for (int i = 0; i < cell->faceCount; ++i) {
        cDist = getDistFace(v, &cell->faces[i]);
        if (cDist < sDist) sDist = cDist;
        if (isAlmostZero(sDist)) return 0;
    }

    return sDist;
}

double getDistCells(Vec* v, Cell** cells, int n) {
    double sDist = DBL_MAX;
    double cDist;

    for (int i = 0; i < n; ++i) {
        cDist = getDistCell(v, cells[i]);
        if (cDist < sDist) sDist = cDist;
    }

    return sDist;
}

int equalFaces(Face* f1, Face* f2) {
    if (f1->count != f2->count) return 0;
    if (!equalVecs(f1->normalVec, f2->normalVec) && !negEqualVecs(f1->normalVec, f2->normalVec)) return 0;
    int common = 0;
    for (int i = 0; i < f1->count; ++i) {
        for (int j = 0; j < f2->count; ++j) {
            if (equalVecs(&f1->nodes[i], &f2->nodes[j])) ++common;
        }
    }
    return common == f1->count ? 1 : 0;
}

Node* getNodeCopy(Node* node) {
    Node* res = malloc(sizeof(Node));
    res->x = node->x;
    res->y = node->y;
    res->z = node->z;
    return res;
}

Node* getNodesCopy(Node* nodes, int count) {
    Node* res = malloc(sizeof(Node) * count);
    for (int i = 0; i < count; ++i) {
        res[i].x = nodes[i].x;
        res[i].y = nodes[i].y;
        res[i].z = nodes[i].z;
    }
    return res;
}

Face* getFaceCopy(Face* f) {
    Face* res = malloc(sizeof(Face));
    res->nodes = getNodesCopy(f->nodes, f->count);
    res->count = f->count;
    res->normalVec = getNodeCopy(f->normalVec);
    return res;
}

void printCells(Cell** cells, int count) {
    for (int i = 0; i < count; ++i) {
        printf("Cell %d ID=%d  nodeCount=%d  faceCount=%d  neighborCount=%d\n", i, cells[i]->id, cells[i]->nodeCount,
               cells[i]->faceCount, cells[i]->neighborCount);
        if (cells[i]->particle)
            printf("Particle:  x=%lf  y=%lf  z=%lf\n", cells[i]->particle->x, cells[i]->particle->y, cells[i]->particle->z);
        if (cells[i]->nIds) {
            printf("Cell IDs: ");
            for (int j = 0; j < cells[i]->nIds; ++j) {
                printf(" %d", cells[i]->ids[j]);
            }
            printf("\n");
        }
        for (int j = 0; j < cells[i]->nodeCount; ++j)
            printf("   Node %d  x=%lf  y=%lf  z=%lf\n", j, cells[i]->nodes[j].x, cells[i]->nodes[j].y, cells[i]->nodes[j].z);
        for (int j = 0; j < cells[i]->faceCount; ++j) {
            printf("   Face %d  count=%d  normVec={ x=%lf  y=%lf  z=%lf }\n", j, cells[i]->faces[j].count,
                   cells[i]->faces[j].normalVec->x, cells[i]->faces[j].normalVec->y, cells[i]->faces[j].normalVec->z);
            for (int k = 0; k < cells[i]->faces[j].count; ++k)
                printf("      Node %d  x=%lf  y=%lf  z=%lf\n", k,
                       cells[i]->faces[j].nodes[k].x, cells[i]->faces[j].nodes[k].y, cells[i]->faces[j].nodes[k].z);
        }
        for (int j = 0; j < cells[i]->neighborCount; ++j) {
            printf("   Neighbor %d  ID=%d\n", j, cells[i]->neighbors[j]);
        }
        printf("\n\n");
    }
}
