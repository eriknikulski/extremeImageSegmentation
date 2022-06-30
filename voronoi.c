//
// Created by Erik Nikulski on 19.04.22.
//

#include "voronoi.h"

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
        particles[i].x = getRandD();
        particles[i].y = getRandD();
        particles[i].z = getRandD();
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
    char* command = malloc(sizeof(char) * (45 + strlen(fname) + 1));
    strcpy(command, "voro++ -c \"%i:%w:%P:%s:%t:%l:%n\" 0 1 0 1 0 1 ");
    strcat(command, fname);
    strcat(command, "\0");

    createParticleFile(n, fname);
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

        c->faces[i].nFaceCalcs = 0;
        c->faces[i].faceCalcs = NULL;
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
        line = realloc(line, sizeof(char) * strlen(linetmp));
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

Cell** getCells(int n, char* fname) {
    char* fout = malloc(sizeof(char) * (strlen(fname) + 5 + 1));
    strcpy(fout, fname);
    strcat(fout, ".vol\0");

    createCells(n, fname);
    Cell** cells = readCells(n, fout);
    free(fout);
    return cells;
}

FaceCalc* getFaceCalc(Face* f, int i) {
    for (int j = 0; j < f->nFaceCalcs; ++j) {
        if (f->faceCalcs[j].n1 == &f->nodes[i] &&
            f->faceCalcs[j].n2 == &f->nodes[i + 1] &&
            f->faceCalcs[j].n3 == &f->nodes[i + 2])
            return &f->faceCalcs[j];
    }
    return NULL;
}

double getDistFaceNode(Face* f, Node* p) {
    // TODO: use hashmap or something for nodes
    // https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle/544947
    double alpha;
    double beta;
    double gamma;
    Vec* w;
    Vec* tmp;
    Vec res;

    double cDist;
    double sDist = MAXFLOAT;

    for (int i = 0; i < f->count - 2; ++i) {
        FaceCalc* faceCalc = getFaceCalc(f, i);
        if (!faceCalc) {
            f->faceCalcs = realloc(f->faceCalcs, sizeof(FaceCalc) * (f->nFaceCalcs + 1));
            f->nFaceCalcs++;
            faceCalc = &f->faceCalcs[f->nFaceCalcs - 1];
            faceCalc->n1 = &f->nodes[i];
            faceCalc->n2 = &f->nodes[i + 1];
            faceCalc->n3 = &f->nodes[i + 2];

            faceCalc->u = getSubVec(&f->nodes[i + 1], &f->nodes[i]);
            faceCalc->v = getSubVec(&f->nodes[i + 2], &f->nodes[i]);
            faceCalc->n = getCrossProduct(faceCalc->u, faceCalc->v);
            faceCalc->nn = getDotProd(faceCalc->n, faceCalc->n);
        }

        w = getSubVec(p, &f->nodes[i]);

        tmp = getCrossProduct(faceCalc->u, w);
        gamma = getDotProd(tmp, faceCalc->n) / faceCalc->nn;

        tmp = getCrossProduct(w, faceCalc->v);
        beta = getDotProd(tmp, faceCalc->n) / faceCalc->nn;

        alpha = 1 - gamma - beta;

        if (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1) {
            res.x = alpha * f->nodes[i].x + beta * f->nodes[i + 1].x + gamma * f->nodes[i + 2].x;
            res.y = alpha * f->nodes[i].y + beta * f->nodes[i + 1].y + gamma * f->nodes[i + 2].y;
            res.z = alpha * f->nodes[i].z + beta * f->nodes[i + 1].z + gamma * f->nodes[i + 2].z;

            cDist = getDist(&res, p);
            if (cDist < sDist) sDist = cDist;
            if (isZero(sDist)) return 0;
        }
    }
    return sDist;
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
    for (int i = 0; i < count; ++i) {
        for (int j = i + 1; j < count; ++j) {
            if (equalFaces(&c->faces[i], &c->faces[j])) {
                if (j == count - 1) {
                    --count;
                    break;
                }
                c->faces[j] = c->faces[count - 1];
                --count;
                --j;
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

        cells[i] = malloc(sizeof(Cell));
        Cell* cell = cells[i];

        int seed = seeds[i];
        cell->id = seed;

        for (int j = 0; j < nOld; ++j) {
            if (assign[j] != seed) continue;
            nodeCount += cs[j]->nodeCount;
            faceCount += cs[j]->faceCount;
        }

        cell->faces = malloc(sizeof(Face) * faceCount);
        cell->nodes = malloc(sizeof(Node) * nodeCount);
        cell->neighbors = malloc(sizeof(int) * nNew);
        cell->faceCount = faceCount;
        cell->nodeCount = nodeCount;
        cell->neighborCount = nNew;
        int actualNeighborCount = 0;

        for (int j = 0; j < nOld; ++j) {
            if (assign[j] != seed) continue;
            for (int k = 0; k < cs[j]->nodeCount; ++k) {
                cell->nodes[cNode] = cs[j]->nodes[k];
                ++cNode;
            }
            for (int k = 0; k < cs[j]->faceCount; ++k) {
                cell->faces[cFace] = cs[j]->faces[k];
                ++cFace;

                if (cs[j]->neighbors[k] < 0 || assign[cs[j]->neighbors[k]] == seed) continue;
                int hasNeighbor = 0;
                for (int l = 0; l < actualNeighborCount; ++l) {
                    if (assign[cs[j]->neighbors[k]] == cell->neighbors[l]) {
                        hasNeighbor = 1;
                        break;
                    }
                }
                if (!hasNeighbor) {
                    cell->neighbors[actualNeighborCount] = assign[cs[j]->neighbors[k]];
                    ++actualNeighborCount;
                }
            }
        }

        cell->neighbors = realloc(cell->neighbors, sizeof(int) * actualNeighborCount);
        cell->neighborCount = actualNeighborCount;

        removeDupNodes(cell);
        removeDupFaces(cell);
    }
    return cells;
}

Cell** mergeCells(Cell** cells, int nOld, int nNew) {
    assert(nOld >= nNew);
    if (nOld == nNew) return cells;

    Graph* graph = toGraph(cells, nOld);
    int* seeds = getNRandExc(nNew, 0, (int)graph->count - 1);
    int* assign = growSeeds(graph, seeds, nNew);

    // merge cells based on assignments
    return actualizeAssignment(cells, seeds, assign, nOld, nNew);
}

void discretizeFace(Face* face, int size) {
    for (int i = 0; i < face->count; ++i) {
        face->nodes[i].x = round(face->nodes[i].x * (double)size);
        face->nodes[i].y = round(face->nodes[i].y * (double)size);
        face->nodes[i].z = round(face->nodes[i].z * (double)size);
    }
}

void discretizeCell(Cell* cell, int size) {
    for (int i = 0; i < cell->nodeCount; ++i) {
        cell->nodes[i].x = round(cell->nodes[i].x * (double)size);
        cell->nodes[i].y = round(cell->nodes[i].y * (double)size);
        cell->nodes[i].z = round(cell->nodes[i].z * (double)size);
    }
    for (int i = 0; i < cell->faceCount; ++i) {
        discretizeFace(&cell->faces[i], size);
    }
}

void discretizeCells(Cell** cells, int n, int size) {
    for (int i = 0; i < n; ++i) {
        discretizeCell(cells[i], size);
    }
}

double getDistLineSeg(Vec* v1, Vec* v2, Vec* v) {
    // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    double lineDist = getDist(v1, v2);
    if (isZero(lineDist)) return getDist(v1, v);
    Vec* tmp = getSubVec(v2, v1);
    Vec* tmpV = getSubVec(v, v1);
    double t = getDotProd(tmp, tmpV) / lineDist;

    if (t <= 0.0) {
        t = 0.0;
    } else if (t >= 1.0) {
        t = 1.0;
    }

    multVecs(t, tmp, tmp);
    addVecs(v1, tmp, tmp);
    double dist = getDist(v, tmp);

    free(tmp);
    free(tmpV);
    return dist;
}

double getDistFace(Face* f, Vec* v) {
    // https://stackoverflow.com/questions/10983872/distance-from-a-point-to-a-polygon

    // check if vec on face plane
    // if on plane -> calc if inside face
    // else:

    double sDist = DBL_MAX;
    double cDist;

    double d = getDistFaceNode(f, v);
    if (d < MAXFLOAT) return d;

    for (int i = 0; i < f->count - 1; ++i) {
        cDist = getDistLineSeg(&f->nodes[i], &f->nodes[i + 1], v);
        if (cDist < sDist) sDist = cDist;
        if (isZero(sDist)) return 0.0;   // TODO: should in theory never occur; check and remove
    }
    return sDist;
}

double voronoiDist(Vec* v, Cell** cells, int nCells) {
    double sDist = DBL_MAX;
    double cDist;

    for (int i = 0; i < nCells; ++i) {
        for (int j = 0; j < cells[i]->faceCount; ++j) {
            cDist = getDistFace(&cells[i]->faces[j], v);
            if (cDist < sDist) sDist = cDist;
            if (isZero(sDist)) return 0.0;
        }
    }
    return sDist;
}

int equalFaces(Face* f1, Face* f2) {
    if (f1->normalVec != f2->normalVec || f1->count != f2->count) return 0;
    int common = 0;
    for (int i = 0; i < f1->count; ++i) {
        for (int j = 0; j < f2->count; ++j) {
            if (equalVecs(&f1->nodes[i], &f2->nodes[j])) ++common;
        }
    }
    return common == f1->count ? 1 : 0;
}

void printCells(Cell** cells, int count) {
    for (int i = 0; i < count; ++i) {
        printf("Cell %d ID=%d  nodeCount=%d  faceCount=%d  neighborCount=%d\n", i, cells[i]->id, cells[i]->nodeCount,
               cells[i]->faceCount, cells[i]->neighborCount);
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
