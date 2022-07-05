//
// Created by Erik Nikulski on 01.07.22.
//

#include "srg.h"

#include "png.h"
#include "vec.h"
#include "image.h"

#include "math.h"
#include "stdlib.h"
#include "stdio.h"


double calcDelta(Seed* s, uint8_t value) {
    return fabs(value - s->sum / s->count);
}

SSL* insertSSL(SSL* head, Pixel* p, double delta) {
    p->inSSL = 1;
    if (head == NULL) {
        head = malloc(sizeof(SSL));
        head->p = p;
        head->delta = delta;
        head->next = NULL;
        return head;
    }

    SSL* current = head;
    SSL* prev = NULL;
    SSL* elem;
    elem = malloc(sizeof(SSL));
    elem->p = p;
    elem->delta = delta;
    elem->next = NULL;

    while (current->next) {
        if (current->delta > delta) {
            elem->next = current;
            if (prev) {
                prev->next = elem;
            }
            return head;
        }
        prev = current;
        current = current->next;
    }

    current->next = elem;
    return head;
}

Bitmap* srg(char* path, Vec* vSeeds, int nSeeds) {
    // TODO: use border label
    int neighborCount = 0;
    Pixel** neighbors;

    Seed* seeds = malloc(sizeof(Seed) * nSeeds);

    double delta;
    SSL* head = NULL;
    SSL* current = NULL;
    Seed* label = NULL;

    Bitmap* bitmap = read_pngs(path);
    // Label seed points according their initial grouping.
    for (int i = 0; i < nSeeds; ++i) {
        Pixel* p = getPixel(bitmap, (int)vSeeds[i].x, (int)vSeeds[i].y, (int)vSeeds[i].z);
        seeds[i].id = i;
        seeds[i].count = 1;
        seeds[i].sum = p->value;
        seeds[i].p = p;
        p->grouping = &seeds[i];
        p->inSSL = 1;
    }


    int sslCount = 0;


    // Put neighbors of seed points (the initial T ) in the SSL.
    for (int i = 0; i < nSeeds; ++i) {
        neighbors = getNeighbors(bitmap, seeds[i].p, &neighborCount);
        for (int j = 0; j < neighborCount; ++j) {
            delta = calcDelta(&seeds[i], neighbors[j]->value);
            head = insertSSL(head, neighbors[j], delta);
        }
        sslCount += neighborCount;
        free(neighbors);
    }

    while (head) {
        label = NULL;
        current = head;
        head = head->next;

        neighbors = getNeighbors(bitmap, current->p, &neighborCount);
        for (int i = 0; i < neighborCount; ++i) {
            if (neighbors[i]->grouping) {
                if (label == NULL) {
                    label = neighbors[i]->grouping;
                }
                if (label != neighbors[i]->grouping) {
                    label = NULL;
                    break;
                }
            }
        }

        if (!label)
            continue;

        current->p->grouping = label;
        ++(label->count);
        label->sum += current->p->value;

        for (int i = 0; i < neighborCount; ++i) {
            if (neighbors[i]->inSSL || neighbors[i]->grouping) continue;
            delta = calcDelta(label, neighbors[i]->value);
            insertSSL(head, neighbors[i], delta);
            sslCount++;
        }
        free(neighbors);
    }

    printf("SSL COUNT: %d\n", sslCount);

    setGroupings(bitmap);

    return bitmap;
}

void setGroupings(Bitmap* bitmap) {
    Pixel* pixel;
    for (int z = 0; z < bitmap->size; ++z) {
        for (int y = 0; y < bitmap->size; ++y) {
            for (int x = 0; x < bitmap->size; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                if (pixel->grouping) {
                    pixel->value = pixel->grouping->p->value;
                } else {
                    pixel->value = 0;
                }
            }
        }
    }
}