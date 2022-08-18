//
// Created by Erik Nikulski on 01.07.22.
//

#include "srg.h"

#include "params.h"
#include "png.h"
#include "vec.h"
#include "image.h"

#include "math.h"
#include "stdlib.h"
#include "stdio.h"


double calcDelta(Seed* s, uint8_t value, int precision) {
    return round(fabs(value - s->sum / s->count) * precision) / precision;
}

SSL* insertSSL(SSL* head, Pixel* p, double delta) {
    p->inSSL = 1;

    SSL* current = head;
    SSL* prev = NULL;
    SSL* elem = malloc(sizeof(SSL));
    elem->p = p;
    elem->delta = delta;
    elem->next = NULL;
    elem->nextHigher = NULL;

    if (head == NULL)
        return elem;

    while (current) {
        if (delta == current->delta) {
            elem->next = current->next;
            current->next = elem;
            return head;
        }
        if (delta < current->delta) {
            elem->nextHigher = current;
            elem->next = current;
            if (prev) {
                prev->nextHigher = elem;
                while (prev->next != current) {
                    prev = prev->next;
                }
                prev->next = elem;
                return head;
            }
            return elem;
        }
        prev = current;
        current = current->nextHigher;
    }

    current = prev;
    current->nextHigher= elem;
    while (current->next) {
        current = current->next;
    }
    current->next = elem;
    return head;
}

Bitmap* srg(Bitmap* bitmap, Vec* vSeeds, int nSeeds, int precision, char* imagePath, ImageParams* imageParams) {
    int neighborCount = 0;
    Pixel** neighbors;

    Seed* seeds = malloc(sizeof(Seed) * nSeeds);

    double delta;
    SSL* head = NULL;
    SSL* current = NULL;
    Seed* label = NULL;

    // Label seed points according their initial grouping.
    for (int i = 0; i < nSeeds; ++i) {
        Pixel* p = getPixel(bitmap, (int)vSeeds[i].x, (int)vSeeds[i].y, (int)vSeeds[i].z);
        seeds[i].id = i;
        seeds[i].count = 1;
        seeds[i].sum = p->value;
        seeds[i].p = p;
        p->particle = &vSeeds[i];
        p->grouping = &seeds[i];
        p->inSSL = 1;
    }

    // Put neighbors of seed points (the initial T ) in the SSL.
    for (int i = 0; i < nSeeds; ++i) {
        neighbors = getNeighbors(bitmap, seeds[i].p, &neighborCount);
        for (int j = 0; j < neighborCount; ++j) {
            delta = calcDelta(&seeds[i], neighbors[j]->value, precision);
            head = insertSSL(head, neighbors[j], delta);
        }
        free(neighbors);
    }

    while (head) {
        label = NULL;
        current = head;
        head = head->next;

        if (head && head->delta == current->delta)
            head->nextHigher = current->nextHigher;

        neighbors = getNeighbors(bitmap, current->p, &neighborCount);
        for (int i = 0; i < neighborCount; ++i) {
            if (neighbors[i]->grouping) {
                if (label == NULL) {
                    label = neighbors[i]->grouping;
                    continue;
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
        current->p->particle = label->p->particle;
        ++(label->count);
        label->sum += current->p->value;

        for (int i = 0; i < neighborCount; ++i) {
            if (neighbors[i]->inSSL) continue;
            delta = calcDelta(label, neighbors[i]->value, precision);
            head = insertSSL(head, neighbors[i], delta);
        }
        free(neighbors);
        free(current);
    }

    printf("Setting grouping values\n");
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