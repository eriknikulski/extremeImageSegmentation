//
// Created by Erik Nikulski on 19.06.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_PNG_H
#define EXTREMEIMAGESEGMENTATION_PNG_H

#include "libpng-1.6.37/png.h"
#include "stdint.h"

typedef struct {
    uint8_t value;
} pixel_t;

typedef struct {
    pixel_t *pixels;
    size_t width;
    size_t height;
} bitmap_t;

pixel_t* pixel_at(bitmap_t* bitmap, int x, int y);

int save_png_to_file(bitmap_t *bitmap, const char *path);

#endif //EXTREMEIMAGESEGMENTATION_PNG_H
