//
// Created by Erik Nikulski on 19.06.22.
//

#ifndef EXTREMEIMAGESEGMENTATION_PNG_H
#define EXTREMEIMAGESEGMENTATION_PNG_H

#include "image.h"

int save_bitmap_slice_to_file(Bitmap* bitmap, int z, const char *path);

Bitmap* read_pngs(const char* path);

#endif //EXTREMEIMAGESEGMENTATION_PNG_H
