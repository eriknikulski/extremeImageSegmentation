//
// Created by Erik Nikulski on 19.06.22.
//

#include "png.h"

#include "dirent.h"
#include "libpng-1.6.37/png.h"
#include "stdio.h"
#include "stdint.h"
#include "errno.h"
#include "string.h"
#include "stdlib.h"
#include "stdarg.h"


int save_bitmap_slice_to_file(Bitmap* bitmap, int z, const char *path) {
    FILE* fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_byte** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int pixel_size = 1;
    int depth = 8;

    fp = fopen(path, "wb");
    if (!fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }

    /* Set up error handling. */

    if (setjmp(png_jmpbuf(png_ptr))) {
        goto png_failure;
    }

    /* Set image attributes. */

    png_set_IHDR(png_ptr,
                 info_ptr,
                 bitmap->size,
                 bitmap->size,
                 depth,
                 PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    /* Initialize rows of PNG. */

    row_pointers = png_malloc(png_ptr, bitmap->size * sizeof(png_byte*));
    for (int y = 0; y < bitmap->size; ++y) {
        png_byte* row = png_malloc(png_ptr, sizeof(uint8_t) * bitmap->size * pixel_size);
        row_pointers[y] = row;
        for (int x = 0; x < bitmap->size; ++x) {
            Pixel* pixel = getPixel(bitmap, x, y, z);
            *row++ = pixel->value;
        }
    }

    /* Write the image data to "fp". */

    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */

    status = 0;

    for (int y = 0; y < bitmap->size; ++y) {
        png_free(png_ptr, row_pointers[y]);
    }
    png_free(png_ptr, row_pointers);

    png_failure:
    png_create_info_struct_failed:
    png_destroy_write_struct(&png_ptr, &info_ptr);
    png_create_write_struct_failed:
    fclose(fp);
    fopen_failed:
    return status;
}

static void fatal_error (const char * message, ...) {
    va_list args;
    va_start(args, message);
    vfprintf (stderr, message, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

int get_file_count(const char* dir_path) {
    int counter = 0;
    struct dirent* entry;
    DIR* dir = opendir(dir_path);

    if(dir == NULL) {
        printf("Error! Unable to read directory");
        return 0;
    }

    while((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_REG)
            ++counter;
    }

    closedir(dir);
    return counter;
}

Bitmap* read_pngs(const char* path, int size) {
    png_structp	png_ptr;
    png_infop info_ptr;
    FILE * fp;
    png_uint_32 width;
    png_uint_32 height;
    int bit_depth;
    int color_type;
    int interlace_method;
    int compression_method;
    int filter_method;
    png_bytepp rows;

    Bitmap* bitmap = malloc(sizeof(Bitmap));
    Pixel* pixel;

    // assumes image is a cube
    bitmap->size = size;
    bitmap->size2 = size * size;
    bitmap->pixels = calloc(size * size * size, sizeof(Pixel));

    for (int z = 0; z < size; ++z) {
        size_t len = snprintf(NULL, 0, "%simg_%d.png", path, z) + 1;
        char* activeName = malloc(len);
        snprintf(activeName, len, "%simg_%d.png", path, z);

        fp = fopen(activeName, "rb");
        free(activeName);

        if (!fp)
            fatal_error("Cannot open '%s': %s\n", path, strerror(errno));

        png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png_ptr)
            fatal_error("Cannot create PNG read structure");

        info_ptr = png_create_info_struct(png_ptr);
        if (!png_ptr)
            fatal_error("Cannot create PNG info structure");

        png_init_io(png_ptr, fp);
        png_read_png(png_ptr, info_ptr, 0, 0);
        png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
                     &interlace_method, &compression_method, &filter_method);
        rows = png_get_rows(png_ptr, info_ptr);
        int rowbytes;
        rowbytes = png_get_rowbytes(png_ptr, info_ptr);
        for (int y = 0; y < height; ++y) {
            png_bytep row;
            row = rows[y];
            for (int x = 0; x < rowbytes; ++x) {
                pixel = getPixel(bitmap, x, y, z);
                pixel->value = row[x];
                pixel->v = malloc(sizeof(Vec));
                pixel->v->x = x;
                pixel->v->y = y;
                pixel->v->z = z;
                pixel->grouping = NULL;
                pixel->inSSL = 0;
            }
        }
    }
    return bitmap;
}