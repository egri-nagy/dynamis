#ifndef PNG_UTIL_H
#define PNG_UTIL_H

#include <png.h>

int write_indexed_png(char *file_name, png_byte *image, png_uint_32 width, png_uint_32 height, png_byte *palette, png_byte palette_length);

#endif //PNG_UTIL_H
