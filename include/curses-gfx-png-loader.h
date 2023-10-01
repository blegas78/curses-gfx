#ifndef CURSES_GFX_IMAGE_LOADER_H
#define CURSES_GFX_IMAGE_LOADER_H

#include <png.h>

class PngLoader {
public:
    PngLoader();
    ~PngLoader();
    int width, height;
    png_byte color_type;
    png_byte bit_depth;
    png_bytep *row_pointers;
} ;

int read_png_file(const char *filename, PngLoader& mPngLoader);

#endif
