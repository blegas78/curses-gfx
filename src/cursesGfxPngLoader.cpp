#include "curses-gfx-png-loader.h"


#include <cstdio>
#include <cstdlib>

// thanks to https://gist.github.com/niw/5963798
int read_png_file(const char *filename, PngLoader& mPngLoader) {
  FILE *fp = fopen(filename, "rb");
    if(fp == 0) {
        return 1;
    }
  png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png) {
        fclose(fp);
        return 1;
    }

  png_infop info = png_create_info_struct(png);
    if(!info) {
        png_destroy_read_struct(&png, NULL, NULL);  // I'm guess on this call, may segfault
        fclose(fp);
        return 1;
    }

    if(setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return 1;
    }

  png_init_io(png, fp);

  png_read_info(png, info);

    mPngLoader.width      = png_get_image_width(png, info);
    mPngLoader.height     = png_get_image_height(png, info);
    mPngLoader.color_type = png_get_color_type(png, info);
    mPngLoader.bit_depth  = png_get_bit_depth(png, info);

  // Read any color_type into 8bit depth, RGBA format.
  // See http://www.libpng.org/pub/png/libpng-manual.txt

  if(mPngLoader.bit_depth == 16)
    png_set_strip_16(png);

  if(mPngLoader.color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(png);

  // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
  if(mPngLoader.color_type == PNG_COLOR_TYPE_GRAY && mPngLoader.bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png);

  if(png_get_valid(png, info, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png);

  // These color_type don't have an alpha channel then fill it with 0xff.
  if(mPngLoader.color_type == PNG_COLOR_TYPE_RGB ||
     mPngLoader.color_type == PNG_COLOR_TYPE_GRAY ||
     mPngLoader.color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

  if(mPngLoader.color_type == PNG_COLOR_TYPE_GRAY ||
     mPngLoader.color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png);

  png_read_update_info(png, info);

  if (mPngLoader.row_pointers)
      delete [] mPngLoader.row_pointers;

    mPngLoader.row_pointers = new png_bytep[mPngLoader.height];// (png_bytep*)malloc(sizeof(png_bytep) * );
  for(int y = 0; y < mPngLoader.height; y++) {
      mPngLoader.row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
  }

  png_read_image(png, mPngLoader.row_pointers);

  fclose(fp);

  png_destroy_read_struct(&png, &info, NULL);
    return 0;
}

