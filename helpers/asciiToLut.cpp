/*
 * A simple libpng example program
 * http://zarb.org/~gc/html/libpng.html
 *
 * Modified by Yoshimasa Niwa to make it much simpler
 * and support all defined color_type.
 *
 * To build, use the next instruction on OS X.
 * $ brew install libpng
 * $ clang -lz -lpng16 libpng_test.c
 *
 * Copyright 2002-2010 Guillaume Cottenceau.
 *
 * This software may be freely redistributed under the terms
 * of the X11 license.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <png.h>
#include <algorithm>
#include <vector>

int width, height;
png_byte color_type;
png_byte bit_depth;
png_bytep *row_pointers = NULL;

void read_png_file(char *filename) {
  FILE *fp = fopen(filename, "rb");

  png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if(!png) abort();

  png_infop info = png_create_info_struct(png);
  if(!info) abort();

  if(setjmp(png_jmpbuf(png))) abort();

  png_init_io(png, fp);

  png_read_info(png, info);

  width      = png_get_image_width(png, info);
  height     = png_get_image_height(png, info);
  color_type = png_get_color_type(png, info);
  bit_depth  = png_get_bit_depth(png, info);

  // Read any color_type into 8bit depth, RGBA format.
  // See http://www.libpng.org/pub/png/libpng-manual.txt

  if(bit_depth == 16)
    png_set_strip_16(png);

  if(color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(png);

  // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
  if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png);

  if(png_get_valid(png, info, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png);

  // These color_type don't have an alpha channel then fill it with 0xff.
  if(color_type == PNG_COLOR_TYPE_RGB ||
     color_type == PNG_COLOR_TYPE_GRAY ||
     color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

  if(color_type == PNG_COLOR_TYPE_GRAY ||
     color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png);

  png_read_update_info(png, info);

    if (row_pointers) free(row_pointers);

  row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
  for(int y = 0; y < height; y++) {
    row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
  }

  png_read_image(png, row_pointers);

  fclose(fp);

  png_destroy_read_struct(&png, &info, NULL);
}

void write_png_file(char *filename) {
  int y;

  FILE *fp = fopen(filename, "wb");
  if(!fp) abort();

  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png) abort();

  png_infop info = png_create_info_struct(png);
  if (!info) abort();

  if (setjmp(png_jmpbuf(png))) abort();

  png_init_io(png, fp);

  // Output is 8bit depth, RGBA format.
  png_set_IHDR(
    png,
    info,
    width, height,
    8,
    PNG_COLOR_TYPE_RGBA,
    PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_DEFAULT,
    PNG_FILTER_TYPE_DEFAULT
  );
  png_write_info(png, info);

  // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
  // Use png_set_filler().
  //png_set_filler(png, 0, PNG_FILLER_AFTER);

  if (!row_pointers) abort();

  png_write_image(png, row_pointers);
  png_write_end(png, NULL);

  for(int y = 0; y < height; y++) {
    free(row_pointers[y]);
  }
  free(row_pointers);

  fclose(fp);

  png_destroy_write_struct(&png, &info);
}
struct CharIntense {
    unsigned char c;
    double intensity;
    int type;
};
bool compareByI(const CharIntense &a, const CharIntense &b)
{
    return a.intensity < b.intensity;
}
void process_png_file(std::vector<CharIntense>& pairs, int type) {
    int characterWidth = 12;
    int characterCount = floor(width/characterWidth);
    int characterHeight = 28;
    printf("There are %d characters\n", characterCount-1);
    printf("Height of file %d \n", height);
    
    
    
    for(int c = 0; c < characterCount; c++) {
        double totalBrightness = 0;
        for(int x = characterWidth*c; x < (characterWidth*c + characterWidth); x++) {
            for(int y = 0; y < characterHeight; y++) {
                png_bytep row = row_pointers[y];
                png_bytep px = &(row[x * 4]);
                
                for (int j = 0; j < 3; j++) {
                    totalBrightness += 255*pow((double)px[j]/255.0,1.);
                }
                // Do something awesome for each pixel here...
//                printf("%4d, %4d = RGBA(%3d, %3d, %3d, %3d)\n", x, y, px[0], px[1], px[2], px[3]);
            }
        }
        totalBrightness /= characterWidth*characterHeight*3;
//        totalBrightness *= 255/96.660714;
        if(c == 0) totalBrightness = 0; // HACK first character was the cursor(space)
//        printf("c: % 3d %c %02f\n", c, c+' ', totalBrightness);
        
        CharIntense ci;
        ci.c = c+' ';
        ci.intensity = totalBrightness;
        ci.type = type;
        pairs.push_back( ci );
    }
    
    
    
}

int main(int argc, char *argv[]) {
//    if(argc != 1) abort();
    std::vector<CharIntense> pairs;
    
//    read_png_file("macOsGreyTerminal.png");
//    process_png_file(pairs, 1);
    
    read_png_file("macOsTerminal.png");
    process_png_file(pairs, 2);
    
    double max = 0;
    
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
        if(max < it->intensity)
            max = it->intensity;
    }
    std::sort(pairs.begin(), pairs.end(), compareByI);
    int i = 0;
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
        CharIntense* ci = &(*it);
        ci->intensity *= 255/max;
        if(ci->type == 1)
            printf("\x1b[2;37;40m");
        printf("%d % 2d: %c %02f\n", i++, ci->type, ci->c, ci->intensity);
        if(ci->type ==1 )
            printf("\x1b[0m");
    }
    
    
    printf("\e[1;1H\e[2J");
    for(int i = 0; i < 5; i++) {
        int luminance = 0;
        for (auto it = pairs.begin(); (it != pairs.end()); ) {
            CharIntense ci = *it;
            CharIntense ciNext = *(it+1);
            while(luminance > (ci.intensity + ciNext.intensity)/2 && it != pairs.end()-1 ) {
                it++;
                ci = *it;
                ciNext = *(it+1);
            }
            
            if(ci.type == 1)
                printf("\x1b[2;37;40m");
            printf("%c", ci.c);
            
            
            if(ci.type == 1 )
                printf("\x1b[0m");
            luminance++;
            if (luminance < 256) {
//                printf(",");
            } else {
                break;
            }
            
        }
        printf("\n");
    }
    
    printf("struct <type>[256] = {");
    //    for(int i = 0; i < 5; i++) {
    int luminance = 0;
        for (auto it = pairs.begin(); (it != pairs.end()); ) {
            CharIntense ci = *it;
            CharIntense ciNext = *(it+1);
            while(luminance > (ci.intensity + ciNext.intensity)/2 && it != pairs.end()-1 ) {
                it++;
                ci = *it;
                ciNext = *(it+1);
            }
            
                if(ci.type == 1)
                    printf("\x1b[2;37;40m");
//                printf("%c", ci.c);
            
            if(ci.c == '\\') {
                printf("{'\\\\', %d}", ci.type);
            } else if(ci.c == '\'') {
                printf("{'\\'', %d}", ci.type);
            } else  {
                printf("{'%c', %d}", ci.c, ci.type);
            }
            
            
            if(ci.type == 1 )
                printf("\x1b[0m");
            luminance++;
            if (luminance < 256) {
                printf(",");
            } else {
                break;
            }
            
        }
        printf("};\n");
//    }
    
    
    return 0;
}
