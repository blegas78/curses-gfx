#include "curses-gfx-texture.h"
#include "curses-gfx-png-loader.h"
//#include <png.h>


Texture::Texture() {
    data = new ColorRGBA[1];
    width = 1;
    height = 1;
}
Texture::Texture(int width, int height)
: width(width), height(height) {
    data = new ColorRGBA[width*height];
    
}

Texture::Texture(const char* pngFilename)
:data(NULL), width(0), height(0) {
    loadPng(pngFilename);
}

Texture::~Texture() {
    width = 0;
    height = 0;
    delete [] data;
}

void Texture::set(const double& x, const double& y, const ColorRGBA& value) {
    int xPart = mod(x*(double)width,width);
    int yPart = mod(y*(double)height,height);
    printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
    data[ xPart + width*yPart] = value;
}

void Texture::setPixel(const int& x, const int& y, const ColorRGBA& value) {
    int xPart = mod(x,width);
    int yPart = mod(y,height);
//        printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
    data[ xPart + width*yPart] = value;
}

void Texture::monochromize() {
    for(int i = 0;  i < width*height; i++) {
        double avg = ((double)data[i].r + (double)data[i].g + (double)data[i].b)/3.0; // MONOCHROME
        data[i].r = avg;
        data[i].g = avg;
        data[i].b = avg;
    }
}

int Texture::min(int a, int b) {
    return a < b ? a : b;
}
int Texture::max(int a, int b) {
    return a > b ? a : b;
}

void Texture::normalize(double factor) {
    int maxi = 0, mini = 255;
    for(int i = 0;  i < width*height; i++) {
        int L = max(max(data[i].r,data[i].g),data[i].b);
        maxi = max(L,maxi);
        mini = min(L,mini);// tricky one here (min of maxes)
    }
    double scale = (255.0/(double)(maxi - mini)) * factor + (1.0-factor)*1.0;
    mini = mini * (factor);
    for(int i = 0;  i < width*height; i++) {
        data[i].r = min(max(round(scale * (double)(data[i].r - mini)), 0), 255);
        data[i].g = min(max(round(scale * (double)(data[i].g - mini)), 0), 255);
        data[i].b = min(max(round(scale * (double)(data[i].b - mini)), 0), 255);
        
    }
}

void Texture::offsetAvergageToCenter(double factor) {
    int avg = 0;
    for(int i = 0;  i < width*height; i++) {
        avg += max(max(data[i].r,data[i].g),data[i].b);
    }
    avg = round((double)avg/(double)(width*height) * factor);
//    for(int i = 0;  i < width*height; i++) {
//        data[i].r = min(max((int)data[i].r - avg, 0), 255);
//        data[i].g = min(max((int)data[i].g - avg, 0), 255);
//        data[i].b = min(max((int)data[i].b - avg, 0), 255);
//    }
    offset(-avg);
}

void Texture::offset(int o) {
    for(int i = 0;  i < width*height; i++) {
        data[i].r = min(max((int)data[i].r + o, 0), 255);
        data[i].g = min(max((int)data[i].g + o, 0), 255);
        data[i].b = min(max((int)data[i].b + o, 0), 255);
    }
}
void Texture::scale(int s) {
    for(int i = 0;  i < width*height; i++) {
        data[i].r = min(max(round(s * (double)data[i].r), 0), 255);
        data[i].g = min(max(round(s * (double)data[i].g), 0), 255);
        data[i].b = min(max(round(s * (double)data[i].b), 0), 255);
        
    }
}

void Texture::invert() {
    for(int i = 0;  i < width*height; i++) {
        data[i].r = 255 - data[i].r;
        data[i].g = 255 - data[i].g;
        data[i].b = 255 - data[i].b;
    }
}

ColorRGBA Texture::sample( const double& x, const double& y) {
    int xPart = mod(x*(double)width,width);
    int yPart = mod(y*(double)height,height);
    return data[ xPart + width*yPart];
}

void Texture::resize(int dimX, int dimY) {
    if(data) {
        width = height = 0; // prevent race conditions
        delete [] data;
    }
    data = new ColorRGBA[dimX * dimY];
    width = dimX;
    height = dimY;
}

bool Texture::loadPng(const char* filename) {
    PngLoader mPngLoader;
    if(read_png_file(filename, mPngLoader)) {
        if (data == NULL) {
            data = new ColorRGBA[1];
            width = 1;
            height = 1;
            data[0] = {127,127,127,0};  // for testing
        }
        return 1;
    }
    
    if(mPngLoader.height*mPngLoader.width != width*height) {
        if (data) {
            width = height = 0; // prevent race conditions
            delete [] data;
        }
        data = new ColorRGBA[mPngLoader.height*mPngLoader.width];
        width = mPngLoader.width;
        height = mPngLoader.height;
    }
    for(int y = 0; y < mPngLoader.height; y++) {
        png_bytep row = mPngLoader.row_pointers[mPngLoader.height - y - 1]; // HACK reverse y
        for(int x = 0; x < mPngLoader.width; x++) {
            png_bytep px = &(row[x * 4]);
            setPixel(y, x, {px[0], px[1], px[2], px[3]});   // HACK swap x/y
        }
    }
    
    delete [] mPngLoader.row_pointers;
    
    return 0;
}
