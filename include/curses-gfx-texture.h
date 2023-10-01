#ifndef CURSES_GFX_TEXTURE_H
#define CURSES_GFX_TEXTURE_H

#include "curses-gfx-3d.h"


class Texture {
public:
    ColorRGBA* data;
    
    int width, height;
    
    Texture();
    Texture(int width, int height);
    Texture(const char* pngFilename);
    ~Texture();
    
    void resize(int dimX, int dimY);
    void set(const double& x, const double& y, const ColorRGBA& value);
    
    void setPixel(const int& x, const int& y, const ColorRGBA& value);
    
    void monochromize();
    
    int min(int a, int b);
    int max(int a, int b);
    
    void normalize(double factor);
    void offsetAvergageToCenter(double factor);
    
    void offset(int o);
    void scale(int s);
    void invert();
    
    ColorRGBA sample( const double& x, const double& y);
    
    bool loadPng(const char* filename); // 0 success, 1 failure;
};


#endif
