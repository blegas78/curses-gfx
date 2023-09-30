#ifndef CURSES_GFX_TEXTURE_H
#define CURSES_GFX_TEXTURE_H

#include "curses-gfx-3d.h"
#include "curses-gfx-png-loader.h"

class Texture {
public:
    ColorRGBA* data;
    
    int width, height;
    
    Texture(int width, int height);
    ~Texture();
    
    void set(const double& x, const double& y, const ColorRGBA& value);
    
    void setPixel(const int& x, const int& y, const ColorRGBA& value);
    
    void monochromize();
    
    int min(int a, int b);
    int max(int a, int b);
    
    void normalize(double factor);
    
    void offsetAvergageToCenter(double factor);
    
    ColorRGBA sample( const double& x, const double& y);
};


#endif
