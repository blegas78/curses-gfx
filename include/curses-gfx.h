#ifndef CURSES_GFX_H
#define CURSES_GFX_H

#include "curses-gfx-types.h"




void drawDotFloat(double x, double y);

void ln2( Coordinates2D a, Coordinates2D b);

void set(const Coordinates2D& pt, const char& c);
char getp( Coordinates2D* pts, double err);
void ln( Coordinates2D a, Coordinates2D b);

class CursesGfxTerminal {
private:
    bool configured;
    static int numColors;
    static int numSatLevels;
    static int numHueLevels;
    static int satShift;
    static int satMask;
    static int hueMask;
    static int enabledColor;
    short* restoreR;
    short* restoreG;
    short* restoreB;
    
public:
    void setupTerminal();
    void cleanupTerminal();
    
    static int rgbToColorIndex(const Coordinates3D& rgb, double& outputLevel);
    static void setRGB( const Coordinates2D& pixel, const Coordinates3D& rgb);
    
    static void enableColor(const Coordinates3D& rgb, double& outputLevel);
    static void disableColor();
    
    CursesGfxTerminal();
    ~CursesGfxTerminal();
    
    static const int& hueLevels() { return numHueLevels; }
    static const int& satLevels() { return numSatLevels; }
    
};

#endif
