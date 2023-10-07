#ifndef CURSES_GFX_H
#define CURSES_GFX_H

#include "curses-gfx-types.h"

typedef struct _Coordinates2D {
    int x;
    int y;
} Coordinates2D;


void drawDotFloat(double x, double y);

void ln2( Coordinates2D a, Coordinates2D b);

void set( Coordinates2D pt, char c);
char getp( Coordinates2D* pts, double err);
void ln( Coordinates2D a, Coordinates2D b);

class CursesGfxTerminal {
private:
    bool configured;
    static int numColors;
    short* restoreR;
    short* restoreG;
    short* restoreB;
    
public:
    void setupTerminal();
    void cleanupTerminal();
    
    static void setRGB( const Coordinates2D& pixel, const Coordinates3D& rgb);
    
    CursesGfxTerminal();
    ~CursesGfxTerminal();
    
};

#endif
