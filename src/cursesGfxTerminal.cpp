#include "curses-gfx.h"
#include "curses-gfx-3d.h"

#include <ncurses.h>
#include <locale.h>
#include <cmath>
#include <iostream>


int CursesGfxTerminal::numColors = 1;
int CursesGfxTerminal::numHueLevels = 16;
int CursesGfxTerminal::numSatLevels = 16;
int CursesGfxTerminal::satShift = 4;
int CursesGfxTerminal::satMask = 0xF0;
int CursesGfxTerminal::hueMask = 0x0F;
int CursesGfxTerminal::enabledColor = 0;

CursesGfxTerminal::CursesGfxTerminal()
: configured(false) {
    restoreR = restoreG = restoreB = NULL;
    
}
CursesGfxTerminal::~CursesGfxTerminal() {
    if(configured) {
        cleanupTerminal();
    }
}

void CursesGfxTerminal::setupTerminal() {
    setlocale(LC_ALL, "");
    
    // Start up Curses window
    initscr();
    cbreak();
    noecho();
    nodelay(stdscr, 1);    // Don't wait at the getch() function if the user hasn't hit a key
    keypad(stdscr, 1); // Allow Function key input and arrow key input

    start_color();
    use_default_colors();
//
//    init_pair(1, COLOR_RED, COLOR_BLACK);
//    init_pair(2, COLOR_GREEN, COLOR_BLACK);
//    init_pair(3, COLOR_CYAN, COLOR_BLACK);
//    init_pair(4, COLOR_BLUE, COLOR_WHITE);
//
//    init_pair(5, COLOR_BLACK, COLOR_RED );
//    init_pair(6, COLOR_BLACK, COLOR_GREEN );
//    init_pair(7, COLOR_BLACK, COLOR_CYAN );
//    init_pair(8, COLOR_WHITE, COLOR_BLUE );
    int line = 0;
    mvprintw(line++, 0, "ncurses COLORS: %d\n", COLORS);
    mvprintw(line++, 0, "ncurses COLOR_PAIRS: %d\n", COLOR_PAIRS  );
    mvprintw(line++, 0, "ncurses has_colors: %d\n", has_colors());
    mvprintw(line++, 0, "ncurses can_change_color: %d\n", can_change_color());
    
    numColors = COLORS;
//    short restoreR[COLORS];
//    short restoreG[COLORS];
//    short restoreB[COLORS];
    restoreR = new short[COLORS];
    restoreG = new short[COLORS];
    restoreB = new short[COLORS];
    for(short i = 0; i < COLORS; i++) {
        color_content(i, &restoreR[i], &restoreG[i], &restoreB[i]);
    }
    
//    init_color(0, 1000,1000,1000);
    init_pair(0, 0, -1);    // White
    
//    int numSatLevels = 8;
    int numSatBits = log2(numColors)*3/8;
    int numHueBits = log2(numColors) - numSatBits;
    satShift = numHueBits;
    hueMask = ((1 << numHueBits) - 1) & 0xFF;
    satMask = ((numColors-1) - hueMask) & 0xFF;
    numSatLevels = (1 << numSatBits);
    numHueLevels = numColors/numSatLevels;

    for(int i = 1; i < COLORS; i++) {
////        int hueIndex = mod((i-1)*numSatLevels,COLORS);
////        int satIndex = ceil((double)(i-1)/((double)COLORS/(double)(numSatLevels)));
//        int hueIndex = mod(i*(numSatLevels),COLORS);
//        int satIndex = ceil((double)i/((double)COLORS/(double)(numSatLevels)));
//        double hue = (double)(hueIndex)/(double)(COLORS) * 360.0;
//        double sat = ((double)satIndex)/((double)(numSatLevels));
        
        
        double hue = ((double)(((uint8_t)i) & hueMask ) - 0.5)*360.0/(double)numHueLevels;
        if(hue < 0) {
            hue += 360;
        }
        double sat = ((double)((((uint8_t)i) & satMask ) >> satShift) + 1) /(double)numSatLevels;
        
        Coordinates3D rgb = hsvToRgb({hue, sat, 1});
//        if(i == 0)
//            init_color(i, 1000,1000,1000);
//        else
            init_color(i, (int)(rgb.x/255.0*1000.0), (int)(rgb.y/255.0*1000.0), (int)(rgb.z/255.0*1000.0));
        init_pair(i, i, -1);
    }
    
    numSatLevels += 1; // HACK we techincally have 9 levels since sat of 0 == white
    
//    for(int i = 0; i < COLORS; i++) {
//
//        attron(COLOR_PAIR(i));
////        attron(A_DIM);
//        mvaddch(line, i, 'W');
////        attroff(COLOR_PAIR(i));
//    }
    curs_set(0);    // no cursor

//    getch();
    
    
    
    configured = true;
}
void CursesGfxTerminal::cleanupTerminal() {
    if(restoreR) {
        for(int i = 1; i < numColors; i++) {
            init_color(i, restoreR[i], restoreG[i], restoreB[i]);
        }
        delete [] restoreR;
        delete [] restoreG;
        delete [] restoreB;

    }
    
    use_default_colors();
    clear();
    endwin();

    std::cout << "CursesGfxTerminal::cleanupTerminal(): Console has been cleaned!" << std::endl;
    
    configured = false;
}



struct AsciiLutInfo {
    char c;
    uint8_t type;
};
//#ifdef __APPLE__

struct AsciiLutInfo asciiLUT[256] = {{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'.', 2},{'.', 2},{'.', 2},{'.', 2},{'.', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'\'', 2},{'\'', 2},{'\'', 2},{'\'', 2},{'\'', 2},{':', 2},{':', 2},{':', 2},{',', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'~', 2},{'~', 2},{'~', 2},{'~', 2},{';', 2},{';', 2},{';', 2},{';', 2},{';', 2},{';', 2},{'|', 2},{'|', 2},{'|', 2},{'|', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'/', 2},{'/', 2},{'\\', 2},{'\\', 2},{'^', 2},{'=', 2},{'i', 2},{'i', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'l', 2},{'l', 2},{'l', 2},{'l', 2},{'L', 2},{'L', 2},{'c', 2},{'*', 2},{'*', 2},{'*', 2},{'v', 2},{'v', 2},{'v', 2},{'J', 2},{'J', 2},{'J', 2},{'T', 2},{'T', 2},{'7', 2},{'t', 2},{'t', 2},{'f', 2},{'1', 2},{'j', 2},{'j', 2},{'Y', 2},{'s', 2},{'u', 2},{'n', 2},{'(', 2},{'I', 2},{'C', 2},{'C', 2},{'}', 2},{'}', 2},{'{', 2},{'o', 2},{'o', 2},{'o', 2},{'e', 2},{'e', 2},{']', 2},{'y', 2},{'2', 2},{'2', 2},{'2', 2},{'a', 2},{'a', 2},{'S', 2},{'E', 2},{'h', 2},{'k', 2},{'V', 2},{'Z', 2},{'Z', 2},{'Z', 2},{'5', 2},{'5', 2},{'5', 2},{'5', 2},{'U', 2},{'U', 2},{'U', 2},{'p', 2},{'G', 2},{'q', 2},{'w', 2},{'d', 2},{'d', 2},{'A', 2},{'A', 2},{'H', 2},{'H', 2},{'K', 2},{'K', 2},{'K', 2},{'O', 2},{'O', 2},{'O', 2},{'D', 2},{'#', 2},{'9', 2},{'6', 2},{'$', 2},{'$', 2},{'$', 2},{'R', 2},{'R', 2},{'R', 2},{'R', 2},{'B', 2},{'B', 2},{'B', 2},{'g', 2},{'8', 2},{'8', 2},{'8', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'%', 2},{'%', 2},{'%', 2},{'%', 2},{'@', 2},{'0', 2},{'0', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2}};

//#else
//
//struct AsciiLutInfo asciiLUT[256] = {{' ', 1},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'.', 1},{'.', 1},{'.', 1},{'-', 1},{'-', 1},{'-', 1},{'\'', 1},{'\'', 1},{':', 1},{',', 1},{'_', 1},{'`', 2},{'`', 2},{'.', 2},{'"', 1},{'"', 1},{';', 1},{';', 1},{';', 1},{'|', 1},{'|', 1},{'!', 1},{'!', 1},{'-', 2},{'+', 1},{'\\', 1},{'^', 1},{'=', 1},{'i', 1},{'\'', 2},{'r', 1},{'r', 1},{'?', 1},{'l', 1},{'L', 1},{'_', 2},{'v', 1},{'J', 1},{'J', 1},{'T', 1},{'f', 1},{'j', 1},{'x', 1},{'z', 1},{'C', 1},{'I', 1},{'o', 1},{'~', 2},{'y', 1},{'2', 1},{'[', 1},{'a', 1},{'3', 1},{'V', 1},{'P', 1},{';', 2},{'5', 1},{'U', 1},{'G', 1},{'m', 1},{'d', 1},{'A', 1},{'K', 1},{'!', 2},{'O', 1},{'D', 1},{'$', 1},{'$', 1},{'R', 1},{'g', 1},{'+', 2},{'8', 1},{'8', 1},{'N', 1},{'/', 2},{'Q', 1},{'&', 1},{'@', 1},{'=', 2},{'i', 2},{'i', 2},{'0', 1},{'r', 2},{'M', 1},{'M', 1},{'W', 1},{'W', 1},{'W', 1},{'W', 1},{'l', 2},{'l', 2},{'L', 2},{'L', 2},{'c', 2},{'*', 2},{'*', 2},{'*', 2},{'v', 2},{'v', 2},{'v', 2},{'J', 2},{'J', 2},{'J', 2},{'T', 2},{'T', 2},{'7', 2},{'t', 2},{'t', 2},{'f', 2},{'1', 2},{'j', 2},{'j', 2},{'Y', 2},{'F', 2},{'u', 2},{'n', 2},{'(', 2},{'I', 2},{'C', 2},{'C', 2},{'}', 2},{'}', 2},{'{', 2},{'o', 2},{'o', 2},{'o', 2},{'e', 2},{'e', 2},{']', 2},{'y', 2},{'2', 2},{'2', 2},{'2', 2},{'a', 2},{'a', 2},{'S', 2},{'E', 2},{'h', 2},{'k', 2},{'V', 2},{'Z', 2},{'Z', 2},{'Z', 2},{'5', 2},{'5', 2},{'5', 2},{'5', 2},{'U', 2},{'U', 2},{'U', 2},{'p', 2},{'G', 2},{'q', 2},{'w', 2},{'d', 2},{'d', 2},{'A', 2},{'A', 2},{'H', 2},{'H', 2},{'K', 2},{'K', 2},{'K', 2},{'O', 2},{'O', 2},{'O', 2},{'D', 2},{'#', 2},{'9', 2},{'6', 2},{'$', 2},{'$', 2},{'$', 2},{'R', 2},{'R', 2},{'R', 2},{'R', 2},{'B', 2},{'B', 2},{'B', 2},{'g', 2},{'8', 2},{'8', 2},{'8', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'%', 2},{'%', 2},{'%', 2},{'%', 2},{'@', 2},{'0', 2},{'0', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2}};
//#endif //__APPLE__

int CursesGfxTerminal::rgbToColorIndex(const Coordinates3D& rgb, double& outputLevel) {
    Coordinates3D clippedRGB = clipRGB(rgb);
    Coordinates3D hsl = rgbToHsv(clippedRGB);

    int colorIndex = 0;
    if(hsl.y != 0 ) {
        colorIndex = floor(hsl.x/360.0*(double)(numHueLevels));
        colorIndex |= (((uint8_t)floor(hsl.y*((double)numSatLevels-2))) << satShift) & satMask;
    }
    outputLevel = hsl.z;
    return colorIndex;
}

void CursesGfxTerminal::setRGB( const Coordinates2D& pixel, const Coordinates3D& rgb) {

    double level;
//    int colorIndex = rgbToColorIndex(rgb, level);
//    attron(COLOR_PAIR(colorIndex));
    enableColor(rgb, level);
    AsciiLutInfo aLut = asciiLUT[(int)(level*255)];
    set(pixel, aLut.c);
    
//    attroff(COLOR_PAIR(colorIndex));
    disableColor();
}

void CursesGfxTerminal::enableColor(const Coordinates3D& rgb, double& outputLevel) {
    enabledColor = rgbToColorIndex(rgb, outputLevel);
//    if(enabledColor) {
    attron(A_BOLD);
    attron(COLOR_PAIR(enabledColor));
//    }
}

void CursesGfxTerminal::disableColor() {
//    if(enabledColor) {
    attroff(COLOR_PAIR(enabledColor));
    attroff(A_BOLD);
//    }
}
