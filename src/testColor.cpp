#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <chrono>

#include <time.h>

#include "curses-gfx.h"
#include "curses-clock.h"
#include "curses-gfx-3d.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}


int main(void) {

    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();


	
	
	auto now = std::chrono::high_resolution_clock::now();
	auto before = now;
	int line;
	
	double angle = 0.0000;
	double modelAngle = 0;
	double distance =5;
	bool autoRotate = false;
	while (keepRunning == true) {
		line = 0;
		usleep(1000000.0/2.0);
		erase();
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
		mvprintw(line++, 0, "FPS: %f", 1000.0/float_ms.count());
		
        
        int x = -1;
        for(int i = 0; i < COLORS + CursesGfxTerminal::hueLevels(); i++) {
            if((i % CursesGfxTerminal::hueLevels()) == 0) {
            x++;
        }
            
            if(i < CursesGfxTerminal::hueLevels())
                attron(COLOR_PAIR(0));
            else
                attron(COLOR_PAIR(i-CursesGfxTerminal::hueLevels()));
            //        attron(A_DIM);
            mvaddch(10, i + x, 'W');
            //        attroff(COLOR_PAIR(i));
            
            if(i < CursesGfxTerminal::hueLevels())
                attroff(COLOR_PAIR(0));
            else
                attroff(COLOR_PAIR(i-CursesGfxTerminal::hueLevels()));
            
        }
        
        
        x = 0;
//        for(int j = 0; j < 32; j++) {
//            CursesGfxTerminal::setRGB({y++, 11}, {255, 255, 255});
//        }
        
        for(int i = 0; i < CursesGfxTerminal::satLevels(); i++) {
            double saturation = (double)(i)/(CursesGfxTerminal::satLevels()-1);
            attron(COLOR_PAIR(0));
            mvprintw(9, x+i, "s=%0.2f", saturation);
            attroff(COLOR_PAIR(0));
            for(int j = 0; j < CursesGfxTerminal::hueLevels(); j++)
            {
                double hue = ((double)(j)/(double)CursesGfxTerminal::hueLevels())*360.0;
                Coordinates3D rgb;
                rgb = hsvToRgb({ hue, saturation, 1});
                rgb.x /= 255;
                rgb.y /= 255;
                rgb.z /= 255;
                CursesGfxTerminal::setRGB({x++ + i, 11}, rgb);
            }
        }
        
        
        // Show cube-ified slices of the HSV cube
//        double hueLevels = 8;
//        double satLevels = 8;
//        double valLevels = 8;
//        for(int h = 0; h < hueLevels; h++) {
//            for(int s = 0; s < satLevels; s++) {
//                for(int v = 0; v < valLevels; v++) {
//
//                    Coordinates3D rgb;
//                    rgb = hsvToRgb({ (double)h/hueLevels*360,  (double)s/satLevels, (double)v/valLevels});
//                    rgb.x /= 255;
//                    rgb.y /= 255;
//                    rgb.z /= 255;
//                    CursesGfxTerminal::setRGB({h*(1+satLevels) , 13}, rgb);
//                }
//            }
//        }
        int cubeSlices = 15;
        for(int r = 0; r <= cubeSlices; r++) {
            for(int g = 0; g <= cubeSlices; g++) {
                for(int b = 0; b <= cubeSlices; b++) {
                    Coordinates3D rgb = {(double)r/(double)cubeSlices, (double)g/(double)cubeSlices, (double)b/(double)cubeSlices};
                    CursesGfxTerminal::setRGB({b + g*(cubeSlices+2), 13 + r}, rgb);
                }
            }
        }
        
        
        
		int ch;
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		}

	}
	
    mCursesGfxTerminal.cleanupTerminal();
	
	printf("angle = %f\n", angle);
	printf("distance = %f\n", distance);
	printf("modelAngle = %f\n", modelAngle);
	
	return 0;
};


