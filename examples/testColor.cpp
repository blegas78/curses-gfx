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
		
        int line = 3;
        int x = -1;
        for(int i = 0; i < COLORS + CursesGfxTerminal::hueLevels(); i++) {
            if((i % CursesGfxTerminal::hueLevels()) == 0) {
            x++;
        }
            
            if(i < CursesGfxTerminal::hueLevels())
                attron(COLOR_PAIR(0));
            else
                attron(COLOR_PAIR(i-CursesGfxTerminal::hueLevels()));
            
            mvaddch(line, i + x, 'W');
            
            if(i < CursesGfxTerminal::hueLevels())
                attroff(COLOR_PAIR(0));
            else
                attroff(COLOR_PAIR(i-CursesGfxTerminal::hueLevels()));
            
        }
        
        line++;
        
        
        x = 0;
//        for(int j = 0; j < 32; j++) {
//            CursesGfxTerminal::setRGB({y++, 11}, {255, 255, 255});
//        }
        
        for(int i = 0; i < CursesGfxTerminal::satLevels(); i++) {
            double saturation = (double)(i)/(CursesGfxTerminal::satLevels()-1);
            attron(COLOR_PAIR(0));
            mvprintw(line, x+i, "s=%0.2f", saturation);
            attroff(COLOR_PAIR(0));
            for(int j = 0; j < CursesGfxTerminal::hueLevels(); j++)
            {
                double hue = ((double)(j)/(double)CursesGfxTerminal::hueLevels())*360.0;
                Coordinates3D rgb;
                rgb = hsvToRgb({ hue, saturation, 1});
                rgb.x /= 255;
                rgb.y /= 255;
                rgb.z /= 255;
                CursesGfxTerminal::setRGB({x + i, line+1}, rgb);
                
                CursesGfxTerminal::enableThinAscii();
                CursesGfxTerminal::setRGB({x + i, line+2}, rgb);
                CursesGfxTerminal::disableThinAscii();
                
                attron(A_DIM);
                CursesGfxTerminal::setRGB({x + i, line+3}, rgb);
                attroff(A_DIM);
                
                x++;
            }
        }
        
        line += 5;

        // Show cube-ified slices of the RGB cube
        mvprintw(line++, 0, "RGB cube slices:");
        int cubeSlices = 15;
        for(int r = 0; r <= cubeSlices; r++) {
            for(int g = 0; g <= cubeSlices; g++) {
                for(int b = 0; b <= cubeSlices; b++) {
                    Coordinates3D rgb = {(double)r/(double)cubeSlices, (double)g/(double)cubeSlices, (double)b/(double)cubeSlices};
                    
                    
                    mvprintw(line, g*(cubeSlices+2) + 7, "g = %d",  (int)((double)g/(double)cubeSlices*255));
                    mvprintw(line+r+1, 0, "r=%d", (int)((double)r/(double)cubeSlices*255));
                    CursesGfxTerminal::setRGB({7+ b + g*(cubeSlices+2), line + 1 + r}, rgb);
                }
            }
        }
        line += cubeSlices + 3;
        
        
        
        
        // Show cube-ified slices of the HSV cube, but rotated
        mvprintw(line++, 0, "HSV cube slices:");
        for(int h = 0; h <= cubeSlices; h++) {
            for(int s = 0; s <= cubeSlices; s++) {
                for(int v = 0; v <= cubeSlices; v++) {

                    Coordinates3D rgb;
                    rgb = hsvToRgb({ (double)h/(double)(cubeSlices+1)*360,  (double)s/(double)cubeSlices, (double)v/(double)cubeSlices});
                    rgb.x /= 255;
                    rgb.y /= 255;
                    rgb.z /= 255;
//                    CursesGfxTerminal::setRGB({h*(1+satLevels) , 13}, rgb);
                    mvprintw(line, h*(cubeSlices+2) + 7, "h = %0.2f", (double)h/(double)(cubeSlices+1)*360);
                    mvprintw(line+v+1, 0, "v=%0.2f", (double)v/(double)cubeSlices);
                    CursesGfxTerminal::setRGB({7 + s + h*(cubeSlices+2), line + 1 + v}, rgb);
                }
            }
        }
        
        line += cubeSlices + 3;
        
        
        // Again show cube-ified slices of the HSV cube, but rotated
        mvprintw(line++, 0, "HSV cube slices:");
        for(int h = 0; h <= cubeSlices; h++) {
            for(int s = 0; s <= cubeSlices; s++) {
                for(int v = 0; v <= cubeSlices; v++) {

                    Coordinates3D rgb;
                    rgb = hsvToRgb({ (double)h/(double)(cubeSlices+1)*360,  (double)s/(double)cubeSlices, (double)v/(double)cubeSlices});
                    rgb.x /= 255;
                    rgb.y /= 255;
                    rgb.z /= 255;
//                    CursesGfxTerminal::setRGB({h*(1+satLevels) , 13}, rgb);
                    mvprintw(line, s*(cubeSlices+2) + 7, "s = %0.2f", (double)s/(double)cubeSlices);
                    mvprintw(line+h+1, 0, "h=%d", (int)((double)h/(double)(cubeSlices+1)*360));
                    CursesGfxTerminal::setRGB({7 + v + s*(cubeSlices+2), line + 1 + h}, rgb);
                }
            }
        }
        
        line += cubeSlices + 3;
        
        
        // Again show cube-ified slices of the HSV cube
        mvprintw(line++, 0, "HSV cube slices:");
        for(int h = 0; h <= cubeSlices; h++) {
            for(int s = 0; s <= cubeSlices; s++) {
                for(int v = 0; v <= cubeSlices; v++) {

                    Coordinates3D rgb;
                    rgb = hsvToRgb({ (double)h/(double)(cubeSlices+1)*360,  (double)s/(double)cubeSlices, (double)v/(double)cubeSlices});
                    rgb.x /= 255;
                    rgb.y /= 255;
                    rgb.z /= 255;
//                    CursesGfxTerminal::setRGB({h*(1+satLevels) , 13}, rgb);
                    mvprintw(line, v*(cubeSlices+2) + 7, "v = %0.2f", (double)v/(double)cubeSlices);
                    mvprintw(line+s+1, 0, "s=%0.2f", (double)s/(double)cubeSlices);
                    CursesGfxTerminal::setRGB({7 + h + v*(cubeSlices+2), line + 1 + s}, rgb);
                }
            }
        }
        
        line += cubeSlices + 2;
        
        
        
        
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


