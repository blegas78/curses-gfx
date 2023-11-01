#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>


#include "curses-gfx.h"
#include "curses-clock.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}


void setupTerminal()
{
	// Start up Curses window
	initscr();
	cbreak();
	noecho();
	nodelay(stdscr, 1);	// Don't wait at the getch() function if the user hasn't hit a key
	keypad(stdscr, 1); // Allow Function key input and arrow key input

	start_color();
	init_pair(1, COLOR_RED, COLOR_BLACK);
	init_pair(2, COLOR_GREEN, COLOR_BLACK);
	init_pair(3, COLOR_CYAN, COLOR_BLACK);
	init_pair(4, COLOR_BLUE, COLOR_WHITE);

	init_pair(5, COLOR_BLACK, COLOR_RED );
	init_pair(6, COLOR_BLACK, COLOR_GREEN );
	init_pair(7, COLOR_BLACK, COLOR_CYAN );
	init_pair(8, COLOR_WHITE, COLOR_BLUE );

	curs_set(0);	// no cursor

//	atexit(destroy);
}

void cleanupConsole() {
	clear();
	endwin();

	std::cout << "Console has been cleaned!" << std::endl;
}


int main(void) {
	
	setupTerminal();
	
	struct timeval time;
	struct timeval startTime;
	struct timeval timeStopwatch = {0,0};
	
	long double runningTime = 0;
	gettimeofday(&startTime, NULL);
//	startTime.tv_sec += 60*60*5;
	
	while (keepRunning == true) {
		
		usleep(1000000.0/30.0);
		erase();
		
//		oldtime = time;
		gettimeofday(&time, NULL);
		runningTime = (double)(time.tv_sec-startTime.tv_sec) +  ((double)(time.tv_usec-startTime.tv_usec))/1000000.0;
		runningTime *= 60;
		double remainder;
		timeStopwatch.tv_usec = modf(runningTime, &remainder) * 1000000;
		timeStopwatch.tv_sec = remainder-60*60*5;
		
		Coordinates2D center = {51, 22};
		drawclock(50*0.75, 21.5*0.75, center, true, time);
		
		for (double i = 0; i < 4; i++) {
			Coordinates2D center2;
			center2.x = center.x + cos((45+i*90) * M_PI/180.0) * 50;
			center2.y = center.y + sin((45+i*90) * M_PI/180.0) * 21.5;
			drawclock(50*0.25, 21.5*0.25, center2, true, timeStopwatch);
		}
		
		if (getch() == 0x1B) {	// Escape
			keepRunning = false;
		}
	}
	
	cleanupConsole();
	return 0;
}
