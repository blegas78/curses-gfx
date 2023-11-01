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

enum StopwatchState {
	SWS_IDLE,
	SWS_RUNNING,
	SWS_PAUSED
};

int main(void) {
	
	setupTerminal();
	
	
	struct timeval currentTime;
	struct timeval time;
	struct timeval startTime;
	struct timeval resumeTime;
	struct timeval timeStopwatch = {0,0};
	struct timeval lapTime = {7*60*60, 0};
	
	long double runningTime = 0;
	gettimeofday(&startTime, NULL);
//	startTime.tv_sec += 60*60*5;
	
	StopwatchState state = SWS_IDLE;
	char ch;
	
	
	double aspect = 12.0/28.0;	// MacOS
//	double aspect = 14.0/28.0;	// raspbian
//	double aspect = 4.0/6.0;	// ZipitZ2

	int screenSizeX, screenSizeY;
	getmaxyx(stdscr, screenSizeY, screenSizeX);
	double size = (screenSizeY-1)/(2*aspect);	// radius in x
	Coordinates2D center = {screenSizeX/2, (int)(size*aspect)};
	
	
	while (keepRunning == true) {
		
		usleep(1000000.0/1000.0);
		erase();
		
		// Current time
		Coordinates2D centerMiniClock;
		centerMiniClock.x = center.x - size/2.4;
		centerMiniClock.y = center.y - size*aspect/2.4+1;
		
		gettimeofday(&currentTime, NULL);
		drawclock(size/3.5, size*aspect/3.5, centerMiniClock, true, currentTime);
		
		// Lap time
		Coordinates2D centerLapTime;
		centerLapTime.x = center.x - size/2.4;
		centerLapTime.y = center.y + size*aspect/2.4;
		drawclock(size/3.5, size*aspect/3.5, centerLapTime, true, lapTime);
		
//		oldtime = time;
		if (state == SWS_RUNNING) {
//			gettimeofday(&time, NULL);
			time = currentTime;
			runningTime = (double)(time.tv_sec-startTime.tv_sec) +  ((double)(time.tv_usec-startTime.tv_usec))/1000000.0;
			double remainder;
			
			timeStopwatch.tv_usec = modf(runningTime, &remainder) * 1000000;
			timeStopwatch.tv_sec = remainder;
		}
		
		struct timeval timeStopwatchCorrected = timeStopwatch;
		timeStopwatchCorrected.tv_sec += 60*60*7;	// remove 5  hours
		
		Coordinates2D centerSeconds;
		centerSeconds.x = center.x + size/2.4 +1;
		centerSeconds.y = center.y - size/2.4*aspect+1;
		drawSecondsDial(size/3.5, size*aspect/3.5, centerSeconds, true, timeStopwatchCorrected);
		
		drawclock(size, size*aspect, center, true, timeStopwatchCorrected);
		
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == ' ') {
			switch (state) {
				case SWS_IDLE:
					gettimeofday(&startTime, NULL);
					state = SWS_RUNNING;
					break;
					
				case SWS_RUNNING:
					state = SWS_PAUSED;
					break;
					
				case SWS_PAUSED:
					gettimeofday(&startTime, NULL);
					startTime.tv_sec -= timeStopwatch.tv_sec;
					if (startTime.tv_usec < timeStopwatch.tv_sec) {
						startTime.tv_sec--;
						startTime.tv_usec = 1000000 - timeStopwatch.tv_usec;
					} else {
						startTime.tv_usec -= timeStopwatch.tv_usec;
					}
					state = SWS_RUNNING;
					break;
					
			}
		} else if (ch == 'r' || ch == 'R') {
			switch (state) {
				case SWS_RUNNING:
					lapTime = timeStopwatchCorrected;
					break;
					
				default:
					timeStopwatch.tv_usec = 0;
					timeStopwatch.tv_sec = 0;
					lapTime.tv_sec =  60*60*7;
					lapTime.tv_usec = 0;
					state = SWS_IDLE;
					break;
			}
		}
	}
	
	cleanupConsole();
	return 0;
}
