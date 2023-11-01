#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>


#include "curses-gfx.h"

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

	double time = M_PI_2;
	
		
//		printf("====== New line plot\n");
		Coordinates2D a = {0,0};
//		Coordinates2D b = {40 + 40*cos(time),20+20*sin(time)};
		Coordinates2D b = {0+5, 0+4};
		time += 0.005;
	
	while (keepRunning == true) {
		
		
		erase();
//
////		printf("b =  %d, %d ", b.x, b.y);
//		ln2(a, b);
//
//		a.y +=2;
//		b.y += 2;
//		ln(a, b);
//
//
//
//
//
//		a.x = 40;
//		a.y = 30;
//		b.x = 300;
//		b.y = 31;
//		ln2(a, b);
//		a.y += 2;
//		b.y += 2;
//		ln(a, b);
		
//		a.x = 1;
//		a.y = 1;
//		b.x = 38*3;
//		b.y = 20*3;
//		ln2(a, b);
//		a.y += 2;
//		b.y += 2;
		ln2(a, b);
		
		a.x = 0;
		a.y = 10;
		b.x = 0 + 5*2;
		b.y = 10 - 4*2;
		ln2(a, b);
		
		a.x = 10;
		a.y = 10;
		b.x = 10 + 5*2;
		b.y = 10 - 4*2;
		ln2(b, a); // good
		
		a.x = 0;
		a.y = 10;
		b.x = 0 + 5*2;
		b.y = 10 + 4*2;
		ln2(a, b); // good
		
		a.x = 10;
		a.y = 10;
		b.x = 10 + 5*2;
		b.y = 10 + 4*2;
		ln2(b, a);
		
//	while(keepRunning) {
		usleep(1000000.0/30.0);
		
		if (getch() == 0x1B) {	// Escape
			keepRunning = false;
		}

	}
	sleep(5);
	cleanupConsole();

	a.x = 0;
	a.y = 10;
	b.x = 0 + 5*2;
	b.y = 10 - 4*2;
	ln2(a, b);
	
	a.x = 10;
	a.y = 10;
	b.x = 10 + 5*2;
	b.y = 10 - 4*2;
	ln2(b, a); // good
	
	a.x = 0;
	a.y = 10;
	b.x = 0 + 5*2;
	b.y = 10 + 4*2;
	ln2(a, b); // good
	
	a.x = 10;
	a.y = 10;
	b.x = 10 + 5*2;
	b.y = 10 + 4*2;
	ln2(b, a);
	
	return 0;
}
