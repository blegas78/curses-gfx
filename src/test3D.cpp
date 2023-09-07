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


void setupTerminal()
{
	setlocale(LC_ALL, "");
	
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

//	Coordinates4D vertices[] = {
//		{8, 2, 0, 1},
//		{-4, 0.0, 0, 1},
//		{0,  -2, 0, 1}
//	};
	
	Coordinates4D vertices[] = {
		{-1, 0, 0, 1},
		{0, 0.0, 0.0, 1},
		{1,  0, 0, 1}
	};
	
	int edgeIndices[][3] = {
		{0, 1, 2}
	};
	
	Polygon4D polygon[1];
	memset(polygon, 0, sizeof(polygon));
	polygon->vertices[0].x = 2;
	polygon->vertices[0].y = 0;
	
	polygon->vertices[1].x = 1;
	polygon->vertices[1].y = 1.73205;
	polygon->vertices[2].x = -1;
	polygon->vertices[2].y = 1.73205;
	
	polygon->vertices[3].x = -2;
	polygon->vertices[3].y = 0;
	polygon->vertices[4].x = -1;
	polygon->vertices[4].y = -1.73205;
	polygon->vertices[5].x = 1;
	polygon->vertices[5].y = -1.73205;
	
	polygon->numVertices = 6;
	for(int i = 0; i < polygon->numVertices; i++) {
		polygon->vertices[i].z = 0;
		polygon->vertices[i].w = 1;
	}
	
	double characterAspect = 28.0/12.0; // macOs terminal
	
	int screenSizeX, screenSizeY;
	getmaxyx(stdscr, screenSizeY, screenSizeX);
	
	double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
	
	// Depth buffer
	DepthBuffer depthBuffer;
	depthBuffer.setSize(screenSizeX, screenSizeY);
	
	// Model
	Mat4D scaleMat = scaleMatrix(1, 1, 1);
	
	// View
	Mat4D cameraTranslation = translationMatrix(0, 0, -5);
	Coordinates3D cameraAxis = {0, 1, 0};
	cameraAxis = normalizeVector(cameraAxis);
	Mat4D cameraOrientation = rotationFromAngleAndUnitAxis(-M_PI_4, cameraAxis);
	Mat4D viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation );
	
	// Projection
	Mat4D projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, 100, 2);
	
	// Window
	Mat4D windowScale = scaleMatrix(screenSizeX/2, screenSizeY/2, 1);
	Mat4D translationScreen = translationMatrix(screenSizeX/2, screenSizeY/2, 0);
	Mat4D windowFull = matrixMultiply(translationScreen, windowScale);
	
	int numEdges;
	bool usePerspective = true;
	
	
	auto now = std::chrono::high_resolution_clock::now();
	auto before = now;
	int line;
	
	double angle = 0.0000;
	double modelAngle = 0;
	double distance =5;
	bool autoRotate = false;
	while (keepRunning == true) {
		line = 0;
		usleep(1000000.0/60.0);
		depthBuffer.reset();
		erase();
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
		mvprintw(line++, 0, "Time: %f", float_ms.count());
		/*
		 Make a model matrix
		 */
		Coordinates3D modelAxis = {0,1,0};
		Mat4D modelMatrix = rotationFromAngleAndUnitAxis(modelAngle, modelAxis);
		
		/*
		 Build the view/camera matrix:
		 */
		viewMatrix = translationMatrix(0, 0, -distance);
//		angle = 0.6;
		if (autoRotate)
			angle += 0.01;
		Coordinates3D axis = {0,0,1};
		Mat4D camRotation = rotationFromAngleAndUnitAxis(angle, axis);
		viewMatrix = matrixMultiply(viewMatrix,  camRotation);
		
		// Cube
		Coordinates4D verticesCopy[sizeof(vertices)/sizeof(vertices[0])];
		for (int i = 0; i < sizeof(vertices)/sizeof(vertices[0]); i++) {
			// Model
			verticesCopy[i] = matrixVectorMultiply(modelMatrix, vertices[i]);
			
			// View
			verticesCopy[i] = matrixVectorMultiply(viewMatrix, verticesCopy[i]);
			
			// Projection (final step)
			verticesCopy[i] = matrixVectorMultiply(projection, verticesCopy[i]);
		}
		numEdges = sizeof(edgeIndices)/sizeof(edgeIndices[0]);
		attron(COLOR_PAIR(3));
		attron(A_BOLD);
//		rasterizeTriangle(verticesCopy, edgeIndices, numEdges, windowFull, &depthBuffer, 'O', 1);
		attroff(A_BOLD);
		attroff(COLOR_PAIR(3));
		

		
		drawHorizonalLineWithDepthBuffer(48, 50, 0, '0', 0, 0, NULL);
//		drawHorizonalLineWithDepthBuffer(50, 70, 0, '1', 0, 0, NULL);
		Coordinates4D dummy;
		Coordinates3D dummy3D;
		drawHorizonalLineWithShader(50, 70, 0, 1, 1, dummy, dummy, dummy3D, dummy3D, NULL, &depthBuffer, defaultFragment);
		drawHorizonalLineWithDepthBuffer(70, 72, 0, '0', 0, 0, NULL);
		
		
		// Polygon
		int numPolys = sizeof(polygon)/sizeof(polygon[0]);
		Polygon4D polyCopy[numPolys];
		for (int p = 0; p < 1; p++) {
			polyCopy[p] = polygon[p];
			for (int i = 0; i < polyCopy[p].numVertices; i++) {
				// Model
				polyCopy[p].vertices[i] = matrixVectorMultiply(modelMatrix, polyCopy[p].vertices[i]);
				
				// View
				polyCopy[p].vertices[i] = matrixVectorMultiply(viewMatrix, polyCopy[p].vertices[i]);
				
				// Projection (final step)
				polyCopy[p].vertices[i] = matrixVectorMultiply(projection, polyCopy[p].vertices[i]);
			}
		}
		attron(COLOR_PAIR(2));
		attron(A_BOLD);
//		rasterizePolygon(polyCopy, numPolys, windowFull, &depthBuffer, '*', line);
		Mat4D polyModelView = matrixMultiply(viewMatrix, modelMatrix);
		rasterizePolygonsShader(polygon, numPolys, polyModelView, projection, windowFull, NULL, &depthBuffer, defaultFragment, line);
		attroff(A_BOLD);
		attroff(COLOR_PAIR(2));
		
		int ch;
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == KEY_RESIZE) {
			getmaxyx(stdscr, screenSizeY, screenSizeX);
			
			screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
			windowScale = scaleMatrix(screenSizeX/2, screenSizeY/2, 1);
			translationScreen = translationMatrix(screenSizeX/2, screenSizeY/2, 0);
			windowFull = matrixMultiply(translationScreen, windowScale);
			
			if (usePerspective) {
				projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, 100, 2);
			} else {
				projection = projectionMatrixOrtho(5*screenAspect, 5, 100, 0);
			}
			depthBuffer.setSize(screenSizeX, screenSizeY);
		} else if ( ch == KEY_LEFT) {
			angle += 0.05;
		} else if ( ch == KEY_RIGHT) {
			angle -= 0.05;
		} else if ( ch == KEY_UP) {
			distance -= 0.05;
		} else if ( ch == KEY_DOWN) {
			distance += 0.05;
		} else if ( ch == 'w' || ch == 'W') {
			modelAngle -= 0.05;
		} else if ( ch == 's' || ch == 'S') {
			modelAngle += 0.05;
		} else if ( ch == ' ' ) {
			autoRotate = !autoRotate;
		}

	}
	
	cleanupConsole();
	
	printf("angle = %f\n", angle);
	printf("distance = %f\n", distance);
	printf("modelAngle = %f\n", modelAngle);
	
	return 0;
};


