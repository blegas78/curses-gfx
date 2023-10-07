#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>

#include <chrono>

#include "curses-gfx.h"
#include "curses-clock.h"
#include "curses-gfx-3d.h"
#include "curses-gfx-handler.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}

typedef struct _LightParams {
	Coordinates4D modelView[10];
	Coordinates3D color[10];
	int numLights;
} LightParams;


void lightModelFs(const FragmentInfo& fInfo) {
	Coordinates3D* colorRGB = (Coordinates3D*)fInfo.data;
	//setRGB(fInfo.pixel, *colorRGB);
	
	Coordinates3D clippedRGB = clipRGB(*colorRGB);
	fInfo.colorOutput->r = clippedRGB.x*255;
	fInfo.colorOutput->g = clippedRGB.y*255;
	fInfo.colorOutput->b = clippedRGB.z*255;
	fInfo.colorOutput->a = 0;
}

void lightFs2(const FragmentInfo& fInfo) {
	LightParams* lights = (LightParams*)fInfo.data;
	
	Coordinates3D colorRGB = {0,0,0};
	
	Coordinates4D lightToFrag;
	Coordinates3D lightToFragNoramlized;
	double intensity2;
//	double lightMagnitudeSquared;
	double lightMagnitude;
	Coordinates3D lightReflected;
	double intensitySpecular;
	double intensity;
	
	for (int i = 0; i < lights->numLights; i++) {
		lightToFrag = vectorSubtract(lights->modelView[i], fInfo.location3D);
		if(dotProduct(lightToFrag, fInfo.normal) < 0) {
			continue;
		}
		lightToFragNoramlized = normalizeVectorFast(lightToFrag);
		intensity2 = 5.0*1.0;//dotProduct(lightToFragNoramlized, fInfo.normal) *0.75;
		
//		lightMagnitudeSquared = dotProduct(lightToFrag, lightToFrag) ;
		lightMagnitude = Q_rsqrt( dotProduct(lightToFrag, lightToFrag) );
		
		lightReflected = vectorScale(fInfo.normal, 2*dotProduct(fInfo.normal, lightToFragNoramlized));
		lightReflected = vectorSubtract(lightToFragNoramlized, lightReflected);
		
		intensitySpecular = dotProduct(lightReflected, fInfo.location3D);
//		intensitySpecular /= sqrt(dotProduct(fInfo.location3D, 	fInfo.location3D));
		intensitySpecular *= Q_rsqrt(dotProduct(fInfo.location3D, fInfo.location3D));
		intensitySpecular = pow(intensitySpecular, 32)*0.9;
		
//		intensity = (1/(lightMagnitudeSquared) * intensity2 +  intensitySpecular);
		intensity = (lightMagnitude*lightMagnitude * intensity2 +  intensitySpecular);
		
		
		colorRGB.x += intensity*lights->color[i].x;
		colorRGB.y += intensity*lights->color[i].y;
		colorRGB.z += intensity*lights->color[i].z;
	}
	
	
	
	//setRGB(fInfo.pixel, colorRGB);
	
	//	ColorRGBA result;
	Coordinates3D clippedRGB = clipRGB(colorRGB);
	fInfo.colorOutput->r = clippedRGB.x*255;
	fInfo.colorOutput->g = clippedRGB.y*255;
	fInfo.colorOutput->b = clippedRGB.z*255;
	fInfo.colorOutput->a = 0;
	//	fInfo.colorOutput = result;
}

void lightFs(const FragmentInfo& fInfo) {
	char c = ' ';
	
	Coordinates4D* lightModelView = (Coordinates4D*)fInfo.data;
	
	Coordinates4D lightToFrag = vectorSubtract(*lightModelView, fInfo.location3D);
	if(dotProduct(lightToFrag, fInfo.normal)< 0) {
		set(fInfo.pixel, c);
		return;
	}
	Coordinates3D lightToFragNoramlized = normalizeVector(lightToFrag);
	double intensity2 = 12.0;//dotProduct(lightToFragNoramlized, fInfo.normal) *0.75;
	
	double lightMagnitudeSquared = dotProduct(lightToFrag, lightToFrag) ;
	
	Coordinates3D lightReflected = vectorScale(fInfo.normal, 2*dotProduct(fInfo.normal, lightToFragNoramlized));
	lightReflected = vectorSubtract(lightToFragNoramlized, lightReflected);
	
	double intensitySpecular = dotProduct(lightReflected, fInfo.location3D);
	intensitySpecular /= sqrt(dotProduct(fInfo.location3D, 	fInfo.location3D));
	intensitySpecular = pow(intensitySpecular, 32)*0.9;
	
	double intensity = (1/(lightMagnitudeSquared) * intensity2 +  intensitySpecular);
	
	Coordinates3D colorRGB = {0,0,0};
	colorRGB.x = intensity*1.25 +0.1;
	colorRGB.y = intensity*0.9;
	colorRGB.z = intensity*3;
	
	//setRGB(fInfo.pixel, colorRGB);
	
//	ColorRGBA result;
	Coordinates3D clippedRGB = clipRGB(colorRGB);
	fInfo.colorOutput->r = clippedRGB.x*255;
	fInfo.colorOutput->g = clippedRGB.y*255;
	fInfo.colorOutput->b = clippedRGB.z*255;
	fInfo.colorOutput->a = 0;
//	fInfo.colorOutput = result;
	
}


int main(int argc, char** argv) {
	
    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();

	Coordinates4D cube[] = {
		{-1, -1, -1, 1},
		{-1, -1,  1, 1},
		{-1,  1,  1, 1},
		{-1,  1, -1, 1},
		{ 1, -1, -1, 1},
		{ 1, -1,  1, 1},
		{ 1,  1,  1, 1},
		{ 1,  1, -1, 1}
	};
	
	int edgeIndices[12][2] = {
		{0, 1},
		{1, 2},
		{2, 3},
		{3, 0},
		{0, 4},
		{1, 5},
		{2, 6},
		{3, 7},
		{4, 5},
		{5, 6},
		{6, 7},
		{7, 4}
	};
	
	int cubeTriangleIndices[][3] = {
		{0, 1, 2},	// left
		{2, 3, 0},
		{4, 6, 5},	// right
		{6, 4, 7},
		{6, 2, 5}, 	// top
		{2, 1, 5},
		{0, 3, 4}, 	// bottom
		{4, 3, 7},
		{2, 6, 3}, 	// front
		{7, 3, 6},
		{0, 4, 5}, 	// back
		{5, 1, 0}
	};
	int cubeQuadIndices[][4] = {
		{0, 1, 2, 3},	// left
		{7, 6, 5, 4},	// right
		{6, 2, 1, 5}, 	// top
		{0, 3, 7, 4}, 	// bottom
		{2, 6, 7, 3}, 	// front
		{0, 4, 5, 1} 	// back
	};
	
	Coordinates4D pyramid[] = {
		{-0.5, -0.5, 0, 1},
		{ 0.5, -0.5, 0, 1},
		{ 0.5,  0.5, 0, 1},
		{-0.5,  0.5, 0, 1},
		{0, 0, 1, 1}
	};
	int pyramidEdgeIndices[12][2] = {
		{0, 1},
		{1, 2},
		{2, 3},
		{3, 0},
		{0, 4},
		{1, 4},
		{2, 4},
		{3, 4}
	};
	
	Coordinates4D triangle[] = {
		{3, 0, 0, 1},
		{0, 3, 0, 1},
		{-3,  0, 0, 1},
		{0,  -3, 0, 1}
	};
	
	int triangleEdgeIndices[][3] = {
		{0, 1, 2},
		{3, 0, 2},
	};
	
#define NUM_GRID_LINES (11)
#define GRID_LINE_UNIT (2)
	Coordinates4D grid[NUM_GRID_LINES*2*2];	// 5x5 grid
	int gridEdgeIndices[NUM_GRID_LINES*2][2];
	for (int i = 0; i < NUM_GRID_LINES*2*2; i++) {
		
		grid[i].z = -1;
		grid[i].w = 1;
	}
	for (int i = 0; i < NUM_GRID_LINES; i++) {
		grid[i].x = -(NUM_GRID_LINES-1)/2 + i;
		grid[i].y = -(NUM_GRID_LINES-1)/2;
		
		grid[i+NUM_GRID_LINES].x = -(NUM_GRID_LINES-1)/2 + i;
		grid[i+NUM_GRID_LINES].y = (NUM_GRID_LINES-1)/2;
		
		grid[i+NUM_GRID_LINES*2].x = -(NUM_GRID_LINES-1)/2;
		grid[i+NUM_GRID_LINES*2].y = -(NUM_GRID_LINES-1)/2 + i;
		
		grid[i+NUM_GRID_LINES*3].x = (NUM_GRID_LINES-1)/2;
		grid[i+NUM_GRID_LINES*3].y = -(NUM_GRID_LINES-1)/2 + i;
		
		gridEdgeIndices[i][0] = i;
		gridEdgeIndices[i][1] = i+NUM_GRID_LINES;
		
		gridEdgeIndices[i+NUM_GRID_LINES][0] = i+NUM_GRID_LINES*2;
		gridEdgeIndices[i+NUM_GRID_LINES][1] = i+NUM_GRID_LINES*3;
	}
	
	for (int i = 0; i < NUM_GRID_LINES*2*2; i++) {
		grid[i].x *= GRID_LINE_UNIT;
		grid[i].y *= GRID_LINE_UNIT;
	}
	
	
	Coordinates4D floor[] = {
		{3, 0, 0, 1},
		{0, 3, 0, 1},
		{-3,  0, 0, 1},
		{0,  -3, 0, 1}
	};
	
	Polygon4D ground[5];
	ground[0].vertices[0].x = 10;
	ground[0].vertices[0].y = 10;
	ground[0].vertices[1].x = -10;
	ground[0].vertices[1].y = 10;
	ground[0].vertices[2].x = -10;
	ground[0].vertices[2].y = -10;
	ground[0].vertices[3].x = 10;
	ground[0].vertices[3].y = -10;
	for (int p = 0; p < 5; p++) {
		ground[p].numVertices = 4;
		for (int i = 0; i < 4; i++) {
			ground[p].vertices[i].z = 0;
			ground[p].vertices[i].w = 1;
//			ground[p].normals[i].x = 0;
//			ground[p].normals[i].y = 0;
//			ground[p].normals[i].z = 1;
			
			Coordinates3D hsv, rgb;
			hsv.x = (i*5+p)*360/20;
			hsv.y = 0.5;
			hsv.z = 0.5;
			rgb = hslToRgb(hsv);
			
			ground[p].colors[i].r = rgb.x;
			ground[p].colors[i].g = rgb.y;
			ground[p].colors[i].b = rgb.z;
			ground[p].colors[i].a = 0;
		}
	}
	Coordinates4D wallHeight = {0,0,10,0};
	for (int p = 1; p < 5; p++) {
		ground[p].vertices[0] = ground[0].vertices[p%4];
		ground[p].vertices[1] = vectorAdd(ground[0].vertices[p%4], wallHeight);
		ground[p].vertices[2] = vectorAdd(ground[0].vertices[(p+1)%4], wallHeight);
		ground[p].vertices[3] = ground[0].vertices[(p+1)%4];
		
//		Coordinates4D normal = crossProduct(vectorSubtract(ground[p].vertices[1], ground[p].vertices[0]), vectorSubtract(ground[p].vertices[2], ground[p].vertices[0]));
//		ground[p].normals[0].x = normal.x;
//		ground[p].normals[0].y = normal.y;
//		ground[p].normals[0].z = normal.z;
//		ground[p].normals[0] = normalizeVector(ground[p].normals[0] );
//		ground[p].normals[1] = ground[p].normals[0];
//		ground[p].normals[2] = ground[p].normals[0];
//		ground[p].normals[3] = ground[p].normals[0];
	}
	fillPolygonNormals(ground, sizeof(ground)/sizeof(ground[0]));
	
	
	
	
	double characterAspect = 28.0/12.0; // macOs terminal
//	double characterAspect = 28.0/14.0; // raspbian terminal
//	double characterAspect = 6.0/4.0; // zipitZ2
	
	int screenSizeX, screenSizeY;
	getmaxyx(stdscr, screenSizeY, screenSizeX);
	
#ifdef FB_SUPPORT	// HACK
	characterAspect = 1.0;	// RGB panel
#define PANEL_X_RES (64)
#define PANEL_Y_RES (64)
	screenSizeY = PANEL_Y_RES;
	screenSizeX = PANEL_X_RES;
#endif
	
	RenderPipeline mRenderPipeline;
	mRenderPipeline.resize(screenSizeX, screenSizeY);
	
	double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
	
	// Depth buffer
//	DepthBuffer depthBuffer;
//	depthBuffer.setSize(screenSizeX, screenSizeY);
	
	// Model
	Mat4D scaleMat = scaleMatrix(1, 1, 1);
	
	// View
	Mat4D cameraTranslation = translationMatrix(0, 0, -5);
	Coordinates3D cameraAxis = {0, 1, 0};
	cameraAxis = normalizeVector(cameraAxis);
	Mat4D cameraOrientation = rotationFromAngleAndUnitAxis(-M_PI_4, cameraAxis);
	Mat4D viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation );
	
	// Projection
	double zFar = 30;
	double zNear = 0.1;
	Mat4D projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, zFar, zNear);
	
	// Viewport
//	Mat4D windowScale = scaleMatrix((double)screenSizeX/2, (double)screenSizeY/2, 1);
//	Mat4D translationScreen = translationMatrix((double)screenSizeX/2 -0.5, (double)screenSizeY/2 -0.5, 0);
//	Mat4D windowFull = matrixMultiply(translationScreen, windowScale);
	Mat4D windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
	
	// Light
	Coordinates4D light[3];// = {3, 3, 3, 1};
	Coordinates4D lightModelView = {3, 3, 3, 1};
	
	auto now = std::chrono::high_resolution_clock::now();
	auto before = now;
	
	int debugLine = 0;
	int numEdges;
	double lightAngle = 0;
	double angle = 0;
	double cube2angle = 0;
	double tilt = M_PI/4;
	bool usePerspective = true;
	bool showGrid = false;
	bool autoRotate = true;
	bool showDepth = false;
	double delayTime = 0;//1.0/60;
	while (keepRunning == true) {
		debugLine = 0;
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
		
		delayTime += 0.001*(1.0/60.0 - float_ms.count()/1000.0);
        if(delayTime > 0 && delayTime < 1.0/60.0){
            usleep(1000000.0*delayTime);
        }
//		depthBuffer.reset();
		mRenderPipeline.reset();
		//erase();
		mvprintw(debugLine++, 0, "FPS: %f", 1000.0/float_ms.count());
		mvprintw(debugLine++, 0, "Delay time %f", delayTime);
//		mvprintw(debugLine++, 0, "tilt: %f", tilt*180.0/M_PI);
		
		
		
		/*
		 Build the view/camera matrix:
		 */
		cameraAxis.x = 0;
		cameraAxis.y = 0;
		cameraAxis.z = 1;
		cameraOrientation = rotationFromAngleAndUnitAxis(angle, cameraAxis);
		cameraOrientation = transpose(cameraOrientation);
		
		double distance = 2*sin(angle/2);
		cameraTranslation = translationMatrix(-(5+distance)*sin(angle), (5+distance) * cos(angle), -(2*cos(angle/5)+distance)-4);

		viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation);
		
		cameraAxis.x = 1;
		cameraAxis.y = 0;
		cameraAxis.z = 0;
		cameraOrientation = rotationFromAngleAndUnitAxis(tilt, cameraAxis);
		cameraOrientation = transpose(cameraOrientation);

		viewMatrix = matrixMultiply( cameraOrientation, viewMatrix);
		
		lightAngle+=0.01;
		light[0].x = 6*cos(lightAngle);
		light[0].y = 6*sin(lightAngle);
		light[0].z = 3 + 1.0*sin(lightAngle*5.0);
		light[0].w = 1;
		lightModelView = matrixVectorMultiply(viewMatrix, light[0]);
		
		// Light Model
		LightParams mLightParams;
		mLightParams.numLights = 3;
		mLightParams.modelView[0] = lightModelView;
		mLightParams.color[0].x = 0.0;//1.25/3.0;
		mLightParams.color[0].y = 0.0;//0.9/3.0;
		mLightParams.color[0].z = 3.0/3.0;
		
		
		light[1].x = 6*cos(lightAngle*1.25);
		light[1].y = 6*sin(lightAngle*1.25);
		light[1].z = 3 + 1.0*sin(lightAngle*6.0);
		light[1].w = 1;
		mLightParams.modelView[1] = matrixVectorMultiply(viewMatrix, light[1]);
		mLightParams.color[1].x = 1.0;
		mLightParams.color[1].y = 0.0;//0.1;
		mLightParams.color[1].z = 0.0;
		
		light[2].x = 6*cos(lightAngle*1.75);
		light[2].y = 6*sin(lightAngle*1.75);
		light[2].z = 3 + 1.0*sin(lightAngle*7.0);
		light[2].w = 1;
		mLightParams.modelView[2] = matrixVectorMultiply(viewMatrix, light[2]);
		mLightParams.color[2].x = 0.0;
		mLightParams.color[2].y = 1.0;
		mLightParams.color[2].z = 0.0;
		
		for (int i = 0; i < mLightParams.numLights; i++) {
			Mat4D lightScale = scaleMatrix(0.15, 0.15, 0.15);
			Mat4D lightTranslation = translationMatrix(light[i].x, light[i].y, light[i].z);
			Mat4D lightModel = matrixMultiply(lightTranslation, lightScale);
			Mat4D lightCubeModelView = matrixMultiply(viewMatrix, lightModel);
			numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
			
			lightModelView = matrixVectorMultiply(viewMatrix, light[i]);
//			rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, windowFull, (void*)&mLightParams.color[i], &depthBuffer, lightModelFs, debugLine);
			mRenderPipeline.setFragmentShader(lightModelFs);

			mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, windowFull, (void*)&mLightParams.color[i], debugLine);
		}
		
//		rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, windowFull, (void*)&mLightParams, &depthBuffer, lightFs2, debugLine);
		
//		// Cube
//		Coordinates4D cubeRotated[sizeof(cube)/sizeof(cube[0])];
//		for (int i = 0; i < sizeof(cube)/sizeof(cube[0]); i++) {
//			// Copy
//			cubeRotated[i] = cube[i];
//
//			// Model
//
//			// View
//			cubeRotated[i] = matrixVectorMultiply(viewMatrix, cubeRotated[i]);
//
//			// Projection (final step)
//			cubeRotated[i] = matrixVectorMultiply(projection, cubeRotated[i]);
//		}
//		numEdges = sizeof(edgeIndices)/sizeof(edgeIndices[0]);
//		attron(COLOR_PAIR(3));
//		attron(A_BOLD);
//		rasterize(cubeRotated, edgeIndices, numEdges, windowFull, &depthBuffer);
//		attroff(A_BOLD);
//		attroff(COLOR_PAIR(3));
		
		
		// Cube 2
		cube2angle += 0.02;
		Coordinates3D cube2roationAxis = {1+sin(0.9*cube2angle), sin(0.8*cube2angle), sin(0.7*cube2angle)};
		cube2roationAxis = normalizeVector(cube2roationAxis);
		Mat4D cube2rotation = rotationFromAngleAndUnitAxis(cube2angle*1.1, cube2roationAxis);
		Mat4D cube2translation = translationMatrix(1, -1, 0.1);
		Mat4D modelViewCube2 = matrixMultiply(cube2translation, cube2rotation);
		modelViewCube2 = matrixMultiply(viewMatrix, modelViewCube2);
//		for (int i = 0; i < sizeof(cube)/sizeof(cube[0]); i++) {
//			// Copy
//			cubeRotated[i] = cube[i];
//
//			// Model
//			cubeRotated[i] = matrixVectorMultiply(cube2rotation, cube[i]);
//			cubeRotated[i] = matrixVectorMultiply(cube2translation, cubeRotated[i]);
//
//			// View
//			cubeRotated[i] = matrixVectorMultiply(viewMatrix, cubeRotated[i]);
//
//			// Projection (final step)
//			cubeRotated[i] = matrixVectorMultiply(projection, cubeRotated[i]);
//		}
//		attron(COLOR_PAIR(2));
//		attron(A_BOLD);
//		rasterize(cubeRotated, edgeIndices, numEdges, windowFull, &depthBuffer);
//		numEdges = sizeof(cubeTriangleIndices)/sizeof(cubeTriangleIndices[0]);
//		rasterizeTriangle(cubeRotated, cubeTriangleIndices, numEdges, windowFull, &depthBuffer, '`', 12);
		numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
//		rasterizeQuads(cubeRotated, cubeQuadIndices, numEdges, windowFull, &depthBuffer, '`', debugLine);
//		rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewCube2, projection, windowFull, NULL, &depthBuffer, defaultFragment, debugLine);
		mRenderPipeline.setFragmentShader(defaultFragment);
		mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewCube2, projection, windowFull, NULL, debugLine);
//		attroff(A_BOLD);
//		attroff(COLOR_PAIR(2));
		
//		// Pyramid
//		Coordinates4D pyramid2[sizeof(pyramid)/sizeof(pyramid[0])];
////		Mat4D pyramidTranslation = translationMatrix(0, 0, sin(cube2angle*2));
//		Mat4D pyramidTranslation = translationMatrix(0, 0, -1);
//		Mat4D pyramidScale = scaleMatrix(1.5, 1.5, 1.5);
//		for (int i = 0; i < sizeof(pyramid2)/sizeof(pyramid2[0]); i++) {
//			// Copy
//			pyramid2[i] = pyramid[i];
//
//			// Model
//			pyramid2[i] = matrixVectorMultiply(pyramidScale, pyramid2[i]);
//			pyramid2[i] = matrixVectorMultiply(pyramidTranslation, pyramid2[i]);
//
//			// View
//			pyramid2[i] = matrixVectorMultiply(viewMatrix, pyramid2[i]);
//
//			// Projection (final step)
//			pyramid2[i] = matrixVectorMultiply(projection, pyramid2[i]);
//		}
//		numEdges = sizeof(pyramidEdgeIndices)/sizeof(pyramidEdgeIndices[0]);
//		attron(COLOR_PAIR(1));
//		attron(A_BOLD);
//		rasterize(pyramid2, pyramidEdgeIndices, numEdges, windowFull, &depthBuffer);
//		attroff(A_BOLD);
//		attroff(COLOR_PAIR(1));
//
//
//		// Grid
//		if (showGrid) {
//			Coordinates4D grid2[sizeof(grid)/sizeof(grid[0])];
//			for (int i = 0; i < sizeof(grid2)/sizeof(grid2[0]); i++) {
//				grid2[i] = grid[i];
//				// Model
//
//				// View
//				grid2[i] = matrixVectorMultiply(viewMatrix, grid2[i]);
//
//				// Projection
//				grid2[i] = matrixVectorMultiply(projection, grid2[i]);
//			}
//			numEdges = sizeof(gridEdgeIndices)/sizeof(gridEdgeIndices[0]);
//			rasterize(grid2, gridEdgeIndices, numEdges, windowFull, &depthBuffer);
//		}
		
//		// triangle
//		Coordinates4D triangle2[sizeof(triangle)/sizeof(triangle[0])];
//		for (int i = 0; i < sizeof(triangle2)/sizeof(triangle2[0]); i++) {
//			// Model
//			triangle2[i] = matrixVectorMultiply(cube2translation, triangle[i]);
//
//			// View
//			triangle2[i] = matrixVectorMultiply(viewMatrix, triangle2[i]);
//
//			// Projection
//			triangle2[i] = matrixVectorMultiply(projection, triangle2[i]);
//		}
//		numEdges = sizeof(triangleEdgeIndices)/sizeof(triangleEdgeIndices[0]);
//		rasterizeTriangle(triangle2, triangleEdgeIndices, numEdges, windowFull, &depthBuffer, '`');
		
		// triangle
//		Coordinates4D cubeTriangle[sizeof(triangle)/sizeof(triangle[0])];
		Mat4D solidCubeScale = scaleMatrix(1.05, 1.05, 1.05);
		Coordinates3D solidAxis = {0,0,1};
		Mat4D solidRotation = rotationFromAngleAndUnitAxis(M_PI_4, solidAxis);
//		for (int i = 0; i < sizeof(cube)/sizeof(cube[0]); i++) {
//			// Model
//
//			cubeRotated[i] = matrixVectorMultiply(solidCubeScale, cube[i]);
//			cubeRotated[i] = matrixVectorMultiply(solidRotation, cube[i]);
//			cubeRotated[i] = matrixVectorMultiply(cube2translation, cubeRotated[i]);
//
//			// View
//			cubeRotated[i] = matrixVectorMultiply(viewMatrix, cubeRotated[i]);
//
//			// Projection
//			cubeRotated[i] = matrixVectorMultiply(projection, cubeRotated[i]);
//		}
////		numEdges = sizeof(cubeTriangleIndices)/sizeof(cubeTriangleIndices[0]);
////		rasterizeTriangle(cubeRotated, cubeTriangleIndices, numEdges, windowFull, &depthBuffer, ' ', 0);
//		numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
//		rasterizeQuads(cubeRotated, cubeQuadIndices, numEdges, windowFull, &depthBuffer, ' ', debugLine);
		
		Mat4D modelViewBlackCube = matrixMultiply(solidRotation, solidCubeScale);
		modelViewBlackCube = matrixMultiply( cube2translation, modelViewBlackCube);
		modelViewBlackCube = matrixMultiply(viewMatrix, modelViewBlackCube);
		numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
//		rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&lightModelView, &depthBuffer, lightFs, debugLine);
//		rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, &depthBuffer, lightFs2, debugLine);
		mRenderPipeline.setFragmentShader(lightFs2);
		mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, debugLine);
		
		
		// Another cube in the middle:
		Mat4D middleCubeScale = scaleMatrix(0.75, 0.75, 0.75);
		Coordinates3D middleCubeAxis = {1+sin(0.82*cube2angle), sin(0.93*cube2angle), sin(01.41*cube2angle)};
		middleCubeAxis = normalizeVector(middleCubeAxis);
		Mat4D middleCubeRotation = rotationFromAngleAndUnitAxis(cube2angle, middleCubeAxis);
		Mat4D middleCubeTranslation = translationMatrix(0, 0, 3);
		Mat4D modelViewMiddleCube = matrixMultiply( middleCubeRotation, middleCubeScale);
		modelViewMiddleCube = matrixMultiply( middleCubeTranslation, modelViewMiddleCube);
		modelViewMiddleCube = matrixMultiply( viewMatrix, modelViewMiddleCube);
		numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
		mRenderPipeline.setFragmentShader(lightFs2);
		mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewMiddleCube, projection, windowFull, (void*)&mLightParams, debugLine);
		
		
		// ground:
		if (!showGrid) {
//			Polygon4D groundCopy;
//			groundCopy.numVertices = ground.numVertices;
//			for (int i = 0; i < groundCopy.numVertices; i++) {
//				// View
//				groundCopy.vertices[i] = matrixVectorMultiply(viewMatrix, ground.vertices[i]);
//				groundCopy.normals[i] = matrixVectorMultiply(viewMatrix, ground.normals[i]);
//
////				 Projection
////				groundCopy.vertices[i] = matrixVectorMultiply(projection, groundCopy.vertices[i]);
//			}
			Mat4D groundTranslation = translationMatrix(0, 0, -1);
			Mat4D groundModelView = matrixMultiply( viewMatrix, groundTranslation);
//			rasterizePolygon(&groundCopy, 1, windowFull, &depthBuffer, ' ', debugLine);
//			attron(COLOR_PAIR(0));
			debugLine += 5;
//			rasterizePolygonsShader(&ground, 1, groundModelView, projection, windowFull, &depthBuffer, defaultFragment, debugLine);
//			rasterizePolygonsShader(ground, 5, groundModelView, projection, windowFull, (void*)&lightModelView, &depthBuffer, lightFs, debugLine);
//			rasterizePolygonsShader(ground, 5, groundModelView, projection, windowFull, (void*)&mLightParams, &depthBuffer, lightFs2, debugLine);
			mRenderPipeline.setFragmentShader(lightFs2);
			
			mRenderPipeline.rasterizePolygonsShader(ground, 5, groundModelView, projection, windowFull, (void*)&mLightParams, debugLine);
//			rasterizePolygonsShader(&ground, 1, groundModelView, projection, windowFull, NULL, &depthBuffer, defaultFragment, debugLine);
//			attroff(COLOR_PAIR(0));
		}
		
		
//		Coordinates4D p;
//		Coordinates3D n;
//		drawHorizonalLineWithShader(30, 210, 1, 0.4, 0.01, p, p, n, n, modelViewBlackCube, NULL, &depthBuffer, defaultFragment);
		
		if (showDepth) {
			mRenderPipeline.depthBufferToTerminal();
		} else {
			mRenderPipeline.renderBufferToTerminal();
		}
		
		if (autoRotate) {
			angle -= 0.01;
		}

		int ch;
//		refresh();
//		continue;
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == KEY_RESIZE) {
#ifdef FB_SUPPORT	// HACK
			//screenSizeY = 64;
			//screenSizeX = 64;
#else
			getmaxyx(stdscr, screenSizeY, screenSizeX);
#endif
			
			mRenderPipeline.resize(screenSizeX, screenSizeY);
			
			screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
			windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
			
			if (usePerspective) {
				projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, zFar, zNear);
			} else {
				projection = projectionMatrixOrtho(5*screenAspect, 5, zFar, zNear);
			}
//			depthBuffer.setSize(screenSizeX, screenSizeY);
		} else if (ch == 'o' || ch == 'O') {
			usePerspective	= !usePerspective;
			if (usePerspective) {
				projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, zFar, zNear);
			} else {
				projection = projectionMatrixOrtho(5*screenAspect, 5, zFar, zNear);
			}
		} else if ( ch == 'g' || ch == 'G') {
			showGrid = !showGrid;
		} else if ( ch == KEY_LEFT) {
			angle -= 0.05;
		} else if ( ch == KEY_RIGHT) {
			angle += 0.05;
		} else if ( ch == KEY_UP) {
			tilt += 0.05;
		} else if ( ch == KEY_DOWN) {
			tilt -= 0.05;
		} else if ( ch == ' ' ) {
			autoRotate = !autoRotate;
		} else if ( ch == 'd' || ch == 'D' ) {
			showDepth = !showDepth;
		}

	}
	
	return 0;
};


