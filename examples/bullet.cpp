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
#include "curses-gfx-texture.h"
#include "curses-gfx-physics.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}

Texture testTexture(10,10);

typedef struct _LightParams {
	Coordinates4D modelView[10];
	Coordinates3D color[10];
	int numLights;
    Coordinates3D* colorModel;
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


typedef struct _UniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
    Mat4D normalMatrix;
} UniformInfo;



template <class T, class U> void myVertexShader(U* uniformInfo, T& output, const T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
    output.normal = matrixVectorMultiply(uniformInfo->normalMatrix, input.normal);
    output.color = input.color;
    output.textureCoord = input.textureCoord;
}

typedef struct _CubeVertexInfo {
    Coordinates4D vertex;
    Coordinates4D location;
    Coordinates3D normal;
    Coordinates2Df textureCoord;
    ColorRGBA color;
} CubeVertexInfo;

REGISTER_VERTEX_LAYOUT(CubeVertexInfo)
    MEMBER(location),
    MEMBER(normal),
    MEMBER(textureCoord),
    MEMBER(color)
END_VERTEX_LAYOUT(CubeVertexInfo)


void lightFs3(const FragmentInfo& fInfo) {
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    
//    *fInfo.colorOutput = vertexInfo->color;
    
//    return;
    
    LightParams* lights = (LightParams*)fInfo.data;
//    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    
    Coordinates3D colorRGB = {0,0,0};
    
    Coordinates4D lightToFrag;
    Coordinates3D lightToFragNoramlized;
    double intensity2;
//    double lightMagnitudeSquared;
    double lightMagnitude;
    Coordinates3D lightReflected;
    double intensitySpecular;
    double intensity;
    
    Coordinates3D vNormal = normalizeVectorFast(vertexInfo->normal);
    Coordinates3D viewDir = normalizeVectorFast(vertexInfo->location);
    
    for (int i = 0; i < lights->numLights; i++) {
//        lightToFrag = vectorSubtract(lights->modelView[i], fInfo.location3D);
        lightToFrag = vectorSubtract(lights->modelView[i], vertexInfo->location);
        if(dotProduct(lightToFrag, vNormal) < 0) {
            continue;
        }
        lightToFragNoramlized = normalizeVectorFast(lightToFrag);
        intensity2 = 5.0*1.0;//dotProduct(lightToFragNoramlized, fInfo.normal) *0.75;
        
//        lightMagnitudeSquared = dotProduct(lightToFrag, lightToFrag) ;
        lightMagnitude = Q_rsqrt( dotProduct(lightToFrag, lightToFrag) );
        
        lightReflected = vectorScale(vNormal, 2*dotProduct(vNormal, lightToFragNoramlized));
        lightReflected = vectorSubtract(lightToFragNoramlized, lightReflected);
        
        intensitySpecular = dotProduct(lightReflected, vertexInfo->location);
//        intensitySpecular /= sqrt(dotProduct(fInfo.location3D,     fInfo.location3D));
//        intensitySpecular *= Q_rsqrt(dotProduct(vertexInfo->location, vertexInfo->location));
        
        lightReflected = normalizeVectorFast(lightReflected);
        intensitySpecular = (fmax(0,dotProduct(viewDir, lightReflected)));
        intensitySpecular = pow(intensitySpecular, 32)*0.9;
        
//        intensity = (1/(lightMagnitudeSquared) * intensity2 +  intensitySpecular);
        intensity = (lightMagnitude*lightMagnitude * intensity2 + lightMagnitude) +  intensitySpecular;
        intensity *= 1;
        
        colorRGB.x += intensity*lights->color[i].x;
        colorRGB.y += intensity*lights->color[i].y;
        colorRGB.z += intensity*lights->color[i].z;
    }
    
    
    
    //setRGB(fInfo.pixel, colorRGB);
    
    //    ColorRGBA result;
//    colorRGB.x *= 0.5;
//    colorRGB.y *= 0.5;
//    colorRGB.z *= 0.5;
//    colorRGB.x += (double)vertexInfo->color.r * (1.0/255.0 * 0.5);
//    colorRGB.y += (double)vertexInfo->color.g * (1.0/255.0 * 0.5);
//    colorRGB.z += (double)vertexInfo->color.b * (1.0/255.0 * 0.5);
    Coordinates3D clippedRGB = clipRGB(colorRGB * *lights->colorModel + *lights->colorModel*0.15);
    fInfo.colorOutput->r = clippedRGB.x*255;
    fInfo.colorOutput->g = clippedRGB.y*255;
    fInfo.colorOutput->b = clippedRGB.z*255;
    fInfo.colorOutput->a = 0;
    
    //    fInfo.colorOutput = result;
}

void lightFs4(const FragmentInfo& fInfo) {
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    
//    *fInfo.colorOutput = vertexInfo->color;
    
//    return;
    
    LightParams* lights = (LightParams*)fInfo.data;
//    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    
    Coordinates3D colorRGB = {0,0,0};
    
    Coordinates4D lightToFrag;
    Coordinates3D lightToFragNoramlized;
    double intensity2;
//    double lightMagnitudeSquared;
    double lightMagnitude;
    Coordinates3D lightReflected;
    double intensitySpecular;
    double intensity;
    
    Coordinates3D vNormal = normalizeVectorFast(vertexInfo->normal);
    
    for (int i = 0; i < lights->numLights; i++) {
//        lightToFrag = vectorSubtract(lights->modelView[i], fInfo.location3D);
        lightToFrag = vectorSubtract(lights->modelView[i], vertexInfo->location);
        if(dotProduct(lightToFrag, vNormal) < 0) {
            continue;
        }
        lightToFragNoramlized = normalizeVectorFast(lightToFrag);
        intensity2 = 5.0*1.0;//dotProduct(lightToFragNoramlized, fInfo.normal) *0.75;
        
//        lightMagnitudeSquared = dotProduct(lightToFrag, lightToFrag) ;
        lightMagnitude = Q_rsqrt( dotProduct(lightToFrag, lightToFrag) );
        
        lightReflected = vectorScale(vNormal, 2*dotProduct(vNormal, lightToFragNoramlized));
        lightReflected = vectorSubtract(lightToFragNoramlized, lightReflected);
        
        intensitySpecular = dotProduct(lightReflected, vertexInfo->location);
//        intensitySpecular /= sqrt(dotProduct(fInfo.location3D,     fInfo.location3D));
        intensitySpecular *= Q_rsqrt(dotProduct(vertexInfo->location, vertexInfo->location));
        intensitySpecular = pow(intensitySpecular, 32)*0.9;
        
//        intensity = (1/(lightMagnitudeSquared) * intensity2 +  intensitySpecular);
        intensity = (lightMagnitude*lightMagnitude * intensity2 +  intensitySpecular);
        
        
        colorRGB.x += intensity*lights->color[i].x;
        colorRGB.y += intensity*lights->color[i].y;
        colorRGB.z += intensity*lights->color[i].z;
    }
    
    
    
    //setRGB(fInfo.pixel, colorRGB);
    
    //    ColorRGBA result;
    colorRGB.x *= 0.85;
    colorRGB.y *= 0.85;
    colorRGB.z *= 0.85;
    colorRGB.x += (double)testTexture.sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y).r * (1.0/255.0 * 0.5);
    colorRGB.y += (double)testTexture.sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y).g * (1.0/255.0 * 0.5);
    colorRGB.z += (double)testTexture.sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y).b * (1.0/255.0 * 0.5);
    Coordinates3D clippedRGB = clipRGB(colorRGB);
    fInfo.colorOutput->r = clippedRGB.x*255;
    fInfo.colorOutput->g = clippedRGB.y*255;
    fInfo.colorOutput->b = clippedRGB.z*255;
    fInfo.colorOutput->a = 0;
    
    //    fInfo.colorOutput = result;
}

float randZeroToOne() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

int main(int argc, char** argv) {

    const int numCubes = 100;
    Mat4D cubeInitial[numCubes];
    Coordinates3D cubeScale[numCubes];
    double cubeMass[numCubes];
    btRigidBody* cubePhysics[numCubes];
    Coordinates3D cubeColor[numCubes];
    
    for(int i = 0; i < numCubes; i++) {
        
        Coordinates3D axis = {randZeroToOne(), randZeroToOne(), randZeroToOne()};
        cubeColor[i] = axis;
//        if(dotProduct(cubeColor[i],cubeColor[i]) < 0.2) {
//            cubeColor[i] = normalizeVector(cubeColor[i]);
//        }
        axis = normalizeVector(axis);
        Mat4D initialCubeRotation = rotationFromAngleAndUnitAxis(i*40, axis);
        cubeInitial[i] = translationMatrix(0, 0, 5*(i+1));
        cubeInitial[i] = matrixMultiply( cubeInitial[i], initialCubeRotation);
//        initialCubeRotation = rotationFromAngleAndUnitAxis(348923, {sqrt(2), 0, sqrt(2)});
//        cubeInitial[1] = translationMatrix(0, 0, 20);
//        cubeInitial[1] = matrixMultiply( cubeInitial[1], initialCubeRotation);
        double scale = i*1.0/(double)numCubes +0.15;
        cubeScale[i] = {scale, scale, scale};
//        cubeScale[1] = 1.5;
//        cubeMass[0] = 1;
        cubeMass[i] = 1.0 * pow(scale, 3);
        
        
    }
    
    const int numStatic = 5;
    btRigidBody* staticCollision[numStatic];
    Coordinates3D staticScale[numStatic];
    Mat4D staticMatrix[numStatic];
    Coordinates3D staticColor[numStatic];
    
    staticColor[0] = {0.75, 0.75, 0.75};
    staticScale[0] = {10, 10, 1};
    for(int i = 1; i < numStatic; i++) {
        staticColor[i] = {0.1, 0.1, 1};
        staticColor[i] = normalizeVector(staticColor[i]);
    }
    
//    staticCollision[0] = mCursesGfxPhysics.addCube(10, 0, cubeInitial[i]);
    staticMatrix[0] = translationMatrix(0, 0, -1);
    staticMatrix[1] = translationMatrix(-10.1, 0, 1);
    staticScale[1] = {.1, 10, 5};
    staticMatrix[2] = translationMatrix(10.1, 0, 1);
    staticScale[2] = {.1, 10, 5};
    staticMatrix[3] = translationMatrix(0, -10.1, 1);
    staticScale[3] = {10, .1, 5};
    staticMatrix[4] = translationMatrix(0, 10.1, 1);
    staticScale[4] = {10, .1, 5};
    
    
    CursesGfxPhysics mCursesGfxPhysics;
    printf("mCursesGfxPhysics.initialize()...\n");
    mCursesGfxPhysics.initialize();
//    printf("mCursesGfxPhysics.addBody()...\n");
//    mCursesGfxPhysics.addBody();
    
//    printf("mCursesGfxPhysics.addSphere()...\n");
//    mCursesGfxPhysics.addSphere();
    for(int i = 0; i < numCubes; i++) {
//        printf("mCursesGfxPhysics.addCube()... %d\n", i);
        cubePhysics[i] = mCursesGfxPhysics.addCube(cubeScale[i], cubeMass[i], cubeInitial[i]);
    }
    
    for(int i = 0; i < numStatic; i++) {
//            printf("mCursesGfxPhysics.addCube()... %d\n", i);
        staticCollision[i] = mCursesGfxPhysics.addCube(staticScale[i], 0, staticMatrix[i]);
    }
//    printf("mCursesGfxPhysics.addCube() again...\n");
//    initial = translationMatrix(0, 0, 20);
//    double scaleCube2 = 1;
//    btRigidBody* cubePhysics2 = mCursesGfxPhysics.addCube(scaleCube2, 4.0, initial);
//    printf("mCursesGfxPhysics.update()...\n");
//    mCursesGfxPhysics.update(0, 1);
	
    for(int i = 0; i < testTexture.width; i++) {
        for(int j = 0; j < testTexture.height; j++) {
            ColorRGBA color = {255,255,255,0};
            testTexture.set(((double)i)/(double)testTexture.width, ((double)j)/(double)testTexture.height, color);
        }
    }
    
    bool skip = true;
    for(int i = 0; i < testTexture.width; i++) {
        for(int j = 0; j < testTexture.height; j++) {

            if(skip) {
                skip = false;
                ColorRGBA color = {0,0,0,0};
                testTexture.set(((double)i)/(double)testTexture.width, ((double)j)/(double)testTexture.height, color);
            } else {
                skip = true;
            }
        }
        skip = !skip;
    }
    
    
    testTexture.set(0.2, 0.2, {255,0,0,0});

    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();
    
    UniformInfo mUniformInfo;

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
	
    CubeVertexInfo cubeVi[12*3];
    int cubeViIndices[12][3];
    for(int t = 0; t < 12; t++) {   // for each triangle in a cube
        for(int v = 0; v < 3; v++) {
            cubeVi[v + 3*t].vertex = cube[cubeTriangleIndices[t][v]];
            cubeViIndices[t][v] = v + 3*t;
        }
    }
    for(int i = 0; i < 6; i++) {    // 6 sides
        cubeVi[i+0].normal =  {-1, 0, 0};     // left
        cubeVi[i+6].normal =  { 1, 0, 0};     // right
        cubeVi[i+12].normal = { 0, 0, 1};     // top
        cubeVi[i+18].normal = { 0, 0,-1};     // bottom
        cubeVi[i+24].normal = { 0, 1, 0};     // front
        cubeVi[i+30].normal = { 0,-1, 0};     // back
    }
    
    
    CubeVertexInfo squareVi[4];
    int squareViIndices[2][3];
    squareVi[0].vertex = {-1, -1, 0, 1};
    squareVi[1].vertex = { 1, -1, 0, 1};
    squareVi[2].vertex = { 1,  1, 0, 1};
    squareVi[3].vertex = {-1,  1, 0, 1};
    
    squareVi[0].textureCoord = {0, 0};
    squareVi[1].textureCoord = {0, 1};
    squareVi[2].textureCoord = {1, 1};
    squareVi[3].textureCoord = {1, 0};
    
    squareVi[0].color = {255, 0, 0, 0};
    squareVi[1].color = {0, 255, 0, 0};
    squareVi[2].color = {0, 0, 255, 0};
    squareVi[3].color = {0, 0, 0, 0};
    
    squareViIndices[0][0] = 0;  // right handed
    squareViIndices[0][1] = 1;
    squareViIndices[0][2] = 2;
    squareViIndices[1][0] = 0;  // left handed
    squareViIndices[1][1] = 2;
    squareViIndices[1][2] = 3;
    
    for(int i = 0; i < 4; i++) {
        squareVi[i].normal = {0,0,1};
    }
	
	double characterAspect = 28.0/12.0; // macOs terminal
//	double characterAspect = 28.0/14.0; // raspbian terminal
//	double characterAspect = 6.0/4.0; // zipitZ2
	
	int screenSizeX, screenSizeY;
	getmaxyx(stdscr, screenSizeY, screenSizeX);
	
//    RasterizerThreadPool::setRenderThreadCount(4);
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
	double zFar = 1000;
	double zNear = 0.1;
	Mat4D projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, zFar, zNear);
	
	// Viewport
//	Mat4D windowScale = scaleMatrix((double)screenSizeX/2, (double)screenSizeY/2, 1);
//	Mat4D translationScreen = translationMatrix((double)screenSizeX/2 -0.5, (double)screenSizeY/2 -0.5, 0);
//	Mat4D windowFull = matrixMultiply(translationScreen, windowScale);
//	Mat4D windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
    mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
	
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
	bool autoRotate = true;
	bool showDepth = false;
	double delayTime = 1.0/60;
    double minuteResetTracker = 0;
	while (keepRunning == true) {
		debugLine = 0;
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
        double dTime = float_ms.count()/1000.0;
        mCursesGfxPhysics.update(dTime);
        
        minuteResetTracker += dTime;
        if(minuteResetTracker > 30) {
            minuteResetTracker = 0;
            for(int c = 0; c < numCubes; c++) {
                
                btTransform trans;
                float mat[16];
                for(int i = 0; i < 4; i++) {
                    for(int j = 0; j < 4; j++) {
                        mat[i + 4*j] = cubeInitial[c].d[i][j];
                    }
                }
                trans.setFromOpenGLMatrix(mat);
                cubePhysics[c]->clearForces();
                //                cubePhysics[c]->setGravity(btVector3(0,0,-9.81));
                cubePhysics[c]->setAngularVelocity(btVector3(0,0,0));
                cubePhysics[c]->setLinearVelocity(btVector3(0,0,0));
                cubePhysics[c]->setWorldTransform(trans);
                //                cubePhysics[c]->set
                cubePhysics[c]->getMotionState()->setWorldTransform(trans);
            }
            
            for(int c = 0; c < numCubes; c++) {
                cubePhysics[c]->forceActivationState(DISABLE_DEACTIVATION);
            }
        }
        
		delayTime += 0.001*(1.0/60.0 - dTime);
        if (delayTime> 0 && delayTime < 1.0/60) {
            usleep(1000000.0*delayTime);
        }
//		depthBuffer.reset();
		mRenderPipeline.reset();
		//erase();
//		mvprintw(debugLine++, 0, "tilt: %f", tilt*180.0/M_PI);
		
		
		
		/*
		 Build the view/camera matrix:
		 */
		cameraAxis.x = 0;
		cameraAxis.y = 0;
		cameraAxis.z = 1;
		cameraOrientation = rotationFromAngleAndUnitAxis(angle, cameraAxis);
		cameraOrientation = transpose(cameraOrientation);
		
		double distance = 2*sin(angle/2)+10;
		cameraTranslation = translationMatrix(-(5+distance)*sin(angle), (5+distance) * cos(angle), -(2*cos(angle/5)+distance)-6);

		viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation);
		
		cameraAxis.x = 1;
		cameraAxis.y = 0;
		cameraAxis.z = 0;
		cameraOrientation = rotationFromAngleAndUnitAxis(tilt, cameraAxis);
		cameraOrientation = transpose(cameraOrientation);

		viewMatrix = matrixMultiply( cameraOrientation, viewMatrix);
		
		lightAngle+=0.75 * dTime;
		light[0].x = 6*cos(lightAngle);
		light[0].y = 6*sin(lightAngle);
		light[0].z = 6 + 1.0*sin(lightAngle*5.0);
		light[0].w = 1;
		lightModelView = matrixVectorMultiply(viewMatrix, light[0]);
		
		// Light Model
		LightParams mLightParams;
		mLightParams.numLights = 2;
		mLightParams.modelView[0] = lightModelView;
		mLightParams.color[0].x = 255.0/255.0;//1.25/3.0;
        mLightParams.color[0].y = 142.0/255.0;//0.9/3.0;
		mLightParams.color[0].z = 41.0/255.0;
		
		
		light[1].x = 6*cos(lightAngle*1.25);
		light[1].y = 6*sin(lightAngle*1.25);
		light[1].z = 6 + 1.0*sin(lightAngle*6.0);
		light[1].w = 1;
		mLightParams.modelView[1] = matrixVectorMultiply(viewMatrix, light[1]);
		mLightParams.color[1].x = 200.0/255.0;
		mLightParams.color[1].y = 180.0/255.0;//0.1;
		mLightParams.color[1].z = 255.0/255.0;
		
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

			mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, mRenderPipeline.viewport, (void*)&mLightParams.color[i], debugLine);
		}

		
		
		// Cube 2
		cube2angle += 0.02;
		Coordinates3D cube2roationAxis = {1+sin(0.9*cube2angle), sin(0.8*cube2angle), sin(0.7*cube2angle)};
		cube2roationAxis = normalizeVector(cube2roationAxis);
		Mat4D cube2rotation = rotationFromAngleAndUnitAxis(cube2angle*1.1, cube2roationAxis);
		Mat4D cube2translation = translationMatrix(1, -1, 0.1);
		Mat4D modelViewCube2 = matrixMultiply(cube2translation, cube2rotation);
		modelViewCube2 = matrixMultiply(viewMatrix, modelViewCube2);
        

		// triangle
//		Coordinates4D cubeTriangle[sizeof(triangle)/sizeof(triangle[0])];
		Mat4D solidCubeScale = scaleMatrix(1.05, 1.05, 1.05);
		Coordinates3D solidAxis = {0,0,1};
		Mat4D solidRotation = rotationFromAngleAndUnitAxis(M_PI_4, solidAxis);
		
//        viewMatrix = translationMatrix(0, 0, -2 + sin(angle*4));
        
		Mat4D modelViewBlackCube = matrixMultiply(solidRotation, solidCubeScale);
        btTransform trans;
//        mLightParams.colorModel = {0.25, 0, 0.15};
        for(int i = 0; i < numCubes; i++) {
            Mat4D solidCubeScale = scaleMatrix(cubeScale[i].x, cubeScale[i].y, cubeScale[i].z);
            
            if(cubePhysics[i] && cubePhysics[i]->getMotionState()) {
                cubePhysics[i]->getMotionState()->getWorldTransform(trans);
                //            modelGround = scaleMatrix(50, 50, 50);
                //            Mat4D modelTranslation = translationMatrix(trans.getOrigin().getX(), trans.getOrigin().getY(), trans.getOrigin().getZ());
                //            Mat4D modelSolidCube = matrixMultiply(modelTranslation, solidCubeScale);
                //
                //            Mat4D modelRotation;
                //            for(int i = 0; i < 3; i++) {
                //                modelRotation.d[i][0] = trans.getBasis().getRow(i).getX();
                //                modelRotation.d[i][1] = trans.getBasis().getRow(i).getY();
                //                modelRotation.d[i][2] = trans.getBasis().getRow(i).getZ();
                //            }
                //            modelRotation.d[3][3] = 1;
                //            modelViewBlackCube = matrixMultiply(modelRotation, modelSolidCube);
                
                float mat[16];
                trans.getOpenGLMatrix(mat);
                for(int i = 0; i < 4; i++) {
                    for(int j = 0; j < 4; j++) {
                        modelViewBlackCube.d[i][j] = mat[i + 4*j];
                    }
                }
            } else {
                
            }
            modelViewBlackCube = matrixMultiply( modelViewBlackCube, solidCubeScale);
            modelViewBlackCube = matrixMultiply(viewMatrix, modelViewBlackCube);
            numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
            mUniformInfo.modelView = modelViewBlackCube;
            mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
            mUniformInfo.normalMatrix = invert3x3Slow(mUniformInfo.modelView);
            mUniformInfo.normalMatrix = transpose(mUniformInfo.normalMatrix);
            mRenderPipeline.setFragmentShader(lightFs3);
            mRenderPipeline.mRasterizerThreadPool.busyWait();
            mLightParams.colorModel = &cubeColor[i];
            mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParams, myVertexShader);
        }
        
//        if(cubePhysics2 && cubePhysics2->getMotionState()) {
//            cubePhysics2->getMotionState()->getWorldTransform(trans);
//            float mat[16];
//            trans.getOpenGLMatrix(mat);
//            for(int i = 0; i < 4; i++) {
//                for(int j = 0; j < 4; j++) {
//                    modelViewBlackCube.d[i][j] = mat[i + 4*j];
//                }
//            }
//        } else {
//
//        }
//        modelViewBlackCube = matrixMultiply(viewMatrix, modelViewBlackCube);
//        mUniformInfo.modelView = modelViewBlackCube;
//        mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
//        mRenderPipeline.setFragmentShader(lightFs3);
//        mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParams, myVertexShader);
        
        
        
        
        
        
        for(int i = 0; i < numStatic; i++) {
            //        mRenderPipeline.trianglesFill(squareVi, squareViIndices, 2);
            Mat4D modelTranslation = staticMatrix[i];
            Mat4D modelScale = scaleMatrix(staticScale[i].x, staticScale[i].y, staticScale[i].z);
            Mat4D modelGround = matrixMultiply(modelTranslation, modelScale);
            
            //        mCursesGfxPhysics.groundShape->
            //        mCursesGfxPhysics.ground->;
            
            if(staticCollision[i] && staticCollision[i]->getMotionState()) {
                staticCollision[i]->getMotionState()->getWorldTransform(trans);
                //            trans.deSerializeDouble(**modelGround.d);
                
                //            modelGround = scaleMatrix(8, 8, 8);
                //            modelTranslation = translationMatrix(trans.getOrigin().getX(), trans.getOrigin().getY(), trans.getOrigin().getZ());
                //            modelGround = matrixMultiply(modelTranslation, modelGround);
                //
                //            Mat4D modelRotation;
                //            for(int i = 0; i < 3; i++) {
                //                modelRotation.d[i][0] = trans.getBasis().getRow(i).getX();
                //                modelRotation.d[i][1] = trans.getBasis().getRow(i).getY();
                //                modelRotation.d[i][2] = trans.getBasis().getRow(i).getZ();
                //            }
                //            modelRotation.d[3][3] = 1;
                //
                //            modelGround = matrixMultiply(modelRotation, modelGround);
                
                float mat[16];
                trans.getOpenGLMatrix(mat);
                for(int i = 0; i < 4; i++) {
                    for(int j = 0; j < 4; j++) {
                        modelGround.d[i][j] = mat[i + 4*j];
                    }
                }
                
                modelGround = matrixMultiply(modelGround, modelScale);
            } else {
                //            modelGround;
            }
            Mat4D modelViewSquare = matrixMultiply(viewMatrix, modelGround);
            mUniformInfo.modelView = modelViewSquare;
            mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
            mUniformInfo.normalMatrix = invert3x3Slow(mUniformInfo.modelView);
            mUniformInfo.normalMatrix = transpose(mUniformInfo.normalMatrix);
            //        mRenderPipeline.setFragmentShader(lightFs4);
            //        mRenderPipeline.rasterizeShader(squareVi, &mUniformInfo, squareViIndices, 2, (void*)&mLightParams, myVertexShader);
            
            mRenderPipeline.setFragmentShader(lightFs3);
//            Coordinates3D groundColor = {201.0/255.0, 226.0/255.0, 255.0/255.0};
//            groundColor = groundColor * 0.1;
            mRenderPipeline.mRasterizerThreadPool.busyWait();
            mLightParams.colorModel = &staticColor[i];
            mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParams, myVertexShader);
        }

		if (showDepth) {
			mRenderPipeline.depthBufferToTerminal();
		} else {
			mRenderPipeline.renderBufferToTerminal();
		}
		
		if (autoRotate) {
			angle -= dTime*0.4;
		}

        // HUD
        mvprintw(debugLine++, 0, "FPS: %f",1.0/dTime);
        mvprintw(debugLine++, 0, "Delay time %f", delayTime);
        
		int ch;
//		refresh();
//		continue;
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == KEY_RESIZE) {

			getmaxyx(stdscr, screenSizeY, screenSizeX);

			mRenderPipeline.resize(screenSizeX, screenSizeY);
			
			screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
//			windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
            mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
			
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
	
    mCursesGfxTerminal.cleanupTerminal();
    printf("\n\r");
    
    printf("Window (mRenderPipeline.viewport):\n\r");
    for(int i = 0; i < 4; i++) {
        printf("% 0.2f % 0.2f % 0.2f % 0.2f \n\r", mRenderPipeline.viewport.d[i][0], mRenderPipeline.viewport.d[i][1], mRenderPipeline.viewport.d[i][2], mRenderPipeline.viewport.d[i][3]);
    }
    
    printf("Perspective:\n\r");
    for(int i = 0; i < 4; i++) {
        printf("% 0.2f % 0.2f % 0.2f % 0.2f \n\r", projection.d[i][0], projection.d[i][1], projection.d[i][2], projection.d[i][3]);
    }
    
    printf("View:\n\r");
    for(int i = 0; i < 4; i++) {
        printf("% 0.2f % 0.2f % 0.2f % 0.2f \n\r", viewMatrix.d[i][0], viewMatrix.d[i][1], viewMatrix.d[i][2], viewMatrix.d[i][3]);
    }
    
    
    printf("viewProjection:\n\r");
    Mat4D viewProjection = matrixMultiply(projection, viewMatrix);
    
    for(int i = 0; i < 4; i++) {
        printf("% 0.2f % 0.2f % 0.2f % 0.2f \n\r", viewProjection.d[i][0], viewProjection.d[i][1], viewProjection.d[i][2], viewProjection.d[i][3]);
    }
    
    printf("viewProjectionWindow:\n\r");
    Mat4D viewProjectionWindow = matrixMultiply(mRenderPipeline.viewport, viewProjection);
    
    for(int i = 0; i < 4; i++) {
        printf("% 0.2f % 0.2f % 0.2f % 0.2f \n\r", viewProjectionWindow.d[i][0], viewProjectionWindow.d[i][1], viewProjectionWindow.d[i][2], viewProjectionWindow.d[i][3]);
    }
    
    
    printf("vertexViewProjection:\n\r");
    Coordinates4D vertexViewProjection =  matrixVectorMultiply(viewProjection, squareVi[0].vertex);
    printf("% 0.2f \n\r% 0.2f \n\r% 0.2f \n\r% 0.2f \n\r", vertexViewProjection.x, vertexViewProjection.y, vertexViewProjection.z, vertexViewProjection.w);
    
    printf("vertexViewProjectionWindow:\n\r");
    vertexViewProjection =  matrixVectorMultiply(viewProjectionWindow, squareVi[0].vertex);
    printf("% 0.2f \n\r% 0.2f \n\r% 0.2f \n\r% 0.2f \n\r", vertexViewProjection.x, vertexViewProjection.y, vertexViewProjection.z, vertexViewProjection.w);
	
	return 0;
};


