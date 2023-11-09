#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>

#include <chrono>

#include <SDL2/SDL.h>

#include "curses-gfx.h"
//#include "curses-clock.h"
#include "curses-gfx-3d.h"
#include "curses-gfx-handler.h"
//#include "curses-gfx-texture.h"
#include "curses-gfx-resources.h"
#include "curses-gfx-texture.h"
#include "curses-gfx-loader.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}

enum GBUFFER_TYPE {
    GBUFFER_TYPE_POSITION = 0,
    GBUFFER_TYPE_NORMAL = 1,
    GBUFFER_TYPE_COLOR = 2,
    GBUFFER_TYPE_TEST = 3,
    GBUFFER_TYPE_MATERIAL = 4,
    GBUFFER_TYPE_COUNT = 5
};


//typedef struct _LightParams {
//	Coordinates4D modelView[10];
//	Coordinates3D color[10];
//	int numLights;
//    Coordinates3D* colorModel;
//} LightParams;
//
//typedef struct _LightParamsAndTexture {
//    LightParams* lightParams;
//    Texture* texture;
//    Coordinates3D cameraLocation;
//    Coordinates3D colorDiffuse;
//} LightParamsAndTexture;


typedef struct _UniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
    Mat4D normalMatrix;
//    FrameBuffer* fboArray;
} UniformInfo;

typedef struct _VertexInfo {
    Coordinates4D vertex;
    Coordinates4D location;
    Coordinates3D normal;
    Coordinates2Df textureCoord;
    Coordinates3D color;
} VertexInfo;

REGISTER_VERTEX_LAYOUT(VertexInfo)
    MEMBER(location),
    MEMBER(normal),
    MEMBER(textureCoord),
    MEMBER(color)
END_VERTEX_LAYOUT(VertexInfo)

REGISTER_VERTEX_LAYOUT(MeshVertexInfo)
    MEMBER(location),
    MEMBER(normal),
    MEMBER(textureCoord),
    MEMBER(color)
END_VERTEX_LAYOUT(MeshVertexInfo)

typedef struct _PostUniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
    Mat4D normalMatrix;
    
//    FrameBuffer* fboArray;
} PostUniformInfo;

typedef struct _PostVertexInfo {
    Coordinates4D vertex;
    Coordinates2Df textureCoord;
} PostVertexInfo;

REGISTER_VERTEX_LAYOUT(PostVertexInfo)
    MEMBER(textureCoord)
END_VERTEX_LAYOUT(PostVertexInfo)

template <class T, class U> void myVertexShader(U* uniformInfo, T& output, const T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);//*0.5 + Coordinates4D(.5,.5,.5,0);
    output.normal = matrixVectorMultiply(uniformInfo->normalMatrix, input.normal);// *0.5 + Coordinates3D(0.5, 0.5, 0.5);
    output.color = input.color;
    output.textureCoord = input.textureCoord;
}

template <class T, class U> void myPostVertexShader(U* uniformInfo, T& output, const T& input) {
//    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
//    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
//    output.normal = matrixVectorMultiply(uniformInfo->normalMatrix, input.normal);
//    output.color = input.color;
    output.vertex = input.vertex;//matrixVectorMultiply(scaleMatrix(2, 2, 2), input.vertex);

    output.textureCoord = input.textureCoord;
}



void gbufferFragmentShader(const FragmentInfo2& fInfo) {
    VertexInfo *mVertexInfo =  (VertexInfo*)fInfo.interpolated;
    
    FrameBuffer* fboArray = (FrameBuffer*)fInfo.data;

//    *fInfo.colorOutput = {127, 0, 127, 0};
    Coordinates3D normalNormalized = normalizeVectorFast( mVertexInfo->normal );
    //    *fInfo.colorOutput = {normalNormalized.x, normalNormalized.y, normalNormalized.z, 0};
    
    fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = mVertexInfo->location;
//    fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = {0.5, 0.5, 0.1, 0};
    fboArray[GBUFFER_TYPE_NORMAL].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = normalNormalized;
    fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = mVertexInfo->color;
    fboArray[GBUFFER_TYPE_MATERIAL].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y) = {0,0,0,0};
    
////        *fInfo.colorOutput = mVertexInfo->color;
//    fInfo.colorOutput->r = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).x;
//    fInfo.colorOutput->g = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).y;
//    fInfo.colorOutput->b = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).z;
//    fInfo.colorOutput->a = 0;
////    *fInfo.colorOutput = fboArray[GBUFFER_TYPE_COLOR].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y);
}

void gbufferFlatFragmentShader(const FragmentInfo2& fInfo) {
    VertexInfo *mVertexInfo =  (VertexInfo*)fInfo.interpolated;
    
    FrameBuffer* fboArray = (FrameBuffer*)fInfo.data;

//    *fInfo.colorOutput = {127, 0, 127, 0};
    Coordinates3D normalNormalized = normalizeVectorFast( mVertexInfo->normal );
    //    *fInfo.colorOutput = {normalNormalized.x, normalNormalized.y, normalNormalized.z, 0};
    
    fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = mVertexInfo->location;
//    fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = {0.5, 0.1, 0.5, 0};
    fboArray[GBUFFER_TYPE_NORMAL].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = normalNormalized;
    fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = {1,1,1};
    fboArray[GBUFFER_TYPE_MATERIAL].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y) = {0,0,0,0};
//        *fInfo.colorOutput = mVertexInfo->color;
//    *fInfo.colorOutput = {0.5, 0.5, 0.5, 0};
//    *fInfo.colorOutput = fboArray[GBUFFER_TYPE_COLOR].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y);
}

struct LightParams {
    Coordinates3D position;
//    Coordinates3D positionView;
    Mat4D modelView;
    Coordinates3D color;
//    double intensity;
    double constant;
    double linear;
    double quadratic;
    double radius;
    FrameBuffer* fboArray;
};

void gbufferLightSourceFragmentShader(const FragmentInfo2& fInfo) {
    VertexInfo *mVertexInfo =  (VertexInfo*)fInfo.interpolated;
    LightParams* lp = (LightParams*)fInfo.data;
    FrameBuffer* fboArray = lp->fboArray;

//    *fInfo.colorOutput = {127, 0, 127, 0};
    Coordinates3D normalNormalized = normalizeVectorFast( mVertexInfo->normal );
    //    *fInfo.colorOutput = {normalNormalized.x, normalNormalized.y, normalNormalized.z, 0};
    
    fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = mVertexInfo->location;
    fboArray[GBUFFER_TYPE_NORMAL].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = normalNormalized;
    fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y) = lp->color;
    
    fboArray[GBUFFER_TYPE_MATERIAL].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y) = ColorRGBA(Coordinates4D(lp->color,0.5));
    
//        *fInfo.colorOutput = mVertexInfo->color;
//    fInfo.colorOutput->r = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).x;
//    fInfo.colorOutput->g = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).y;
//    fInfo.colorOutput->b = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y).z;
//    fInfo.colorOutput->a = 0;
//    *fInfo.colorOutput = fboArray[GBUFFER_TYPE_COLOR].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y);
}


void testShader(const FragmentInfo2& fInfo) {
    VertexInfo *mVertexInfo =  (VertexInfo*)fInfo.interpolated;
    *fInfo.colorOutput = {127,127,127,0};
}

void lightingShader(const FragmentInfo2& fInfo) {
    VertexInfo *mVertexInfo = (VertexInfo*)fInfo.interpolated;
    LightParams* lp = (LightParams*)fInfo.data;
    FrameBuffer* fboArray = lp->fboArray;
    
    Coordinates3D normal = fboArray[GBUFFER_TYPE_NORMAL].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y);
    Coordinates3D color = fboArray[GBUFFER_TYPE_COLOR].at<Coordinates3D>(fInfo.pixel.x, fInfo.pixel.y);
    Coordinates4D* position4D = &fboArray[GBUFFER_TYPE_POSITION].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y);
    Coordinates3D position = Coordinates3D(position4D->x, position4D->y, position4D->z);
    
    ColorRGBA material = fboArray[GBUFFER_TYPE_MATERIAL].at<ColorRGBA>(fInfo.pixel.x, fInfo.pixel.y);
    if (material.a != 0) {
        fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = Coordinates4D((double)material.r/255.0, (double)material.g/255.0, (double)material.b/255.0, 0);
        *fInfo.colorOutput = fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y);
        return;
    }
//
//    Coordinates3D albedo = color * 0.1;
//
//
//    Coordinates3D viewDir = normalizeVectorFast({-position.x, -position.y, -position.z});
//
//    Coordinates3D posToLight = vectorSubtract(lp->position, position);
//    double invDistPosToLight = Q_rsqrt(dotProduct(posToLight, posToLight));
////    invDistPosToLight *= invDistPosToLight;
//    Coordinates3D lightDir = normalizeVectorFast(posToLight);
//    Coordinates3D diffuse = lp->color * (albedo * fmax(dotProduct(normal, lightDir), 0.0)*invDistPosToLight);
//
//
    
//    fInfo.colorOutput->r += clamp(diffuse.x*lp->intensity*255, 0, 255);
//    fInfo.colorOutput->g += clamp(diffuse.y*lp->intensity*255, 0, 255);
//    fInfo.colorOutput->b += clamp(diffuse.z*lp->intensity*255, 0, 255);
//    fInfo.colorOutput->a = 0;
    
    Coordinates3D ambient = color * 0.1;
    
    Coordinates3D lightPositionView( lp->modelView.d[0][3], lp->modelView.d[1][3], lp->modelView.d[2][3]);
    
    Coordinates3D posToLight = vectorSubtract(lightPositionView, position);
    double distance = 1.0/Q_rsqrt(dotProduct(posToLight, posToLight));
//    double distance = sqrt(dotProduct(posToLight, posToLight));
//    double distance = sqrt(dotProduct(position,position));
    
    Coordinates3D lightDir = normalizeVectorFast(posToLight);
    double diff = fmax(dotProduct(normal, lightDir), 0.0);
    Coordinates3D diffuse = color * diff;// *0.5;
    
    Coordinates3D viewDir = normalizeVectorFast({-position.x, -position.y, -position.z});
//    Coordinates3D viewDir = normalizeVectorFast({position.x, position.y, position.z});
//    Coordinates3D reflectDir = reflect(lightDir * -1, normal);
//    Coordinates3D reflectDir = vectorScale(normal, 2*dotProduct(normal, lightDir*-1));
//    reflectDir = vectorSubtract(lightDir*-1, reflectDir);
    Coordinates3D reflectDir = reflect(lightDir * -1, normal);
    
    double spec = pow(fmax(dotProduct(viewDir, reflectDir), 0.0), 32);
    Coordinates3D specular = Coordinates3D(1,1,1) * (spec);
    
    
    double attenuation = 1.0/(lp->constant + lp->linear*distance + lp->quadratic*distance*distance);
    
    Coordinates3D result = lp->color * (ambient + diffuse + specular) * attenuation ;
        
    fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) += Coordinates4D(result, 0);
//    fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = Coordinates4D(0.5,0.5,0.5,0);
//    fInfo.colorOutput->r = clamp((int)fInfo.colorOutput->r + result.x*255, 0, 255);
//    fInfo.colorOutput->g = clamp((int)fInfo.colorOutput->g + result.y*255, 0, 255);
//    fInfo.colorOutput->b = clamp((int)fInfo.colorOutput->b + result.z*255, 0, 255);
//    fInfo.colorOutput->a = 0;
    
    
    
//    fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y) = {1,0,1,0};
    *fInfo.colorOutput = fboArray[GBUFFER_TYPE_TEST].at<Coordinates4D>(fInfo.pixel.x, fInfo.pixel.y);
}





int main(int argc, char** argv) {
    
    FrameBuffer fboArray[GBUFFER_TYPE_COUNT];
//    FrameBuffer fboArray[3];

    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();
    
    UniformInfo mUniformInfo;
    
    Scene sphereScene;
    sphereScene.load((std::string(CURSES_GFX_RESOURCE_PATH) + "unit-sphere-low-poly.dae").c_str());
    
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
    
    Coordinates3D cubeVertexColor[8];
    for(int i = 0; i < 8; i++) {
        cubeVertexColor[i].x = (cube[i].x+1)/2;
        cubeVertexColor[i].y = (cube[i].y+1)/2;
        cubeVertexColor[i].z = (cube[i].z+1)/2;
//        cubeVertexColor[i].a = 0;
    }
    
    int cubeTriangleIndices[][3] = {
        {0, 1, 2},    // left
        {2, 3, 0},
        {4, 6, 5},    // right
        {6, 4, 7},
        {6, 2, 5},     // top
        {2, 1, 5},
        {0, 3, 4},     // bottom
        {4, 3, 7},
        {2, 6, 3},     // front
        {7, 3, 6},
        {0, 4, 5},     // back
        {5, 1, 0}
    };
    
    VertexInfo cubeVi[12*3];
    int cubeViIndices[12][3];
    for(int t = 0; t < 12; t++) {   // for each triangle in a cube
        for(int v = 0; v < 3; v++) {
            cubeVi[v + 3*t].vertex = cube[cubeTriangleIndices[t][v]];
            cubeViIndices[t][v] = v + 3*t;
            cubeVi[v + 3*t].color = cubeVertexColor[cubeTriangleIndices[t][v]];
//            cubeVi[v + 3*t].color = {1, 1, 1};
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
    
    
    // For post processing, need a basic square to rasterize:
    PostVertexInfo squareVi[4];
    PostVertexInfo squareViOrdered[6];
    int squareViIndices[2][3];
    squareVi[0].vertex = {-1, -1, 1, 1};
    squareVi[1].vertex = { 1, -1, 1, 1};
    squareVi[2].vertex = { 1,  1, 1, 1};
    squareVi[3].vertex = {-1,  1, 1, 1};
    
    squareVi[0].textureCoord = {0, 0};
    squareVi[1].textureCoord = {0, 1};
    squareVi[2].textureCoord = {1, 1};
    squareVi[3].textureCoord = {1, 0};
    
    squareViIndices[0][0] = 0;  // right handed
    squareViIndices[0][1] = 1;
    squareViIndices[0][2] = 2;
    squareViIndices[1][0] = 0;  // left handed
    squareViIndices[1][1] = 2;
    squareViIndices[1][2] = 3;
    
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 3; j++) {
            squareViOrdered[i*3 + j].textureCoord = squareVi[squareViIndices[i][j]].textureCoord;
            squareViOrdered[i*3 + j].vertex = squareVi[squareViIndices[i][j]].vertex;
        }
    }
    
	
	double characterAspect = 28.0/12.0; // macOs terminal
//	double characterAspect = 28.0/14.0; // raspbian terminal
//	double characterAspect = 6.0/4.0; // zipitZ2
	
	int screenSizeX, screenSizeY;
	getmaxyx(stdscr, screenSizeY, screenSizeX);
	
//    RasterizerThreadPool::setRenderThreadCount(4);
	RenderPipeline mRenderPipeline;
//    mRenderPipeline.fbo->type = FBT_COORDINATES4D;
	mRenderPipeline.resize(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_COLOR].setSize<Coordinates3D>(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_COLOR].type = FBT_COORDINATES3D;
    fboArray[GBUFFER_TYPE_NORMAL].setSize<Coordinates3D>(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_NORMAL].type = FBT_COORDINATES3D;
    fboArray[GBUFFER_TYPE_POSITION].setSize<Coordinates4D>(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_POSITION].type = FBT_COORDINATES4D;
    fboArray[GBUFFER_TYPE_TEST].setSize<Coordinates4D>(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_TEST].type = FBT_COORDINATES4D;
    fboArray[GBUFFER_TYPE_MATERIAL].setSize<ColorRGBA>(screenSizeX, screenSizeY);
    fboArray[GBUFFER_TYPE_MATERIAL].type = FBT_RGBA;
	
	double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
	
	// Depth buffer
//	DepthBuffer depthBuffer;
//	depthBuffer.setSize(screenSizeX, screenSizeY);
	
	
	// View
	Mat4D cameraTranslation = translationMatrix(0, 0, -15);
	Coordinates3D cameraAxis = {1, 0, 0};
	cameraAxis = normalizeVector(cameraAxis);
    double cameraTilt = -M_PI_4;
    double cameraTiltVelocity = 0;
	Mat4D cameraOrientation = rotationFromAngleAndUnitAxis(cameraTilt, cameraAxis);
    Mat4D viewMatrix = matrixMultiply(  cameraTranslation, cameraOrientation );
	
	// Projection
	double zFar = 1000;
	double zNear = 0.1;
    
    double viewAngle = M_PI*0.25;
    double orthoScale = 5*5;
	Mat4D projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
	
	// Viewport
    mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
	
	// Light
	Coordinates4D light[3];// = {3, 3, 3, 1};
	Coordinates4D lightModelView = {3, 3, 3, 1};
	
	auto now = std::chrono::high_resolution_clock::now();
	auto before = now;
	
	int debugLine = 0;
	int numEdges;
	double cubeAngle = 0;
    double lightAngle = 134123;
    bool cubeRotateEnable = true;
	bool usePerspective = true;
	bool showDepth = false;
	double delayTime = 1.0/60;
    double camAngularVelocity = 0;
    double minuteResetTracker = 0;
    int bufferToShow = GBUFFER_TYPE_TEST;
	while (keepRunning == true) {
		debugLine = 0;
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
        double dTime = float_ms.count()/1000.0;
        if(dTime < 0) dTime = 0;
            
        
        mRenderPipeline.reset();
        
        fboArray[GBUFFER_TYPE_COLOR].clear(Coordinates3D(0.5,0.5,0.5));
        fboArray[GBUFFER_TYPE_NORMAL].clear(Coordinates3D(0.5,0.5,0.5));
        fboArray[GBUFFER_TYPE_POSITION].clear(Coordinates4D(0,0,-1000,0));
        fboArray[GBUFFER_TYPE_TEST].clear(Coordinates4D(0.,0.,0.,0));
        fboArray[GBUFFER_TYPE_MATERIAL].clear(ColorRGBA(0,0,0,0));
        
        // Get camera information
        cameraOrientation = rotationFromAngleAndUnitAxis(cameraTilt, {1,0,0});
        viewMatrix = matrixMultiply(cameraTranslation, cameraOrientation );
        
        // Light, with rendering:
        int numLights = 5;
        LightParams mLightParams[numLights];
        for(int i = 0; i < numLights; i++) {
            mLightParams[i].color = {.5*(1+sin(i)), .5*(1+cos(i)), .5*(1+sin(i*5.3))};
            mLightParams[i].color = rgbToHsv(mLightParams[i].color);
            mLightParams[i].color.z = 1;
            mLightParams[i].color = hsvToRgb(mLightParams[i].color) * (1.0/255.0);
//            mLightParams[i].color = {1,1,1};
            
            mLightParams[i].constant = 1.0;
            mLightParams[i].linear = 0.22;
            mLightParams[i].quadratic = 0.40;
            float lightMax  = std::fmaxf(std::fmaxf(mLightParams[i].color.x, mLightParams[i].color.y), mLightParams[i].color.z);
            mLightParams[i].radius = (-mLightParams[i].linear +  std::sqrtf(mLightParams[i].linear * mLightParams[i].linear - 4 * mLightParams[i].quadratic * (mLightParams[i].constant - (256.0 / 26.0) * lightMax)))
              / (2 * mLightParams[i].quadratic);

            
            double radius = (double)i/(double)numLights * 10;
            mLightParams[i].position = {radius*sin(lightAngle*2*0.1*(double)(i+1)),radius*cos(lightAngle*2 *0.1*(double)(i+1)),0.5};
//            mLightParams[i].positionView = matrixVectorMultiply(viewMatrix, mLightParams[i].position);
            mLightParams[i].modelView = matrixMultiply(viewMatrix, matrixMultiply(translationMatrix(mLightParams[i].position.x, mLightParams[i].position.y, mLightParams[i].position.z), scaleMatrix(0.2, 0.2, 0.2)));
            
//            mLightParams[i].position = matrixVectorMultiply(viewMatrix, mLightParams[i].position);
//            mLightParams[i].intensity = 4000;
            mLightParams[i].fboArray = fboArray;
            
            mUniformInfo.modelView = mLightParams[i].modelView;
            mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
            mUniformInfo.normalMatrix = invert3x3Slow(mUniformInfo.modelView);
            mUniformInfo.normalMatrix = transpose(mUniformInfo.normalMatrix);
            
            mRenderPipeline.setFragmentShader(gbufferLightSourceFragmentShader);
            mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, 12, &mLightParams[i], myVertexShader);
        }

        
        int numCubes = 1 + 1 + 5*5;
        void (*fragmentShader[numCubes])(const FragmentInfo2&);
//        fragmentShader[0] = gbufferFragmentShader;
        fragmentShader[0] = gbufferFlatFragmentShader;
        fragmentShader[1] = gbufferFragmentShader;
        
        // Cube
        if(cubeRotateEnable)
            cubeAngle += .25*dTime;
        Mat4D cubeModel[numCubes];
        cubeModel[1] = matrixMultiply(matrixMultiply(translationMatrix(0, 0, 6), rotationFromAngleAndUnitAxis(cubeAngle, normalizeVector({sin(cubeAngle/4),cos(cubeAngle/5),2}))), scaleMatrix(3, 3, 3));
        cubeModel[0] = matrixMultiply(translationMatrix(0, 0, -10), scaleMatrix(10, 10, 10));
        
        for(int i = 0; i < 5; i++) {
            for(int j = 0; j < 5; j++) {
                cubeModel[2+i+5*j] = matrixMultiply(rotationFromAngleAndUnitAxis(cubeAngle, {0,0,1}), matrixMultiply(translationMatrix(3*(i-2), 3*(j-2), 1), scaleMatrix(1, 1, 1)));
                fragmentShader[2+i+5*j] = gbufferFlatFragmentShader;
            }
        }
        
        for(int i = 0; i < numCubes; i++) {
            mUniformInfo.modelView = matrixMultiply(viewMatrix, cubeModel[i]);
            mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
            mUniformInfo.normalMatrix = invert3x3Slow(mUniformInfo.modelView);
            mUniformInfo.normalMatrix = transpose(mUniformInfo.normalMatrix);
            //        mUniformInfo.fboArray = fboArray;
            //        mRenderPipeline.backfaceCulling = false;
            mRenderPipeline.setFragmentShader(fragmentShader[i]);
            if(i == 0 ) {
                mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, 12, (void*)fboArray, myVertexShader);
            } else {
                mRenderPipeline.rasterizeShader(sphereScene.meshes[0].vi, &mUniformInfo, sphereScene.meshes[0].numTriangles, (void*)fboArray, myVertexShader);
            }
            
        }
        
//        mRenderPipeline.mRasterizerThreadPool.busyWait();
        
        
        
        mRenderPipeline.depthTestEnable = false;
        mRenderPipeline.backfaceCulling = false;
        mRenderPipeline.setFragmentShader(lightingShader);
//        mRenderPipeline.setFragmentShader(testShader);
        mRenderPipeline.reset();
        
        for(int i = 0; i < numLights; i++) {
//            mRenderPipeline.mRasterizerThreadPool.busyWait();
//            mRenderPipeline.rasterizeShader(squareViOrdered, &mUniformInfo, 2, &mLightParams[i], myPostVertexShader);
            Mat4D model = matrixMultiply(translationMatrix(mLightParams[i].position.x, mLightParams[i].position.y, mLightParams[i].position.z), scaleMatrix(mLightParams[i].radius, mLightParams[i].radius, mLightParams[i].radius));
            mUniformInfo.modelView = matrixMultiply(viewMatrix, model);
            mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
            mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, 12, &mLightParams[i], myVertexShader);
//            mRenderPipeline.rasterizeShader(sphereScene.meshes[0].vi, &mUniformInfo, sphereScene.meshes[0].numTriangles, &mLightParams[i], myVertexShader);
        }
        
//        mRenderPipeline.rasterizeShader(squareViOrdered, &mUniformInfo, 2, &mLightParams[0], myPostVertexShader);
//        mRenderPipeline.rasterizeShader(squareViOrdered, &mUniformInfo, 2, &mLightParams[9], myPostVertexShader);

        mRenderPipeline.mRasterizerThreadPool.busyWait();
        mRenderPipeline.depthTestEnable = true;
        mRenderPipeline.backfaceCulling = true;
        

		if (showDepth) {
			mRenderPipeline.depthBufferToTerminal();
//            mRenderPipeline.renderBufferToTerminal(&fboArray[GBUFFER_TYPE_TEST]);
//            mRenderPipeline.renderBufferToTerminal();
		} else {
			mRenderPipeline.renderBufferToTerminal(&fboArray[bufferToShow]);
		}

        // HUD
        mvprintw(debugLine++, 0, "FPS: %f", 1.0/dTime);
        mvprintw(debugLine++, 0, "Delay time %f", delayTime);
        for(int i = 0; i < numLights; i++)
            mvprintw(debugLine++, 0, "light[%d].radius: %0.2f", i, mLightParams[i].radius);
        
        const char* bufferName = "?";
        if(showDepth) {
            bufferName = "Depth Buffer";
        } else {
            switch(bufferToShow) {
                case GBUFFER_TYPE_TEST: bufferName = "GBUFFER_TYPE_TEST"; break;
                case GBUFFER_TYPE_POSITION: bufferName = "GBUFFER_TYPE_POSTION"; break;
                case GBUFFER_TYPE_NORMAL: bufferName = "GBUFFER_TYPE_NORMAL"; break;
                case GBUFFER_TYPE_COLOR: bufferName = "GBUFFER_TYPE_COLOR"; break;
                case GBUFFER_TYPE_MATERIAL: bufferName = "GBUFFER_TYPE_MATERIAL"; break;
                default: break;
            }
        }
        mvprintw(debugLine++, 0, "FBO: %s", bufferName);
        
		int ch = getch();
		if (ch == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == KEY_RESIZE) {

			getmaxyx(stdscr, screenSizeY, screenSizeX);

			mRenderPipeline.resize(screenSizeX, screenSizeY);
            fboArray[GBUFFER_TYPE_COLOR].setSize<Coordinates3D>(screenSizeX, screenSizeY);
            fboArray[GBUFFER_TYPE_NORMAL].setSize<Coordinates3D>(screenSizeX, screenSizeY);
            fboArray[GBUFFER_TYPE_POSITION].setSize<Coordinates4D>(screenSizeX, screenSizeY);
            fboArray[GBUFFER_TYPE_TEST].setSize<Coordinates4D>(screenSizeX, screenSizeY);
            fboArray[GBUFFER_TYPE_MATERIAL].setSize<ColorRGBA>(screenSizeX, screenSizeY);
			
			screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
//			windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
            mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
			
			if (usePerspective) {
				projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
			} else {
				projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
			}
//			depthBuffer.setSize(screenSizeX, screenSizeY);
		} else if (ch == 'o' || ch == 'O') {
			usePerspective	= !usePerspective;
			if (usePerspective) {
				projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
			} else {
				projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
			}
		} else if ( ch == KEY_LEFT) {
//            characterXYvel.x = -1;
            cubeAngle -= 0.1;
		} else if ( ch == KEY_RIGHT) {
//            characterXYvel.x = 1;
            cubeAngle += 0.1;
		} else if ( ch == KEY_UP) {
//            characterXYvel.y = 1;
            cameraTilt -= 0.05;
		} else if ( ch == KEY_DOWN) {
//            characterXYvel.y = -1;
            cameraTilt += 0.05;
		} else if ( ch == ' ' ) {
            cubeRotateEnable = !cubeRotateEnable;
        } else if ( ch == 'b' || ch == 'B' ) {
            if(++bufferToShow == GBUFFER_TYPE_COUNT)
                bufferToShow = 0;
		} else if ( ch == 'd' || ch == 'D' ) {
			showDepth = !showDepth;
		} else if ( ch == 'c' || ch == 'C' ) {
//            showCollision = !showCollision;
        } else if ( ch == 'r' || ch == 'R') {
            if (usePerspective) {
                viewAngle += M_PI * 0.0125;
                if(viewAngle >= M_PI) {
                    viewAngle = M_PI - 0.0001;
                }
                projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
            } else {
                orthoScale += 0.5;  // this can grow as large as we want
                projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
            }
        } else if ( ch == 't' || ch == 'T') {
            if (usePerspective) {
                viewAngle -= M_PI * 0.0125;
                if(viewAngle < 0.0001) {
                    viewAngle = 0.0001;
                }
                projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
            } else {
                orthoScale -= 0.5;
                if(orthoScale  < 0.0100) {
                    orthoScale = 0.01;
                }
                projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
            }
        }
        
	}
	
    mCursesGfxTerminal.cleanupTerminal();

 
	return 0;
};


