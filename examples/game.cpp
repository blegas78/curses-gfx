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
#include "curses-gfx-physics.h"
#include "curses-gfx-loader.h"
#include "curses-gfx-resources.h"

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
    Coordinates3D* colorModel;
} LightParams;

typedef struct _LightParamsAndTexture {
    LightParams* lightParams;
    Texture* texture;
    Coordinates3D cameraLocation;
    Coordinates3D colorDiffuse;
} LightParamsAndTexture;


typedef struct _UniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
    Mat4D normalMatrix;
} UniformInfo;

REGISTER_VERTEX_LAYOUT(MeshVertexInfo)
    MEMBER(location),
    MEMBER(normal),
    MEMBER(textureCoord),
    MEMBER(color)
END_VERTEX_LAYOUT(MeshVertexInfo)


template <class T, class U> void myVertexShader(U* uniformInfo, T& output, const T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
    output.normal = matrixVectorMultiply(uniformInfo->normalMatrix, input.normal);
    output.color = input.color;
    output.textureCoord = input.textureCoord;
}



void textureFs(const FragmentInfo2& fInfo) {
    LightParamsAndTexture* lpt = (LightParamsAndTexture*)fInfo.data;
    MeshVertexInfo* vertexInfo = (MeshVertexInfo*)fInfo.interpolated;
    //setRGB(fInfo.pixel, *colorRGB);
    
    
    Coordinates3D normal = normalizeVectorFast(vertexInfo->normal);
//    // So we could take the dot product:
//    const Coordinates3D viewDir = {0,0,1};
//    double magnitude = dotProduct(normal, viewDir);
    // OR just simply use the z component of the normal given our coordinate space:
    double magnitude = normal.z;
//    magnitude *= magnitude;
    magnitude = fmax(0,fmin(1,magnitude));
//    Coordinates3D clippedRGB = clipRGB(*colorRGB);
    *fInfo.colorOutput = lpt->texture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y) * magnitude;
    
//    *fInfo.colorOutput = userTexture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y);
//    fInfo.colorOutput->r = magnitude*255;
//    fInfo.colorOutput->g = magnitude*255;
//    fInfo.colorOutput->b = magnitude*255;
    fInfo.colorOutput->a = 0;
}

void materialFs(const FragmentInfo2& fInfo) {
    LightParamsAndTexture* lpt = (LightParamsAndTexture*)fInfo.data;
    MeshVertexInfo* vertexInfo = (MeshVertexInfo*)fInfo.interpolated;
    
    Coordinates3D normal = normalizeVectorFast(vertexInfo->normal);
    
    double magnitude = normal.z;
    magnitude *= magnitude;
    magnitude = fmin(1,magnitude);
    fInfo.colorOutput->r = magnitude*255*lpt->colorDiffuse.x;
    fInfo.colorOutput->g = magnitude*255*lpt->colorDiffuse.y;
    fInfo.colorOutput->b = magnitude*255*lpt->colorDiffuse.z;
    fInfo.colorOutput->a = 0;
}

void vertexColorFs(const FragmentInfo2& fInfo) {
//    LightParamsAndTexture* lpt = (LightParamsAndTexture*)fInfo.data;
    MeshVertexInfo* vertexInfo = (MeshVertexInfo*)fInfo.interpolated;
    
    Coordinates3D normal = normalizeVectorFast(vertexInfo->normal);
    
    double magnitude = normal.z;
    magnitude *= magnitude;
    magnitude = fmin(1,magnitude);
    fInfo.colorOutput->r = magnitude*vertexInfo->color.r;
    fInfo.colorOutput->g = magnitude*vertexInfo->color.g;
    fInfo.colorOutput->b = magnitude*vertexInfo->color.b;
    fInfo.colorOutput->a = 0;
}

class WorldObject {
    btTransform worldTrans;
    
public:
    bool hideScene;
    Scene scene;
    Scene sceneCollision;
    Mat4D sceneMatrix;
    btRigidBody* body;
    Coordinates3D scale;
    Mat4D worldTransform;
    LightParamsAndTexture* params;
    LightParamsAndTexture* paramsCollision;
    
    WorldObject() {
        hideScene = false;
        sceneMatrix = scaleMatrix(1, 1, 1);
    }
    void load(const char *file) {
        scene.load(file);
        params = new LightParamsAndTexture[scene.numMeshes];
    }
    
    void loadCollisionScene(const char *file) {
        sceneCollision.load(file);
        paramsCollision = new LightParamsAndTexture[scene.numMeshes];
    }
    
    Mat4D& getWorldTransform() {
        float mat[16];
        body->getMotionState()->getWorldTransform(worldTrans);
        worldTrans.getOpenGLMatrix(mat);
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                worldTransform.d[i][j] = mat[i + 4*j];
            }
        }
        return worldTransform;
    }
    
};



int main(int argc, char** argv) {
    
    putenv("SDL_VIDEODRIVER=dummy");
    SDL_Init(SDL_INIT_VIDEO|SDL_INIT_JOYSTICK);
    SDL_Joystick *joystick;
    joystick = SDL_JoystickOpen(0);
    printf("Name: %s\n", SDL_JoystickNameForIndex(0));
    SDL_Event event;
    
    CursesGfxPhysics mCursesGfxPhysics;
    printf("mCursesGfxPhysics.initialize()...\n");
    mCursesGfxPhysics.initialize();
    
    Mat4D identity4 = translationMatrix(0, 0, 0);
    std::vector<WorldObject*> worldObjects;
    
   
    double characterHeading = 0;
    
    // Load a character:
    WorldObject character;
    character.scale = {1,1,1};
//    character.load((std::string(CURSES_GFX_RESOURCE_PATH) + "capsule.dae").c_str());
    character.load((std::string(CURSES_GFX_RESOURCE_PATH) + "/mario/source/Mario/Mario.obj").c_str());
    character.loadCollisionScene((std::string(CURSES_GFX_RESOURCE_PATH) + "capsule.dae").c_str());
    character.sceneMatrix = rotationFromAngleAndUnitAxis(M_PI_2, {1,0,0});
    character.sceneMatrix = matrixMultiply(scaleMatrix(0.0375, 0.0375, 0.0375), character.sceneMatrix);
    character.sceneMatrix = matrixMultiply(translationMatrix(0, 0, -2), character.sceneMatrix);
    character.body = mCursesGfxPhysics.addUprightCapsule(character.scale, 1, translationMatrix(0, 0, 2));
    character.body->forceActivationState(DISABLE_DEACTIVATION);
    character.body->setAngularFactor(0);
    worldObjects.push_back(&character);
    Coordinates2Df characterXYvel = {0,0};
    
    // Camera character
//    btRigidBody camera;
//    mCursesGfxPhysics.dynamicsWorld->addRigidBody(&camera);
//    btHingeConstraint(camera, *character.body, <#const btVector3 &pivotInA#>, <#const btVector3 &pivotInB#>, <#const btVector3 &axisInA#>, <#const btVector3 &axisInB#>)
//    camera.setCollisionFlags(camera.getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
 
    // Ground:
    WorldObject ground;
    ground.scale = {5,5,1};
    ground.load( (std::string(CURSES_GFX_RESOURCE_PATH) + "ground.dae").c_str());
    ground.loadCollisionScene((std::string(CURSES_GFX_RESOURCE_PATH) + "color-cube.dae").c_str());
    ground.sceneMatrix = translationMatrix(0, 0, 1);
    ground.body = mCursesGfxPhysics.addCube(ground.scale, 0, translationMatrix(0, 0, -1));
    worldObjects.push_back(&ground);
    
    
       // Terrain:
       WorldObject terrain;
       terrain.scale = {1,1,1};
       terrain.load( (std::string(CURSES_GFX_RESOURCE_PATH) + "super-mario-bros-level-1-1/source/supermariobros11/supermariobros11.fbx").c_str());
       terrain.loadCollisionScene((std::string(CURSES_GFX_RESOURCE_PATH) + "color-cube.dae").c_str());
       terrain.sceneMatrix = translationMatrix(0, 0, 0);
//       terrain.body = mCursesGfxPhysics.addCube(terrain.scale, 0, translationMatrix(0, 0, -1));
    for(int i = 0; i < terrain.scene.numMeshes; i++ ) {
        terrain.body = mCursesGfxPhysics.addStaticMesh(terrain.scene.meshes[i].vi, terrain.scene.meshes[i].numTriangles);
        
    }
    worldObjects.push_back(&terrain);
    
    
    WorldObject cube;
    cube.scale = {1,1,1};
    cube.load((std::string(CURSES_GFX_RESOURCE_PATH) + "unit-cube.obj").c_str());
    cube.loadCollisionScene((std::string(CURSES_GFX_RESOURCE_PATH) + "unit-cube.obj").c_str());
    cube.body = mCursesGfxPhysics.addCube(cube.scale, 1, translationMatrix(2, 2, 2));
    worldObjects.push_back(&cube);
    
    WorldObject cube2;
    cube2.scale = {1,1,1};
    cube2.load((std::string(CURSES_GFX_RESOURCE_PATH) + "unit-cube.obj").c_str());
    cube2.loadCollisionScene((std::string(CURSES_GFX_RESOURCE_PATH) + "unit-cube.obj").c_str());
    cube2.body = mCursesGfxPhysics.addCube(cube2.scale, .01, translationMatrix(0, 2, 2));
    cube2.body->setCollisionFlags(cube2.body->getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
//    cube2.body->setCollisionFlags(cube2.body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);
    cube2.body->forceActivationState(DISABLE_DEACTIVATION);
//    cube2.body->setCol
    cube2.hideScene = true;
    
    worldObjects.push_back(&cube2);
    
    btVector3 pivotInA = {0, 0, 0};
    btVector3 pivotInB = {0, 5, 0};
    btVector3 axisInA = {0, 0, 1};
    btVector3 axisInB = {0, 0, 1};
//    btcons
    
//    btHingeConstraint testHinge = btHingeConstraint(*cube.body, *cube2.body, pivotInA, pivotInB, axisInA, axisInB);
    btHingeConstraint testHinge = btHingeConstraint(*character.body, *cube2.body, pivotInA, pivotInB, axisInA, axisInB);
    mCursesGfxPhysics.dynamicsWorld->addConstraint(&testHinge);

    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();
    
    UniformInfo mUniformInfo;

	
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
	double lightAngle = 0;
	bool usePerspective = true;
	bool showDepth = false;
    bool showCollision = false;
	double delayTime = 1.0/60;
    double camAngularVelocity = 0;
    double minuteResetTracker = 0;
	while (keepRunning == true) {
		debugLine = 0;
		
		now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> float_ms = (now - before);
		before = now;
		
        double dTime = float_ms.count()/1000.0;
        if(dTime < 0) dTime = 0;
        if(dTime < 1.0/60.0) {
            usleep((1.0/60.0 - dTime)*1000000.0);
            mCursesGfxPhysics.update(1.0/60.0);
        } else {
            mCursesGfxPhysics.update(dTime);
        }
            
        
        mRenderPipeline.reset();
        
        // Get camera information
        cameraOrientation = rotationFromAngleAndUnitAxis(cameraTilt, {1,0,0});
        
        viewMatrix = cube2.getWorldTransform();
        
//                for(int i = 3; i < 4; i++)
//                    for(int j = 0; j < 3; j++) {
//                        viewMatrix.d[j][i] *= -1;
//                    }
        viewMatrix = invert3x3Slow(viewMatrix);
        character.getWorldTransform();
//        viewMatrix = matrixMultiply( viewMatrix, translationMatrix(-cube2.worldTransform.d[0][3], -cube2.worldTransform.d[1][3], -cube2.worldTransform.d[2][3]));
        viewMatrix = matrixMultiply( viewMatrix, translationMatrix(-character.worldTransform.d[0][3], -character.worldTransform.d[1][3], -character.worldTransform.d[2][3]-1));
//        viewMatrix.d[2][3] = -2;
//        viewMatrix = translationMatrix(0, 0, -2);
        cameraTranslation = translationMatrix(0, 0, -7);
//        Mat4D temp = ;
        viewMatrix = matrixMultiply(matrixMultiply(cameraTranslation, cameraOrientation ), viewMatrix);
//        viewMatrix = matrixMultiply(cameraTranslation, cameraOrientation);
        
//        c
//        for(int i = 3; i < 4; i++)
//            for(int j = 0; j < 3; j++) {
//                viewMatrix.d[j][i] *= -1;
//            }
        
        // character matix modify (HACK)
//        character.sceneMatrix = rotationFromAngleAndUnitAxis(M_PI_2, {1,0,0});
//        character.sceneMatrix = matrixMultiply(scaleMatrix(0.075, 0.075, 0.075), character.sceneMatrix);
//        character.sceneMatrix = matrixMultiply(rotationFromAngleAndUnitAxis(characterHeading+M_PI_2, {0,0,1}), character.sceneMatrix);
        character.sceneMatrix = rotationFromAngleAndUnitAxis(M_PI_2, {1,0,0});
        character.sceneMatrix = matrixMultiply(scaleMatrix(0.067, 0.067, 0.067), character.sceneMatrix);
        character.sceneMatrix = matrixMultiply(translationMatrix(0, 0, -1.5), character.sceneMatrix);
        character.sceneMatrix = matrixMultiply(rotationFromAngleAndUnitAxis(characterHeading+M_PI_2, {0,0,1}), character.sceneMatrix);
        
        
        // Render world
        for(std::vector<WorldObject*>::iterator it = worldObjects.begin(); it != worldObjects.end(); it++) {
            WorldObject* worldObject = *it;
            LightParamsAndTexture* params;
            Scene* scene;
            if(!showCollision) {
                if(worldObject->hideScene) {
                    continue;
                }
                scene = &worldObject->scene;
                params = worldObject->params;
            } else {
                scene = &worldObject->sceneCollision;
                params = worldObject->paramsCollision;
            }
            
            for(int m = 0; m < scene->numMeshes; m++) {
                params[m].cameraLocation = {0,0,0};
                
                Coordinates3D modelColor = {1.0,1.0,1.0};

                Mat4D model;
                if(!showCollision) {
//                    model = scaleMatrix(worldObject->scale.x, worldObject->scale.y, worldObject->scale.z);
//                    model = matrixMultiply( worldObject->sceneMatrix, model);
                    model = worldObject->sceneMatrix;
                    
                } else {
                    model = scaleMatrix(worldObject->scale.x, worldObject->scale.y, worldObject->scale.z);
                }
                model = matrixMultiply(worldObject->getWorldTransform(), model);
                
                UniformInfo mUniformData;
                mUniformData.modelView = matrixMultiply(viewMatrix, model);
                mUniformData.modelViewProjection = matrixMultiply(projection, mUniformData.modelView);
                mUniformData.normalMatrix = invert3x3Slow(mUniformData.modelView);
                mUniformData.normalMatrix = transpose(mUniformData.normalMatrix);
                
                
//                params[m].colorDiffuse = scene->materials[scene->meshes[m].materialId].colorDiffuse;
//                mRenderPipeline.setFragmentShader(materialFs);
                
                //                mRenderPipeline.setFragmentShader(vertexColorFs);
                //                mRenderPipeline.setFragmentShader(vertexColorFs);
                
//                mRenderPipeline.rasterizeShader(scene->meshes[m].vi, &mUniformData, scene->meshes[m].numTriangles, &params[m], myVertexShader);
                
                if(scene->materials[scene->meshes[m].materialId].numTextures > 0) {
                    mRenderPipeline.setFragmentShader(textureFs);
                    params[m].texture = &scene->materials[scene->meshes[m].materialId].textures[0];
                    //            mRenderPipeline.setFragmentShader(textureFs);
                    mRenderPipeline.rasterizeShader(scene->meshes[m].vi, &mUniformData, scene->meshes[m].numTriangles, &params[m], myVertexShader);
                } else {
                    if(scene->materials[scene->meshes[m].materialId].hasDiffuseColor) {
                        params[m].colorDiffuse = scene->materials[scene->meshes[m].materialId].colorDiffuse;
                        mRenderPipeline.setFragmentShader(materialFs);
                    } else {
                        mRenderPipeline.setFragmentShader(vertexColorFs);
                    }
    //                mRenderPipeline.setFragmentShader(vertexColorFs);
    //                mRenderPipeline.setFragmentShader(vertexColorFs);
                    
                    mRenderPipeline.rasterizeShader(scene->meshes[m].vi, &mUniformData, scene->meshes[m].numTriangles, &params[m], myVertexShader);
                }
                
            }
        }
        
        

		if (showDepth) {
			mRenderPipeline.depthBufferToTerminal();
		} else {
			mRenderPipeline.renderBufferToTerminal();
		}

        // HUD
        mvprintw(debugLine++, 0, "FPS: %f",1.0/dTime);
        mvprintw(debugLine++, 0, "Delay time %f", delayTime);
        
		int ch = getch();
//        character
//		refresh();
//		continue;
        
        
        btScalar btYaw, btPitch, btRoll;
        cube2.body->getOrientation().getEulerZYX(btYaw, btPitch, btRoll);
        
        
        bool newHeading = false;
        btVector3 vel;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    keepRunning = false;
                    break;
                    
                case SDL_JOYAXISMOTION:
                    switch (event.jaxis.axis) {
                        case 0:
                            characterXYvel.x = (double)event.jaxis.value/32768;
                            if(fabs(characterXYvel.x) < 0.1) {    // deadzone
                                characterXYvel.x = 0;
                            } else {
                                //                                characterHeading = atan2(characterXYvel.y, characterXYvel.x) + btYaw;
                                newHeading = true;
                            }
                            break;
                            
                        case 1:
                            characterXYvel.y = -(double)event.jaxis.value/32768;
                            if(fabs(characterXYvel.y) < 0.1) {    // deadzone
                                characterXYvel.y = 0;
                            } else {
                                //                                characterHeading = atan2(characterXYvel.y, characterXYvel.x) + btYaw;
                                newHeading = true;
                            }
                            break;
                            
                        case 2:
                            
                            camAngularVelocity = (double)event.jaxis.value/32768;
                            if(fabs(camAngularVelocity) < 0.1) {
                                camAngularVelocity = 0;
                            }
                            
                            break;
                            
                        case 3:
                            //                    characterXYvel.y = -(double)event.jaxis.value/32768;
                            //                    cameraTilt += dTime *(double)event.jaxis.value/32768;
                            cameraTiltVelocity = (double)event.jaxis.value/32768;
                            if(fabs(cameraTiltVelocity) < 0.1) {    // deadzone
                                cameraTiltVelocity = 0;
                            }
                            
                        default:
                            break;
                    }
                    break;
                    
                case SDL_JOYBUTTONDOWN:
                    switch (event.jbutton.button) {
                        case 0:
                            vel = character.body->getLinearVelocity();
                            vel.setZ(10);
                            character.body->setLinearVelocity(vel);
                            break;
                            
                        default:
                            break;
                    }
                    break;
                    
                case SDL_KEYDOWN:
                    switch (event.key.keysym.sym) {
                        case SDLK_UP:
                            characterXYvel.y = 1;
                            newHeading = true;
                            break;
                            
                        case SDLK_DOWN:
                            characterXYvel.y = -1;
                            newHeading = true;
                            break;
                            
                        default:
                            break;
                    }
                    break;
                    
                case SDL_KEYUP:
                    switch (event.key.keysym.sym) {
                        case SDLK_UP:
                            characterXYvel.y = 0;
                            //                            newHeading = true;
                            break;
                            
                        case SDLK_DOWN:
                            characterXYvel.y = 0;
                            //                            newHeading = true;
                            break;
                            
                        default:
                            break;
                    }
                    break;
                    
                default:
                    break;
            }
        }
        if(newHeading) {
            characterHeading = atan2(characterXYvel.y, characterXYvel.x) + btYaw;
        }
        
        //            if (event.type == SDL_QUIT) {
        //                keepRunning = false;
//            } else if (event.type == SDL_JOYAXISMOTION) {
//                if(event.jaxis.axis == 0) {
//                    characterXYvel.x = (double)event.jaxis.value/32768;
//                    if(fabs(characterXYvel.x) < 0.1) {    // deadzone
//                        characterXYvel.x = 0;
//                    } else {
//                        characterHeading = atan2(characterXYvel.y, characterXYvel.x) + btYaw;
//
//                    }
//                } else if(event.jaxis.axis == 1) {
//                    characterXYvel.y = -(double)event.jaxis.value/32768;
//                    if(fabs(characterXYvel.y) < 0.1) {    // deadzone
//                        characterXYvel.y = 0;
//                    } else {
//                        characterHeading = atan2(characterXYvel.y, characterXYvel.x) + btYaw;
//                    }
//                } else if(event.jaxis.axis == 3) {
////                    characterXYvel.y = -(double)event.jaxis.value/32768;
////                    cameraTilt += dTime *(double)event.jaxis.value/32768;
//                    cameraTiltVelocity = (double)event.jaxis.value/32768;
//                    if(fabs(cameraTiltVelocity) < 0.1) {    // deadzone
//                        cameraTiltVelocity = 0;
//                    }
//                } else if(event.jaxis.axis == 2) {
//                    camAngularVelocity = (double)event.jaxis.value/32768;
//                if(fabs(camAngularVelocity) < 0.1) {
//                    camAngularVelocity = 0;
//                }
//
//
//                }
//            } else if (event.type == SDL_JOYBUTTONDOWN) {
//                if(event.jbutton.button == 0) {
//                    btVector3 vel = character.body->getLinearVelocity();
//                    vel.setZ(10);
//                    character.body->setLinearVelocity(vel);
//                }
//            }
//        }
        
        
            
		if (ch == 0x1B) {	// Escape
			keepRunning = false;
		} else if (ch == KEY_RESIZE) {

			getmaxyx(stdscr, screenSizeY, screenSizeX);

			mRenderPipeline.resize(screenSizeX, screenSizeY);
			
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
            characterXYvel.x = -1;
		} else if ( ch == KEY_RIGHT) {
            characterXYvel.x = 1;
		} else if ( ch == KEY_UP) {
            characterXYvel.y = 1;
		} else if ( ch == KEY_DOWN) {
            characterXYvel.y = -1;
		} else if ( ch == ' ' ) {
            btVector3 vel = character.body->getLinearVelocity();
            vel.setZ(4);
            character.body->setLinearVelocity(vel);
		} else if ( ch == 'd' || ch == 'D' ) {
			showDepth = !showDepth;
		} else if ( ch == 'c' || ch == 'C' ) {
            showCollision = !showCollision;
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
        
        
        cube2.body->setAngularVelocity({0,0,-3*camAngularVelocity});
        cameraTilt += dTime * cameraTiltVelocity;
        
//        character.body->getLinearVelocity();
        
        vel = character.body->getLinearVelocity();
        vel.setX(3*(characterXYvel.x * cos(btYaw) + -characterXYvel.y*sin(btYaw)));
        vel.setY(-3*(-characterXYvel.x * sin(btYaw) + -characterXYvel.y*cos(btYaw)));
        character.body->setLinearVelocity(vel);
//
//        character.body->getOrientation().getEulerZYX(btYaw, btPitch, btRoll);
//
//        character.body->setAngularVelocity({0, 0, 5*(characterHeading-btYaw)});
	}
	
    mCursesGfxTerminal.cleanupTerminal();

    SDL_JoystickClose(joystick);
    SDL_Quit();
 
	return 0;
};


