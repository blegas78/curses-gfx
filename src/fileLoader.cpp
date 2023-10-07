//#include <png.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>

#include <chrono>

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

#include "curses-gfx.h"
#include "curses-clock.h"
#include "curses-gfx-3d.h"
#include "curses-gfx-handler.h"
#include "curses-gfx-texture.h"

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
    nodelay(stdscr, 1);    // Don't wait at the getch() function if the user hasn't hit a key
    keypad(stdscr, 1); // Allow Function key input and arrow key input

    start_color();
    use_default_colors();
//
//    init_pair(1, COLOR_RED, COLOR_BLACK);
//    init_pair(2, COLOR_GREEN, COLOR_BLACK);
//    init_pair(3, COLOR_CYAN, COLOR_BLACK);
//    init_pair(4, COLOR_BLUE, COLOR_WHITE);
//
//    init_pair(5, COLOR_BLACK, COLOR_RED );
//    init_pair(6, COLOR_BLACK, COLOR_GREEN );
//    init_pair(7, COLOR_BLACK, COLOR_CYAN );
//    init_pair(8, COLOR_WHITE, COLOR_BLUE );
    int line = 0;
    mvprintw(line++, 0, "ncurses COLORS: %d\n", COLORS);
    mvprintw(line++, 0, "ncurses COLOR_PAIRS: %d\n", COLOR_PAIRS  );
    mvprintw(line++, 0, "ncurses has_colors: %d\n", has_colors());
    mvprintw(line++, 0, "ncurses can_change_color: %d\n", can_change_color());
    
    short restoreR[COLORS];
    short restoreG[COLORS];
    short restoreB[COLORS];
    for(short i = 0; i < COLORS; i++) {
        color_content(i, &restoreR[i], &restoreG[i], &restoreB[i]);
    }
    
    init_pair(0, 0, -1);    // White
    
    int numSatLevels = 8;
    for(int i = 0; i < COLORS; i++) {
//        int hueIndex = mod((i-1)*numSatLevels,COLORS);
//        int satIndex = ceil((double)(i-1)/((double)COLORS/(double)(numSatLevels)));
        int hueIndex = mod(i*numSatLevels,COLORS);
        int satIndex = ceil((double)i/((double)COLORS/(double)(numSatLevels)));
        double hue = (double)(hueIndex)/(double)(COLORS) *360.0;
        double sat = ((double)satIndex)/((double)numSatLevels);
        Coordinates3D rgb = hsvToRgb({hue, sat, 1});
        init_color(i, (int)(rgb.x/255.0*1000.0), (int)(rgb.y/255.0*1000.0), (int)(rgb.z/255.0*1000.0));
        init_pair(i, i, -1);
    }
    
    for(int i = 0; i < COLORS; i++) {
        
        attron(COLOR_PAIR(i));
//        attron(A_DIM);
        mvaddch(line, i, 'W');
//        attroff(COLOR_PAIR(i));
    }
    curs_set(0);    // no cursor
}

void cleanupConsole() {
    clear();
    
//    printf("\e[1;1H\e[2J");
    endwin();

    std::cout << "Console has been cleaned!" << std::endl;
}

typedef struct _CubeVertexInfo {
    Coordinates4D vertex;
    Coordinates4D location;
    Coordinates3D normal;
    Coordinates2Df textureCoord;
    ColorRGBA color;
//    int materialIndex;  // Do not interpolate this
} CubeVertexInfo;

REGISTER_VERTEX_LAYOUT(CubeVertexInfo)
    MEMBER(location),
    MEMBER(normal),
    MEMBER(textureCoord),
    MEMBER(color)
END_VERTEX_LAYOUT(CubeVertexInfo)


typedef struct _LightParams {
    Coordinates4D modelView[10];
    Coordinates3D color[10];
    int numLights;
} LightParams;

typedef struct _LightParamsAndTexture {
    LightParams* lightParams;
    Texture* texture;
    Coordinates3D cameraLocation;
} LightParamsAndTexture;

void lightModelFs(const FragmentInfo& fInfo) {
    Coordinates3D* colorRGB = (Coordinates3D*)fInfo.data;
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    //setRGB(fInfo.pixel, *colorRGB);
    
    
    Coordinates3D normal = normalizeVectorFast(vertexInfo->normal);
//    // So we could take the dot product:
//    const Coordinates3D viewDir = {0,0,1};
//    double magnitude = dotProduct(normal, viewDir);
    // OR just simply use the z component of the normal given our coordinate space:
    double magnitude = normal.z;
    magnitude *= magnitude;
    Coordinates3D clippedRGB = clipRGB(*colorRGB);
    fInfo.colorOutput->r = clippedRGB.x*255*magnitude ;
    fInfo.colorOutput->g = clippedRGB.y*255*magnitude ;
    fInfo.colorOutput->b = clippedRGB.z*255*magnitude ;
    fInfo.colorOutput->a = 0;
}

void textureFs(const FragmentInfo& fInfo) {
    LightParamsAndTexture* lpt = (LightParamsAndTexture*)fInfo.data;
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;
    //setRGB(fInfo.pixel, *colorRGB);
    
    
    Coordinates3D normal = normalizeVectorFast(vertexInfo->normal);
//    // So we could take the dot product:
//    const Coordinates3D viewDir = {0,0,1};
//    double magnitude = dotProduct(normal, viewDir);
    // OR just simply use the z component of the normal given our coordinate space:
    double magnitude = normal.z;
    magnitude *= magnitude;
    magnitude = fmin(1,magnitude);
//    Coordinates3D clippedRGB = clipRGB(*colorRGB);
    *fInfo.colorOutput = lpt->texture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y) * magnitude;
    
//    *fInfo.colorOutput = userTexture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y);
//    fInfo.colorOutput->r = magnitude*255;
//    fInfo.colorOutput->g = magnitude*255;
//    fInfo.colorOutput->b = magnitude*255;
    fInfo.colorOutput->a = 0;
}


typedef struct _UniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
} UniformInfo;



template <class T, class U> void myVertexShader(U* uniformInfo, T& output, const T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
    output.normal = matrixVectorMultiply(uniformInfo->modelView, input.normal);
    output.textureCoord.x = input.textureCoord.x;
    output.textureCoord.y = input.textureCoord.y;
}





void lightAndTextureShader(const FragmentInfo& fInfo) {
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;

    LightParamsAndTexture* lpt = (LightParamsAndTexture*) fInfo.data;
    LightParams* lights = lpt->lightParams;
    Texture* texture = lpt->texture;
    
    Coordinates3D colorRGB = {0,0,0};
    
    Coordinates4D lightToFrag;
    Coordinates3D lightDir;
    double intensity2;
//    double lightMagnitudeSquared;
    double lightMagnitude;
    Coordinates3D lightReflected;
    double intensitySpecular;
    double intensity;
    
//    *fInfo.colorOutput =
    ColorRGBA textureSample = texture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y);
    Coordinates3D textureSamplef;
    textureSamplef.x = (double)textureSample.r /255.0;
    textureSamplef.y = (double)textureSample.g /255.0;
    textureSamplef.z = (double)textureSample.b /255.0;
    Coordinates3D vNormal = normalizeVectorFast(vertexInfo->normal);
    Coordinates3D viewDir ;
    viewDir.x = vertexInfo->location.x;// - lpt->cameraLocation.x;   // already in view space
    viewDir.y = vertexInfo->location.y;// - lpt->cameraLocation.y;
    viewDir.z = vertexInfo->location.z;// - lpt->cameraLocation.z;
    viewDir = normalizeVectorFast(viewDir);
    
    for (int i = 0; i < lights->numLights; i++) {
//        lightToFrag = vectorSubtract(lights->modelView[i], fInfo.location3D);
        lightToFrag = vectorSubtract(lights->modelView[i], vertexInfo->location);
        if(dotProduct(lightToFrag, vNormal) < 0) {
            continue;
        }
        lightDir = normalizeVectorFast(lightToFrag);
        double diffuse = dotProduct(vNormal, lightDir);
        
//        lightMagnitude = Q_rsqrt( dotProduct(lightToFrag, lightToFrag) );
//        lightMagnitude = dotProduct(lightToFrag, vNormal);
        
//        lightReflected = vectorScale(vNormal, 2*dotProduct(vNormal, lightDir));
        lightReflected = vNormal * 2.0*dotProduct(vNormal, lightDir);
        lightReflected = vectorSubtract(lightDir, lightReflected);
        lightReflected = normalizeVectorFast(lightReflected);
        
        intensitySpecular = dotProduct(lightReflected, viewDir);
//        intensitySpecular *= Q_rsqrt(dotProduct(vertexInfo->location, vertexInfo->location));
        intensitySpecular = pow(intensitySpecular, 100);
        double distance = 1.0/Q_rsqrt(dotProduct(lightToFrag,lightToFrag));
        double attenuation = 1.0 / ( 0.0 + 0.025*distance + 0.5 *distance*distance);
//        attenuation *= attenuation;
//        intensity = (lightMagnitude*lightMagnitude * intensity2 +  intensitySpecular);
        intensity = (diffuse*3*0  +  intensitySpecular*10 + 3*0) * attenuation;
        
        diffuse *= attenuation*12;
        intensitySpecular *= attenuation*8;
        double ambient = 12.0*attenuation;
        colorRGB = colorRGB + (lights->color[i]*textureSamplef)*diffuse;
        colorRGB = colorRGB + lights->color[i]*intensitySpecular;
        colorRGB = colorRGB + (lights->color[i]*textureSamplef)*ambient;
//        colorRGB = colorRGB + textureSamplef*diffuse;
        
//        colorRGB.x += intensity*lights->color[i].x;
//        colorRGB.y += intensity*lights->color[i].y;
//        colorRGB.z += intensity*lights->color[i].z;
    }
    
    
    Coordinates3D clippedRGB = clipRGB(colorRGB);
    fInfo.colorOutput->r = clippedRGB.x*255;
    fInfo.colorOutput->g = clippedRGB.y*255;
    fInfo.colorOutput->b = clippedRGB.z*255;
    fInfo.colorOutput->a = 0;
    
}







int main(int argc, char** argv) {
    
    
    if(argc != 2) {
        printf("provide input filename as an argument: %s <image.png>\n", argv[0]);
        abort();
    }
    
    Assimp::Importer importer;
    std::string objectFilePath = argv[1];
    std::string objectFileDirectory;
    const size_t last_slash_idx = objectFilePath.rfind('/');
    if (std::string::npos != last_slash_idx)
    {
        objectFileDirectory = objectFilePath.substr(0, last_slash_idx);
    }

    const aiScene* scene = importer.ReadFile( objectFilePath,
//        aiProcess_CalcTangentSpace       |
        aiProcess_Triangulate            |
       // aiProcess_JoinIdenticalVertices  |
                                    aiProcess_FlipUVs |
        aiProcess_GenSmoothNormals  |
//                                             aiProcess_MakeLeftHanded |
//                                             aiProcess_OptimizeMeshes |
//                                             aiProcess_FlipWindingOrder |
        aiProcess_SortByPType);
    
    if (nullptr == scene) {
        std::cerr << "Error loading file: " << importer.GetErrorString() << std::endl;
        return -1;
      }
    struct Mesh {
        CubeVertexInfo* vi;
        int numTriangles;
        bool hasTexture;
        int materialId;
    };
    
    printf("File has %d meshes\n", scene->mNumMeshes);
//    CubeVertexInfo* fileVertices = NULL;
    //int numFileEdges = 0;
    int numMeshes = scene->mNumMeshes;
    Mesh meshes[numMeshes];
    for(int m = 0; m < scene->mNumMeshes; m++) {
        
        printf(" - mesh[%d] has %d vertices = %d triangles\n", m, scene->mMeshes[m]->mNumVertices, scene->mMeshes[m]->mNumVertices/3);
        meshes[m].numTriangles = scene->mMeshes[m]->mNumVertices/3;
        meshes[m].vi = new CubeVertexInfo[scene->mMeshes[m]->mNumVertices];
        
        printf(" - mesh[%d] has normals:   %d\n", m, scene->mMeshes[m]->HasNormals());
        printf(" - mesh[%d][0] has texcoords: %d\n", m, scene->mMeshes[m]->HasTextureCoords(0));
        //        printf(" - texcooords[0][0] name: %s\n", scene->mMeshes[0]->mTextureCoordsNames[0][0].C_Str());
        
        printf(" - - Allocated!\n");
        for(int i = 0; i < scene->mMeshes[m]->mNumVertices; i++) {
            
            meshes[m].vi[i].vertex.x = scene->mMeshes[m]->mVertices[i].x;
            meshes[m].vi[i].vertex.y = scene->mMeshes[m]->mVertices[i].y;
            meshes[m].vi[i].vertex.z = scene->mMeshes[m]->mVertices[i].z;
            meshes[m].vi[i].vertex.w = 1;
            
            meshes[m].vi[i].normal.x = scene->mMeshes[m]->mNormals[i].x;
            meshes[m].vi[i].normal.y = scene->mMeshes[m]->mNormals[i].y;
            meshes[m].vi[i].normal.z = scene->mMeshes[m]->mNormals[i].z;
            
            if(scene->mMeshes[m]->HasTextureCoords(0)) {
                meshes[m].vi[i].textureCoord.x = scene->mMeshes[m]->mTextureCoords[0][i].x;
                meshes[m].vi[i].textureCoord.y = scene->mMeshes[m]->mTextureCoords[0][i].y;
            }
//            else if(scene->mMeshes[m]->HasTextureCoords(1)) {
//                meshes[m].vi[i].textureCoord.x = scene->mMeshes[m]->mTextureCoords[1][i].x;
//                meshes[m].vi[i].textureCoord.y = scene->mMeshes[m]->mTextureCoords[1][i].y;
//
//            }
//            meshes[m].vi[i].materialIndex = scene->mMeshes[m]->mMaterialIndex;
//            printf("%f %f\n", meshes[m].vi[i].textureCoord.x, meshes[m].vi[i].textureCoord.y);
        }
        
        
        meshes[m].materialId = scene->mMeshes[m]->mMaterialIndex;
        printf("meshes[%d].materialId = %d\n", m, meshes[m].materialId);
        
    }
    struct Material {
        Texture* texture;
        int numTextures;
    };
    int numMaterials = scene->mNumMaterials;
    Material mMaterials[numMaterials];
//    Texture *mTextureDiffuse = new Texture[scene->mNumMaterials];   // Likely incorrect
    printf("Scene has textures: %d\n", scene->HasTextures());
    for(int i = 0; i < scene->mNumMaterials; i++) {
        aiMaterial* material = scene->mMaterials[i];
//        printf(" - Name : %s\n", scene->mMaterials[i]->GetName().C_Str());
        printf(" - Material : %d num props:%d\n", i, material->mNumProperties);

        for(int p = 0; p < material->mNumProperties; p++) {
//            printf(" - - Type: %d\n", scene->mMaterials[i]->mProperties[p]->mType);
//            printf(" - - Name: %s\n", scene->mMaterials[i]->mProperties[p]->mKey.C_Str());
//            printf(" - - Semantic: %d\n", scene->mMaterials[i]->mProperties[p]->mSemantic);
//            printf(" - - Data Length: %d\n", scene->mMaterials[i]->mProperties[p]->mDataLength);
//            printf(" - - mIndex %d\n", scene->mMaterials[i]->mProperties[p]->mIndex);
        }
        printf(" - Texture count : %d\n", material->GetTextureCount(aiTextureType_DIFFUSE));
//        for(int t = 0; t < scene->mMaterials[i]->GetTextureCount(aiTextureType_DIFFUSE); t++) {
//            scene->mMaterials[i]->
//        }
        mMaterials[i].numTextures = material->GetTextureCount(aiTextureType_DIFFUSE);
        if(mMaterials[i].numTextures > 0) {
            mMaterials[i].texture = new Texture[mMaterials[i].numTextures];
            //        if(material->GetTextureCount(aiTextureType_DIFFUSE)) {
            for(int t = 0; t < mMaterials[i].numTextures; t++) {
                aiString textureName;
                material->Get(AI_MATKEY_TEXTURE(aiTextureType_DIFFUSE, t), textureName);
                std::string textureFullPath = objectFileDirectory + "/" + textureName.C_Str();
                printf(" = Texture name: %s\n", textureFullPath.c_str());
                if(mMaterials[i].texture[t].loadPng(textureFullPath.c_str()))
                {
                    printf(" ---- Error loading texture\n");
                }
                
                printf(" Texture size: (%d,%d)\n", mMaterials[i].texture[t].width, mMaterials[i].texture[t].height);
                
            }
        } else {
            mMaterials[i].texture = NULL;
        }
        
    }
    
    printf("Loaded!\n");
//        return 0;
    
//    for(int t = 0; t < meshes[0].numTriangles; t++) {
//        printf("%f %f\n", meshes[0].vi[t].textureCoord.x, meshes[0].vi[t].textureCoord.y);
//    }
    
    double minVertex = 1e100;
    double maxVertex = 0;
    for(int m = 0; m < numMeshes; m++) {
        for(int i = 0; i < meshes[m].numTriangles*3; i++) {
            minVertex = fmin(fabs(meshes[m].vi[i].vertex.x), fmin(fabs(meshes[m].vi[i].vertex.y), fmin(fabs(meshes[m].vi[i].vertex.z), minVertex)));
            maxVertex = fmax(fabs(meshes[m].vi[i].vertex.x), fmax(fabs(meshes[m].vi[i].vertex.y), fmax(fabs(meshes[m].vi[i].vertex.z), maxVertex)));
        }
    }
//    return 0;
   
    
    
    setupTerminal();
    

    
    double characterAspect = 28.0/12.0; // macOs terminal
//    double characterAspect = 28.0/14.0; // raspbian terminal
    
    int screenSizeX, screenSizeY;
    getmaxyx(stdscr, screenSizeY, screenSizeX);
    
    RasterizerThreadPool::setRenderThreadCount(4);
    RenderPipeline mRenderPipeline;
    mRenderPipeline.resize(screenSizeX, screenSizeY);
    
    double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
    
    // Model
    Mat4D scaleMat = scaleMatrix(1, 1, 1);
    
    // View
    Mat4D cameraTranslation = translationMatrix(0, 0, -5);
    Coordinates3D cameraAxis = {0, 1, 0};
    cameraAxis = normalizeVector(cameraAxis);
    Mat4D cameraOrientation = rotationFromAngleAndUnitAxis(-M_PI_4, cameraAxis);
    Mat4D viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation );
    
    // Projection
    double zFar = 100;
    double zNear = .1;
    double viewAngle = M_PI*0.5;
    double orthoScale = 5*5;
    Mat4D projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
    
    // Viewport
    mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
    
    // Light
    Coordinates4D light[3];// = {3, 3, 3, 1};
    Coordinates4D lightModelView = {3, 3, 3, 1};
    
    auto now = std::chrono::high_resolution_clock::now();
    auto before = now;
    auto before2 = now;
    
    
    int debugLine = 0;
    int numEdges;
    double lightAngle = 0;
    double angle = 0;
    double tilt = M_PI/4;
    bool usePerspective = true;
    bool autoRotate = true;
    bool showDepth = false;
    bool showGround = true;
    double delayTime = 0;//1.0/60;
    while (keepRunning == true) {
        debugLine = 0;
        
        now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> float_ms = (now - before);
        before = now;
        double dTime = float_ms.count()/1000.0;
        
        delayTime += 0.001*(1.0/60.0 - dTime);
        if (delayTime> 0 && delayTime < 0.1) {
            usleep(1000000.0*delayTime);
        }
        mRenderPipeline.reset();
        
        
        /*
         Build the view/camera matrix:
         */
        cameraAxis.x = 0;
        cameraAxis.y = 0;
        cameraAxis.z = 1;
        cameraOrientation = rotationFromAngleAndUnitAxis(angle, cameraAxis);
        cameraOrientation = transpose(cameraOrientation);
        
        double distance = 2*sin(angle/2)+1;
        cameraTranslation = translationMatrix(-(5+distance)*sin(angle), (5+distance) * cos(angle), -(2*cos(angle/5)+distance)-4);

        viewMatrix = matrixMultiply( cameraOrientation, cameraTranslation);
        
        cameraAxis.x = 1;
        cameraAxis.y = 0;
        cameraAxis.z = 0;
        cameraOrientation = rotationFromAngleAndUnitAxis(tilt, cameraAxis);
        cameraOrientation = transpose(cameraOrientation);

        viewMatrix = matrixMultiply( cameraOrientation, viewMatrix);
        
        
        lightAngle+=dTime * 0.2;
        light[0].x = 6*cos(lightAngle);
        light[0].y = 6*sin(lightAngle);
        light[0].z = 3 + 1.0*sin(lightAngle*5.0);
        light[0].w = 1;
        lightModelView = matrixVectorMultiply(viewMatrix, light[0]);
        
        // Light Model
        LightParams mLightParams;
        mLightParams.numLights = 3;
        mLightParams.modelView[0] = lightModelView;
        mLightParams.color[0].x = 0.20;//1.25/3.0;
        mLightParams.color[0].y = 0.20;//0.9/3.0;
        mLightParams.color[0].z = 3.0/3.0;
        
        
        light[1].x = 6*cos(lightAngle*1.25);
        light[1].y = 6*sin(lightAngle*1.25);
        light[1].z = 3 + 1.0*sin(lightAngle*6.0);
        light[1].w = 1;
        mLightParams.modelView[1] = matrixVectorMultiply(viewMatrix, light[1]);
        mLightParams.color[1].x = 1.0;
        mLightParams.color[1].y = 0.20;//0.1;
        mLightParams.color[1].z = 0.20;
        
        light[2].x = 6*cos(lightAngle*1.75);
        light[2].y = 6*sin(lightAngle*1.75);
        light[2].z = 3 + 1.0*sin(lightAngle*7.0);
        light[2].w = 1;
        mLightParams.modelView[2] = matrixVectorMultiply(viewMatrix, light[2]);
        mLightParams.color[2].x = 0.20;
        mLightParams.color[2].y = 1.0;
        mLightParams.color[2].z = 0.20;
        
        for (int i = 0; i < mLightParams.numLights; i++) {
            Mat4D lightScale = scaleMatrix(0.25, 0.25, 0.25);
            Mat4D lightTranslation = translationMatrix(light[i].x, light[i].y, light[i].z);
            Mat4D lightModel = matrixMultiply(lightTranslation, lightScale);
            Mat4D lightCubeModelView = matrixMultiply(viewMatrix, lightModel);
//            numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
            
            lightModelView = matrixVectorMultiply(viewMatrix, light[i]);
//            mRenderPipeline.setFragmentShader(lightModelFs);

//            mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, mRenderPipeline.viewport, (void*)&mLightParams.color[i], debugLine);
        }

        LightParamsAndTexture mLightParamsAndTexture;
        mLightParamsAndTexture.lightParams = &mLightParams;
//        mLightParamsAndTexture.texture = &pngTexture;
        Coordinates4D zeroVector = {0,0,0,1};
        Coordinates4D camLocation = matrixVectorMultiply(cameraTranslation, zeroVector);
        mLightParamsAndTexture.cameraLocation = {0,0,0};
        
        now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> float_ms2 = (now - before);
        double timeToPrepare = (double)float_ms2.count()/1000.0;
        before2 = now;
        
        
        RenderStats mRenderStats = {0,0,0};
//        mRenderPipeline.setFragmentShader(lightModelFs);
        for(int m = 0; m < numMeshes; m++) {
            Coordinates3D modelColor = {1.0,1.0,1.0};
            Mat4D model = scaleMatrix(6.0/maxVertex, 6.0/maxVertex, 6.0/maxVertex);
            UniformInfo mUniformData;
            mUniformData.modelView = matrixMultiply(viewMatrix, model);
            mUniformData.modelViewProjection = matrixMultiply(projection, mUniformData.modelView);
            //        mRenderPipeline.rasterizeShader(fileVertices, &mUniformData, fileIndices, numFileEdges, &modelColor, myVertexShader);
//            RenderStats mRenderStats2 = mRenderPipeline.rasterizeShader(meshes[m].vi, &mUniformData, meshes[m].numTriangles, &modelColor, myVertexShader);
            mLightParamsAndTexture.texture = &mMaterials[meshes[m].materialId].texture[0];
            mRenderPipeline.setFragmentShader(textureFs);
            RenderStats mRenderStats2 = mRenderPipeline.rasterizeShader(meshes[m].vi, &mUniformData, meshes[m].numTriangles, &mLightParamsAndTexture, myVertexShader);
            
            mRenderStats.timeDrawing += mRenderStats2.timeDrawing;
            mRenderStats.timeClipping += mRenderStats2.timeClipping;
            mRenderStats.timeVertexShading += mRenderStats2.timeVertexShading;
        }
  
        
        now = std::chrono::high_resolution_clock::now();
        float_ms2 = (now - before2);
        double timeToRasterize = (double)float_ms2.count()/1000.0;
        before2 = now;
        
        if (showDepth) {
            mRenderPipeline.depthBufferToTerminal();
        } else {
            mRenderPipeline.renderBufferToTerminal();
        }
        
        
        if (autoRotate) {
            angle -= dTime*0.4;
        }
        
        now = std::chrono::high_resolution_clock::now();
        float_ms2 = (now - before2);
        double timeToRender = (double)float_ms2.count()/1000.0;
        before2 = now;
        
        
        // Compute some stats:
        double totalTime = float_ms.count()/1000.0;
        double timeForNcursesFill = totalTime - timeToRasterize - timeToPrepare;
        
        // HUD
        mvprintw(debugLine++, 0, "FPS: %f", 1000.0/float_ms.count());
        mvprintw(debugLine++, 0, "Render Threads: %d", RasterizerThreadPool::numberRenderThreads);
//        mvprintw(debugLine++, 0, "Delay time %f", delayTime);
        mvprintw(debugLine++, 0, "Total Time  %0.6f", totalTime);
        mvprintw(debugLine++, 0, "Prep        %0.6f % 3.2f%%", timeToPrepare, timeToPrepare/totalTime*100);
        mvprintw(debugLine++, 0, "Rasterize   %0.6f % 3.2f%%", timeToRasterize, timeToRasterize/totalTime*100);
        mvprintw(debugLine++, 0, " - Vertex   %0.6f % 3.2f%%", mRenderStats.timeVertexShading, mRenderStats.timeVertexShading/totalTime*100);
        mvprintw(debugLine++, 0, " - Clip     %0.6f % 3.2f%%", mRenderStats.timeClipping, mRenderStats.timeClipping/totalTime*100);
        mvprintw(debugLine++, 0, " - Draw     %0.6f % 3.2f%%", mRenderStats.timeDrawing, mRenderStats.timeDrawing/totalTime*100);
        mvprintw(debugLine++, 0, "Window Fill %0.6f % 3.2f%%", timeForNcursesFill, timeForNcursesFill/totalTime*100);
//        refresh();
        int ch;
//        refresh();
//        continue;
        if ((ch = getch()) == 0x1B) {    // Escape
            keepRunning = false;
        } else if (ch == KEY_RESIZE) {

            getmaxyx(stdscr, screenSizeY, screenSizeX);

            mRenderPipeline.resize(screenSizeX, screenSizeY);
            
            screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
//            windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
            mRenderPipeline.viewport = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
            
            if (usePerspective) {
                projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
            } else {
                projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
            }
//            depthBuffer.setSize(screenSizeX, screenSizeY);
        } else if (ch == 'o' || ch == 'O') {
            usePerspective    = !usePerspective;
            if (usePerspective) {
                projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
            } else {
                projection = projectionMatrixOrtho(orthoScale*screenAspect, orthoScale, zFar, zNear);
            }
        } else if ( ch == KEY_LEFT) {
            angle -= 0.05;
        } else if ( ch == KEY_RIGHT) {
            angle += 0.05;
        } else if ( ch == KEY_UP) {
            tilt += 0.05;
        } else if ( ch == KEY_DOWN) {
            tilt -= 0.05;
        } else if ( ch == 'g' || ch == 'G' ) {
            showGround = !showGround;
        } else if ( ch == ' ' ) {
            autoRotate = !autoRotate;
        } else if ( ch == 'd' || ch == 'D' ) {
            showDepth = !showDepth;
        } else if ( ch >= '1' && ch <= '9' ) {
            RasterizerThreadPool::setRenderThreadCount(ch-'0');
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
    
    cleanupConsole();
    
    if(meshes[0].vi) {
        delete [] meshes[0].vi;
    }
    return 0;
};



