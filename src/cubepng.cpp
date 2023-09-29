#include <png.h>

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
    init_pair(1, COLOR_RED, COLOR_BLACK);
    init_pair(2, COLOR_YELLOW, COLOR_BLACK);
    init_pair(3, COLOR_GREEN, COLOR_BLACK);
    init_pair(4, COLOR_CYAN, COLOR_BLACK);
    init_pair(5, COLOR_BLUE, COLOR_BLACK);
    init_pair(6, COLOR_MAGENTA, COLOR_BLACK);
    init_pair(7, COLOR_WHITE, COLOR_BLACK);


//    init_pair(5, COLOR_BLACK, COLOR_RED );
//    init_pair(6, COLOR_BLACK, COLOR_GREEN );
//    init_pair(7, COLOR_BLACK, COLOR_CYAN );
//    init_pair(8, COLOR_WHITE, COLOR_BLUE );

    curs_set(0);    // no cursor

//    atexit(destroy);
}

void cleanupConsole() {
    clear();
    
//    printf("\e[1;1H\e[2J");
    endwin();

    std::cout << "Console has been cleaned!" << std::endl;
}

typedef struct _PngLoader {
    int width, height;
    png_byte color_type;
    png_byte bit_depth;
    png_bytep *row_pointers = NULL;
} PngLoader;


// thanks to https://gist.github.com/niw/5963798
void read_png_file(char *filename, PngLoader& mPngLoader) {
  FILE *fp = fopen(filename, "rb");

  png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if(!png) abort();

  png_infop info = png_create_info_struct(png);
  if(!info) abort();

  if(setjmp(png_jmpbuf(png))) abort();

  png_init_io(png, fp);

  png_read_info(png, info);

    mPngLoader.width      = png_get_image_width(png, info);
    mPngLoader.height     = png_get_image_height(png, info);
    mPngLoader.color_type = png_get_color_type(png, info);
    mPngLoader.bit_depth  = png_get_bit_depth(png, info);

  // Read any color_type into 8bit depth, RGBA format.
  // See http://www.libpng.org/pub/png/libpng-manual.txt

  if(mPngLoader.bit_depth == 16)
    png_set_strip_16(png);

  if(mPngLoader.color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(png);

  // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
  if(mPngLoader.color_type == PNG_COLOR_TYPE_GRAY && mPngLoader.bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png);

  if(png_get_valid(png, info, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png);

  // These color_type don't have an alpha channel then fill it with 0xff.
  if(mPngLoader.color_type == PNG_COLOR_TYPE_RGB ||
     mPngLoader.color_type == PNG_COLOR_TYPE_GRAY ||
     mPngLoader.color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

  if(mPngLoader.color_type == PNG_COLOR_TYPE_GRAY ||
     mPngLoader.color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png);

  png_read_update_info(png, info);

  if (mPngLoader.row_pointers)
      delete [] mPngLoader.row_pointers;

    mPngLoader.row_pointers = new png_bytep[mPngLoader.height];// (png_bytep*)malloc(sizeof(png_bytep) * );
  for(int y = 0; y < mPngLoader.height; y++) {
      mPngLoader.row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
  }

  png_read_image(png, mPngLoader.row_pointers);

  fclose(fp);

  png_destroy_read_struct(&png, &info, NULL);
}

class Texture {
public:
    ColorRGBA* data;
    
    int width, height;
    int widthm1, heightm1;
    Texture(int width, int height)
    : width(width), height(height), widthm1(width-1), heightm1(height-1) {
        
        data = new ColorRGBA[width*height];
        
    }
    ~Texture() {
        delete [] data;
    }
    
    void set(const double& x, const double& y, const ColorRGBA& value) {
        int xPart = mod(x*(double)width,width);
        int yPart = mod(y*(double)height,height);
        printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
        data[ xPart + width*yPart] = value;
    }
    
    void setPixel(const int& x, const int& y, const ColorRGBA& value) {
        int xPart = mod(x,width);
        int yPart = mod(y,height);
//        printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
        data[ xPart + width*yPart] = value;
    }
    
    void monochromize() {
        for(int i = 0;  i < width*height; i++) {
            double avg = ((double)data[i].r + (double)data[i].g + (double)data[i].b)/3.0; // MONOCHROME
            data[i].r = avg;
            data[i].g = avg;
            data[i].b = avg;
        }
    }
    
    int min(int a, int b) {
        return a < b ? a : b;
    }
    int max(int a, int b) {
        return a > b ? a : b;
    }
    
    void normalize(double factor) {
        int maxi = 0, mini = 255;
        for(int i = 0;  i < width*height; i++) {
            int L = max(max(data[i].r,data[i].g),data[i].b);
            maxi = max(L,maxi);
            mini = min(L,mini);// tricky one here
        }
        double scale = (255.0/(double)(maxi - mini)) * factor + (1.0-factor)*1.0;
        mini = mini * (factor);
        for(int i = 0;  i < width*height; i++) {
            data[i].r = min(max(round(scale * (double)(data[i].r - mini)), 0), 255);
            data[i].g = min(max(round(scale * (double)(data[i].g - mini)), 0), 255);
            data[i].b = min(max(round(scale * (double)(data[i].b - mini)), 0), 255);
            
        }
    }
    
    void offsetAvergageToCenter(double factor) {
        int avg = 0;
        for(int i = 0;  i < width*height; i++) {
            avg += max(max(data[i].r,data[i].g),data[i].b);
        }
        avg = round((double)avg/(double)(width*height) * factor);
        for(int i = 0;  i < width*height; i++) {
            data[i].r = min(max((int)data[i].r - avg, 0), 255);
            data[i].g = min(max((int)data[i].g - avg, 0), 255);
            data[i].b = min(max((int)data[i].b - avg, 0), 255);
        }
    }
    
    ColorRGBA sample( const double& x, const double& y) {
        int xPart = mod(x*(double)width,width);
        int yPart = mod(y*(double)height,height);
        return data[ xPart + width*yPart];
    }
};

void process_png_file(const PngLoader& mPngLoader, Texture& texture) {
  for(int y = 0; y < mPngLoader.height; y++) {
    png_bytep row = mPngLoader.row_pointers[y];
    for(int x = 0; x < mPngLoader.width; x++) {
      png_bytep px = &(row[x * 4]);
      // Do something awesome for each pixel here...
//      printf("%4d, %4d = RGBA(%3d, %3d, %3d, %3d)\n", x, y, px[0], px[1], px[2], px[3]);
        texture.setPixel(x, y, {px[0], px[1], px[2], px[3]});
        
//        double avg = ((double)px[0] + (double)px[1] + (double)px[2])/3.0; // MONOCHROME
////        avg = pow(avg/255.0,1.8)*255;
////        avg = (avg - 205)*255/50;
//        texture.setPixel(x, y, {avg,avg,avg, 0});
        
        
//        double R = ((double)px[0] - 205)*255/50; // STREETMAPS range adjust
//        double G = ((double)px[1] - 205)*255/50; // STREETMAPS range adjust
//        double B = ((double)px[2] - 205)*255/50; // STREETMAPS range adjust
//        texture.setPixel(x, y, {R, G, B, 0});
        
    }
  }
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


typedef struct _UniformInfo {
    Mat4D modelView;
    Mat4D modelViewProjection;
} UniformInfo;



template <class T, class U> void myVertexShader(U* uniformInfo, T& output, T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
    output.normal = matrixVectorMultiply(uniformInfo->modelView, input.normal);
}

typedef struct _CubeVertexInfo {
    Coordinates4D vertex;
    Coordinates4D location;
    Coordinates3D normal;
    Coordinates3D textureCoord;
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
//    colorRGB.x *= 0.5;
//    colorRGB.y *= 0.5;
//    colorRGB.z *= 0.5;
//    colorRGB.x += (double)vertexInfo->color.r * (1.0/255.0 * 0.5);
//    colorRGB.y += (double)vertexInfo->color.g * (1.0/255.0 * 0.5);
//    colorRGB.z += (double)vertexInfo->color.b * (1.0/255.0 * 0.5);
    Coordinates3D clippedRGB = clipRGB(colorRGB);
    fInfo.colorOutput->r = clippedRGB.x*255;
    fInfo.colorOutput->g = clippedRGB.y*255;
    fInfo.colorOutput->b = clippedRGB.z*255;
    fInfo.colorOutput->a = 0;
    
    //    fInfo.colorOutput = result;
}

void textureshader(const FragmentInfo& fInfo) {
    CubeVertexInfo* vertexInfo = (CubeVertexInfo*)fInfo.interpolated;

    Texture* userTexture = (Texture*) fInfo.data;
    *fInfo.colorOutput = userTexture->sample(vertexInfo->textureCoord.x, vertexInfo->textureCoord.y);
    fInfo.colorOutput->a = 0;
}







int main(int argc, char** argv) {
    
    
    if(argc != 2) {
        printf("provide PNG filename as an argument: %s <image.png>\n", argv[0]);
        abort();
    }
      PngLoader mPngLoader;
    read_png_file(argv[1], mPngLoader);
      Texture pngTexture(mPngLoader.width, mPngLoader.height);
    process_png_file(mPngLoader, pngTexture);
  //  write_png_file(argv[2]);
    pngTexture.monochromize();
    pngTexture.offsetAvergageToCenter(0.5);
    pngTexture.normalize(1.);
    
//    testTexture.set(0.2, 0.2, {255,0,0,0});

    setupTerminal();
    
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
    int cubeQuadIndices[][4] = {
        {0, 1, 2, 3},    // left
        {7, 6, 5, 4},    // right
        {6, 2, 1, 5},     // top
        {0, 3, 7, 4},     // bottom
        {2, 6, 7, 3},     // front
        {0, 4, 5, 1}     // back
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
    
    squareVi[0].textureCoord = {0, 0, 0};
    squareVi[1].textureCoord = {0, 1, 0};
    squareVi[2].textureCoord = {1, 1, 0};
    squareVi[3].textureCoord = {1, 0, 0};
    
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
//    double characterAspect = 28.0/14.0; // raspbian terminal
//    double characterAspect = 6.0/4.0; // zipitZ2
    
    int screenSizeX, screenSizeY;
    getmaxyx(stdscr, screenSizeY, screenSizeX);
    
    
    RenderPipeline mRenderPipeline;
    mRenderPipeline.resize(screenSizeX, screenSizeY);
    
    double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
    
    // Depth buffer
//    DepthBuffer depthBuffer;
//    depthBuffer.setSize(screenSizeX, screenSizeY);
    
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
//    Mat4D windowScale = scaleMatrix((double)screenSizeX/2, (double)screenSizeY/2, 1);
//    Mat4D translationScreen = translationMatrix((double)screenSizeX/2 -0.5, (double)screenSizeY/2 -0.5, 0);
//    Mat4D windowFull = matrixMultiply(translationScreen, windowScale);
//    Mat4D windowFull = makeWindowTransform(screenSizeX, screenSizeY, characterAspect);
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
    double delayTime = 0;//1.0/60;
    while (keepRunning == true) {
        debugLine = 0;
        
        now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> float_ms = (now - before);
        before = now;
        
        
        delayTime += 0.001*(1.0/60.0 - float_ms.count()/1000.0);
        if (delayTime> 0 && delayTime < 0.1) {
            usleep(1000000.0*delayTime);
        }
//        depthBuffer.reset();
        mRenderPipeline.reset();
        //erase();
//        mvprintw(debugLine++, 0, "tilt: %f", tilt*180.0/M_PI);
        
        
        
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
//            rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, windowFull, (void*)&mLightParams.color[i], &depthBuffer, lightModelFs, debugLine);
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
//        Coordinates4D cubeTriangle[sizeof(triangle)/sizeof(triangle[0])];
        Mat4D solidCubeScale = scaleMatrix(1.05, 1.05, 1.05);
        Coordinates3D solidAxis = {0,0,1};
        Mat4D solidRotation = rotationFromAngleAndUnitAxis(M_PI_4, solidAxis);
        Mat4D solidTranslation = translationMatrix(0, 6, 0);
        
//        viewMatrix = translationMatrix(0, 0, -2 + sin(angle*4));
        
        Mat4D modelViewBlackCube = matrixMultiply(solidRotation, solidCubeScale);
        modelViewBlackCube = matrixMultiply(solidTranslation, modelViewBlackCube);
//        modelViewBlackCube = matrixMultiply( cube2translation, modelViewBlackCube);
        modelViewBlackCube = matrixMultiply(viewMatrix, modelViewBlackCube);
        numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
//        rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&lightModelView, &depthBuffer, lightFs, debugLine);
//        rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, &depthBuffer, lightFs2, debugLine);
//        mRenderPipeline.setFragmentShader(lightFs2);
//        mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, debugLine);
        
        mUniformInfo.modelView = modelViewBlackCube;
        mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
//        mRenderPipeline.userData = (void*)&mLightParams;
        mRenderPipeline.setFragmentShader(lightFs3);
        mRenderPipeline.trianglesFill(cubeVi, cubeViIndices, 12);
        mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParams, myVertexShader);
//        mRenderPipeline.rasterizeShader(cubeVi, cubeViIndices, 12, modelViewBlackCube, projection, (void*)&mLightParams);
        
//        mRenderPipeline.trianglesFill(squareVi, squareViIndices, 2);
        Mat4D modelTranslation = translationMatrix(0, 0, -1 );
        Mat4D modelScale = scaleMatrix(8, 8, 8);
        Mat4D modelSquare = matrixMultiply(modelTranslation, modelScale);
        Mat4D modelViewSquare = matrixMultiply(viewMatrix, modelSquare);
        mUniformInfo.modelView = modelViewSquare;
        mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
        mRenderPipeline.setFragmentShader(textureshader);
//        mRenderPipeline.rasterizeShader(squareVi, &mUniformInfo, squareViIndices, 2, (void*)&mLightParams, myVertexShader);
        mRenderPipeline.rasterizeShader(squareVi, &mUniformInfo, squareViIndices, 2, (void*)&pngTexture, myVertexShader);
        
        
        if (showDepth) {
            mRenderPipeline.depthBufferToTerminal();
        } else {
            mRenderPipeline.renderBufferToTerminal();
        }
        
        if (autoRotate) {
            angle -= float_ms.count()/1000.0*0.4;
        }

        // HUD
        mvprintw(debugLine++, 0, "FPS: %f", 1000.0/float_ms.count());
        mvprintw(debugLine++, 0, "Delay time %f", delayTime);
        
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
                projection = projectionMatrixPerspective(M_PI*0.5, screenAspect, zFar, zNear);
            } else {
                projection = projectionMatrixOrtho(5*screenAspect, 5, zFar, zNear);
            }
//            depthBuffer.setSize(screenSizeX, screenSizeY);
        } else if (ch == 'o' || ch == 'O') {
            usePerspective    = !usePerspective;
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
    
    cleanupConsole();
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



