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


//class Texture {
//public:
//    ColorRGBA* data;
//
//    int width, height;
//    int widthm1, heightm1;
//    Texture(int width, int height)
//    : width(width), height(height), widthm1(width-1), heightm1(height-1) {
//
//        data = new ColorRGBA[width*height];
//
//    }
//    ~Texture() {
//        delete [] data;
//    }
//
//    void set(const double& x, const double& y, const ColorRGBA& value) {
//        int xPart = mod(x*(double)width,width);
//        int yPart = mod(y*(double)height,height);
//        printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
//        data[ xPart + width*yPart] = value;
//    }
//
//    void setPixel(const int& x, const int& y, const ColorRGBA& value) {
//        int xPart = mod(x,width);
//        int yPart = mod(y,height);
////        printf("Index xPart %d yPart %d result%d\n", xPart, yPart, xPart + width*yPart);
//        data[ xPart + width*yPart] = value;
//    }
//
//    void monochromize() {
//        for(int i = 0;  i < width*height; i++) {
//            double avg = ((double)data[i].r + (double)data[i].g + (double)data[i].b)/3.0; // MONOCHROME
//            data[i].r = avg;
//            data[i].g = avg;
//            data[i].b = avg;
//        }
//    }
//
//    int min(int a, int b) {
//        return a < b ? a : b;
//    }
//    int max(int a, int b) {
//        return a > b ? a : b;
//    }
//
//    void normalize(double factor) {
//        int maxi = 0, mini = 255;
//        for(int i = 0;  i < width*height; i++) {
//            int L = max(max(data[i].r,data[i].g),data[i].b);
//            maxi = max(L,maxi);
//            mini = min(L,mini);// tricky one here (min of maxes)
//        }
//        double scale = (255.0/(double)(maxi - mini)) * factor + (1.0-factor)*1.0;
//        mini = mini * (factor);
//        for(int i = 0;  i < width*height; i++) {
//            data[i].r = min(max(round(scale * (double)(data[i].r - mini)), 0), 255);
//            data[i].g = min(max(round(scale * (double)(data[i].g - mini)), 0), 255);
//            data[i].b = min(max(round(scale * (double)(data[i].b - mini)), 0), 255);
//
//        }
//    }
//
//    void offsetAvergageToCenter(double factor) {
//        int avg = 0;
//        for(int i = 0;  i < width*height; i++) {
//            avg += max(max(data[i].r,data[i].g),data[i].b);
//        }
//        avg = round((double)avg/(double)(width*height) * factor);
//        for(int i = 0;  i < width*height; i++) {
//            data[i].r = min(max((int)data[i].r - avg, 0), 255);
//            data[i].g = min(max((int)data[i].g - avg, 0), 255);
//            data[i].b = min(max((int)data[i].b - avg, 0), 255);
//        }
//    }
//
//    ColorRGBA sample( const double& x, const double& y) {
//        int xPart = mod(x*(double)width,width);
//        int yPart = mod(y*(double)height,height);
//        return data[ xPart + width*yPart];
//    }
//};












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



template <class T, class U> void myVertexShader(U* uniformInfo, T& output, const T& input) {
    output.vertex = matrixVectorMultiply(uniformInfo->modelViewProjection, input.vertex);
    output.location = matrixVectorMultiply(uniformInfo->modelView, input.vertex);
    output.normal = matrixVectorMultiply(uniformInfo->modelView, input.normal);
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
    
    
    
    //setRGB(fInfo.pixel, colorRGB);
    
    //    ColorRGBA result;
//    colorRGB.x *= 0.5;
//    colorRGB.y *= 0.5;
//    colorRGB.z *= 0.5;
//    colorRGB.x += (double)vertexInfo->color.r * (1.0/255.0 * 0.5);
//    colorRGB.y += (double)vertexInfo->color.g * (1.0/255.0 * 0.5);
//    colorRGB.z += (double)vertexInfo->color.b * (1.0/255.0 * 0.5);
    
//    colorRGB.x = (colorRGB.x * (double)textureSample.r)*0.8/255 + (double)fInfo.colorOutput->r*0.2/255;
//    colorRGB.y = (colorRGB.y * (double)textureSample.g)*0.8/255 + (double)fInfo.colorOutput->g*0.2/255;
//    colorRGB.z = (colorRGB.z * (double)textureSample.b)*0.8/255 + (double)fInfo.colorOutput->b*0.2/255;
    Coordinates3D clippedRGB = clipRGB(colorRGB);
    fInfo.colorOutput->r = clippedRGB.x*255;
    fInfo.colorOutput->g = clippedRGB.y*255;
    fInfo.colorOutput->b = clippedRGB.z*255;
    fInfo.colorOutput->a = 0;
    
    //    fInfo.colorOutput = result;
}







int main(int argc, char** argv) {
    
    
    if(argc != 2) {
        printf("provide PNG filename as an argument: %s <image.png>\n", argv[0]);
        abort();
    }
//      PngLoader mPngLoader;
//    read_png_file(argv[1], mPngLoader);
//      Texture pngTexture(mPngLoader.width, mPngLoader.height);
//    process_png_file(mPngLoader, pngTexture);
    
//    Texture pngTexture(argv[1]);
    Texture pngTexture;
    pngTexture.loadPng(argv[1]);
//    return 0;
  //  write_png_file(argv[2]);
//    pngTexture.monochromize();
//    pngTexture.offsetAvergageToCenter(0.5);
//    pngTexture.normalize(1.);
//    pngTexture.gaussBlur();
//    pngTexture.gaussBlur();
//    pngTexture.gaussBlur();
//    pngTexture.gaussBlur();
//    pngTexture.gaussBlur();
//    pngTexture.gaussBlur();
    
    
    
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
    cubeVi[0].textureCoord = {0, 0};
    cubeVi[1].textureCoord = {0, 1};
    cubeVi[2].textureCoord = {1, 1};
    cubeVi[3].textureCoord = {1, 1};
    cubeVi[4].textureCoord = {1, 0};
    cubeVi[5].textureCoord = {0, 0};
    
    cubeVi[0+6].textureCoord = {0, 0};
    cubeVi[1+6].textureCoord = {1, 1};
    cubeVi[2+6].textureCoord = {1, 0};
    cubeVi[3+6].textureCoord = {1, 1};
    cubeVi[4+6].textureCoord = {0, 0};
    cubeVi[5+6].textureCoord = {0, 1};
    
    cubeVi[0+12].textureCoord = {1, 0};
    cubeVi[1+12].textureCoord = {0, 0};
    cubeVi[2+12].textureCoord = {1, 1};
    cubeVi[3+12].textureCoord = {0, 0};
    cubeVi[4+12].textureCoord = {0, 1};
    cubeVi[5+12].textureCoord = {1, 1};
    
    cubeVi[0+18].textureCoord = {0, 1};
    cubeVi[1+18].textureCoord = {1, 1};
    cubeVi[2+18].textureCoord = {0, 0};
    cubeVi[3+18].textureCoord = {0, 0};
    cubeVi[4+18].textureCoord = {1, 1};
    cubeVi[5+18].textureCoord = {1, 0};
    
    cubeVi[0+24].textureCoord = {1, 0};
    cubeVi[1+24].textureCoord = {0, 0};
    cubeVi[2+24].textureCoord = {1, 1};
    cubeVi[3+24].textureCoord = {0, 1};
    cubeVi[4+24].textureCoord = {1, 1};
    cubeVi[5+24].textureCoord = {0, 0};
    
    cubeVi[0+30].textureCoord = {0, 0};
    cubeVi[1+30].textureCoord = {0, 1};
    cubeVi[2+30].textureCoord = {1, 1};
    cubeVi[3+30].textureCoord = {1, 1};
    cubeVi[4+30].textureCoord = {1, 0};
    cubeVi[5+30].textureCoord = {0, 0};
    
    // texture wrap:
    for(int r = 0; r < 2; r++) {
        for(int c = 0; c < 3; c++) {
            for(int v = 0; v < 6; v++) {
                Coordinates2Df offset ;
                offset.x = 0.25*(double)c;
                offset.y = 0.5*(double)r;
                cubeVi[v + 6*(r + c*2)].textureCoord.x *= 0.25;
                cubeVi[v + 6*(r + c*2)].textureCoord.y *= 0.5;
                cubeVi[v + 6*(r + c*2)].textureCoord = cubeVi[v + 6*(r + c*2)].textureCoord + offset;
            }
        }
    }
    
    CubeVertexInfo squareVi[4];
    int squareViIndices[2][3];
    squareVi[0].vertex = {-1, -1, 0, 1};
    squareVi[1].vertex = { 1, -1, 0, 1};
    squareVi[2].vertex = { 1,  1, 0, 1};
    squareVi[3].vertex = {-1,  1, 0, 1};
    
    squareVi[0].textureCoord = {0.75, 0};
    squareVi[1].textureCoord = {1, 0};
    squareVi[2].textureCoord = {1, 0.5};
    squareVi[3].textureCoord = {0.75, 0.5};
    
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
    
    CubeVertexInfo squareWall[4];
    squareWall[0].vertex = {-1, -1, 0, 1};
    squareWall[1].vertex = { 1, -1, 0, 1};
    squareWall[2].vertex = { 1,  1, 0, 1};
    squareWall[3].vertex = {-1,  1, 0, 1};
    
    squareWall[0].textureCoord = {0.75, 0.5}; // this makes this square different
    squareWall[1].textureCoord = {1, 0.5};
    squareWall[2].textureCoord = {1, 1};
    squareWall[3].textureCoord = {0.75, 1};
    
    for(int i = 0; i < 4; i++) {
        squareVi[i].normal = {0,0,1};
        squareWall[i].normal = {0,0,1};
    }
    
    double characterAspect = 28.0/12.0; // macOs terminal
//    double characterAspect = 28.0/14.0; // raspbian terminal
//    double characterAspect = 6.0/4.0; // zipitZ2
    
    int screenSizeX, screenSizeY;
    getmaxyx(stdscr, screenSizeY, screenSizeX);
    
    RasterizerThreadPool::setRenderThreadCount(4);
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
    double zFar = 100;
    double zNear = .1;
    double viewAngle = M_PI*0.5;
    double orthoScale = 5*5;
    Mat4D projection = projectionMatrixPerspective(viewAngle, screenAspect, zFar, zNear);
    
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
    auto before2 = now;
    
    RenderStats mRenderStats;
    
    int debugLine = 0;
    int numEdges;
    double lightAngle = 0;
    double angle = 0;
    double cubeAngle = 0;
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
//            usleep(1000000.0*delayTime);
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
            numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
            
            lightModelView = matrixVectorMultiply(viewMatrix, light[i]);
//            rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, windowFull, (void*)&mLightParams.color[i], &depthBuffer, lightModelFs, debugLine);
            mRenderPipeline.setFragmentShader(lightModelFs);

            mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, lightCubeModelView, projection, mRenderPipeline.viewport, (void*)&mLightParams.color[i], debugLine);
        }

        LightParamsAndTexture mLightParamsAndTexture;
        mLightParamsAndTexture.lightParams = &mLightParams;
        mLightParamsAndTexture.texture = &pngTexture;
        Coordinates4D zeroVector = {0,0,0,1};
        Coordinates4D camLocation = matrixVectorMultiply(cameraTranslation, zeroVector);
//        mLightParamsAndTexture.cameraLocation = {camLocation.x, camLocation.y, camLocation.z};
        mLightParamsAndTexture.cameraLocation = {0,0,0};
//        mLightParamsAndTexture.cameraLocation = matrixVectorMultiple(<#Mat3D &rotation#>, -1*mLightParamsAndTexture.cameraLocation);
        
        now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> float_ms2 = (now - before);
        double timeToPrepare = (double)float_ms2.count()/1000.0;
        before2 = now;
        

        // triangle
//        Coordinates4D cubeTriangle[sizeof(triangle)/sizeof(triangle[0])];
        Mat4D solidCubeScale = scaleMatrix(2.05, 2.05, 2.05);
        Coordinates3D solidAxis = {sin(cubeAngle*.092440),sin(cubeAngle*.923840),sin(cubeAngle*1.123)};
        solidAxis = normalizeVector(solidAxis);
        Mat4D solidRotation = rotationFromAngleAndUnitAxis(M_PI_4+cubeAngle, solidAxis);
        Mat4D solidTranslation = translationMatrix(0, 0, 0);
        
//        viewMatrix = translationMatrix(0, 0, -2 + sin(angle*4));
        
        Mat4D modelBlackCube = matrixMultiply(solidRotation, solidCubeScale);
        modelBlackCube = matrixMultiply(solidTranslation, modelBlackCube);
        
        Mat4D modelViewBlackCube = matrixMultiply(viewMatrix, modelBlackCube);
        numEdges = sizeof(cubeQuadIndices)/sizeof(cubeQuadIndices[0]);
//        rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&lightModelView, &depthBuffer, lightFs, debugLine);
//        rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, &depthBuffer, lightFs2, debugLine);
//        mRenderPipeline.setFragmentShader(lightFs2);
//        mRenderPipeline.rasterizeQuadsShader(cube, cubeQuadIndices, numEdges, modelViewBlackCube, projection, windowFull, (void*)&mLightParams, debugLine);
        
        mUniformInfo.modelView = modelViewBlackCube;
//        mUniformInfo.modelView = modelBlackCube;
        mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
//        mRenderPipeline.setFragmentShader(lightFs3);
//        mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParams, myVertexShader);
        mRenderPipeline.setFragmentShader(textureshader);
        mRenderStats = mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&pngTexture, myVertexShader);
//        mRenderPipeline.setFragmentShader(lightAndTextureShader);
//        mRenderPipeline.rasterizeShader(cubeVi, &mUniformInfo, cubeViIndices, 12, (void*)&mLightParamsAndTexture, myVertexShader);
        
        if(showGround) {
//            for (int i = 0; i < 4; i++) {
                
                // Square floor:
                //        mRenderPipeline.trianglesFill(squareVi, squareViIndices, 2);
                double s = 10;
//                Mat4D modelTranslation = translationMatrix(s*cos((double)i/2.0*M_PI + M_PI_4)*sqrt(2), s*sin((double)i/2.0*M_PI + M_PI_4)*sqrt(2), -3 );
            Mat4D modelTranslation = translationMatrix(0, 0, -s );
                Mat4D modelScale = scaleMatrix(s, s, s);
                Mat4D modelSquare = matrixMultiply(modelTranslation, modelScale);
                Mat4D modelViewSquare = matrixMultiply(viewMatrix, modelSquare);
                mUniformInfo.modelView = modelViewSquare;
                //        mUniformInfo.modelView = modelSquare;
                mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
                
                
                
                //        mRenderPipeline.setFragmentShader(textureshader);
                //        mRenderStats = mRenderPipeline.rasterizeShader(squareVi, &mUniformInfo, squareViIndices, 2, (void*)&pngTexture, myVertexShader);
                mRenderPipeline.setFragmentShader(lightAndTextureShader);
            RenderStats mRenderStats2 = mRenderPipeline.rasterizeShader(squareVi, &mUniformInfo, squareViIndices, 2, (void*)&mLightParamsAndTexture, myVertexShader);
                
                mRenderStats.timeDrawing += mRenderStats2.timeDrawing;
                mRenderStats.timeClipping += mRenderStats2.timeClipping;
                mRenderStats.timeVertexShading += mRenderStats2.timeVertexShading;
//            }
            
            
            Coordinates3D wallAxis = {1, 0, 0};
            Mat4D wallScale = scaleMatrix(s, s, s);
            Mat4D wallRotation = rotationFromAngleAndUnitAxis(-M_PI_2, wallAxis);
            Mat4D wallTranslation = translationMatrix(0, -s, 0 );
            Mat4D firstWallModel = matrixMultiply(wallRotation, wallScale);
            firstWallModel = matrixMultiply(wallTranslation, firstWallModel);
            mRenderPipeline.setFragmentShader(lightAndTextureShader);
            wallAxis = {0, 0, 1};
            for(int i = 0; i < 4; i++) {
                Mat4D thisWallRotation = rotationFromAngleAndUnitAxis(M_PI_2 *(double)i, wallAxis);
                Mat4D thisWallModel = matrixMultiply(thisWallRotation, firstWallModel);
                
                Mat4D thisWallModelView = matrixMultiply(viewMatrix, thisWallModel);
                mUniformInfo.modelView = thisWallModelView;
                mUniformInfo.modelViewProjection = matrixMultiply(projection, mUniformInfo.modelView);
                mRenderStats2 = mRenderPipeline.rasterizeShader(squareWall, &mUniformInfo, squareViIndices, 2, (void*)&mLightParamsAndTexture, myVertexShader);
                mRenderStats.timeDrawing += mRenderStats2.timeDrawing;
                mRenderStats.timeClipping += mRenderStats2.timeClipping;
                mRenderStats.timeVertexShading += mRenderStats2.timeVertexShading;
            }
            
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
        cubeAngle -= dTime*0.4;
        
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
    
    return 0;
};



