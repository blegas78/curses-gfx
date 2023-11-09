#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <unistd.h>

#include <ncurses.h>
#include <termios.h>

#include <time.h>

#ifdef FB_SUPPORT
#include <linux/fb.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <sys/mman.h>
#endif


//#include "curses-gfx-3d.h"
#include "curses-gfx-handler.h"
#include "curses-gfx-texture.h"

/*
 Catch ctrl-c for cleaner exits
 */
static volatile bool keepRunning = true;
void killPanda(int killSignal) {
	keepRunning = false;
}


#ifdef FB_SUPPORT
	int fbfd;
	struct fb_var_screeninfo vinfo;
	struct fb_fix_screeninfo finfo;
	int fb_width;
	int fb_height;
	int fb_bpp;
	int fb_bytes;
	int fb_bytes_per_length;
	
	int fb_data_size;

	char *fbdata;
	
void setupLinuxFb() {
	fbfd = open ("/dev/fb0", O_RDWR);
	if (fbfd < 0) {
		printf("Unable to open /dev/fb0!\n");
		//return 1;
	}
	
	ioctl (fbfd, FBIOGET_VSCREENINFO, &vinfo);
	ioctl (fbfd, FBIOGET_FSCREENINFO, &finfo);

	fb_width = vinfo.xres;
	fb_height = vinfo.yres;
	fb_bpp = vinfo.bits_per_pixel;
	fb_bytes = fb_bpp / 8;
	fb_bytes_per_length = finfo.line_length/fb_bytes;
	
//	fb_data_size = fb_width * fb_height * fb_bytes * fb_bytes_per_length;
	fb_data_size = fb_height * fb_bytes_per_length* fb_bytes;

	fbdata = (char*)mmap (0, fb_data_size, PROT_READ | PROT_WRITE, MAP_SHARED, fbfd, (off_t)0);
	
	printf("Framebuffer info: res %dx%d, bpp=%d, Bytes=%d\n", fb_width, fb_height, fb_bpp, fb_bytes);
}
#endif


typedef struct _MyShaderAttributes {
    Coordinates4D vertex;
    ColorRGBA color;
//    ColorRGBA colorOutput;
    Coordinates3D textureCoord;
} MyShaderAttributes;

REGISTER_VERTEX_LAYOUT(MyShaderAttributes)
//    MEMBER(vertex),
    MEMBER(color),
    MEMBER(textureCoord)
END_VERTEX_LAYOUT(MyShaderAttributes)


Texture testTexture(4,4);
void shaderBasic(const FragmentInfo2& fInfo) {
    Coordinates3D* colorRGB = (Coordinates3D*)fInfo.data;
    MyShaderAttributes* mMyShaderAttributes = (MyShaderAttributes*)fInfo.interpolated;
    //setRGB(fInfo.pixel, *colorRGB);
    
    fInfo.colorOutput->r = mMyShaderAttributes->color.r;
    fInfo.colorOutput->g = mMyShaderAttributes->color.g;
    fInfo.colorOutput->b = mMyShaderAttributes->color.b;
    fInfo.colorOutput->a = 0;
    
    *fInfo.colorOutput = testTexture.sample(mMyShaderAttributes->textureCoord.x, mMyShaderAttributes->textureCoord.y)*0.5 + (*fInfo.colorOutput)*0.5;
    
//    fInfo.colorOutput->r = pow((double)fInfo.colorOutput->r/255.0, 2.)*255;
//    fInfo.colorOutput->g = pow((double)fInfo.colorOutput->g/255.0, 2.)*255;
//    fInfo.colorOutput->b = pow((double)fInfo.colorOutput->b/255.0, 2.)*255;
//    *fInfo.colorOutput = mMyShaderAttributes->color;
}

int main(void) {
#ifdef FB_SUPPORT
	setupLinuxFb();
#endif
	
    
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

    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();
	

	Coordinates4D points[6];
    
    MyShaderAttributes mMyShaderAttributes[6];
	
	
	for (int i = 0; i < 6; i++) {
		points[i].x = 63*cos(3.14159/3.0 * (double)i);
		points[i].y = -63*sin(3.14159/3.0 * (double)i);
        
        mMyShaderAttributes[i].textureCoord.x = 10+1.*cos(3.14159/3.0 * (double)i);
        mMyShaderAttributes[i].textureCoord.y = 10+1.*sin(3.14159/3.0 * (double)i);
		
        points[i].z = 0;
        points[i].w = 1;
		
        mMyShaderAttributes[i].vertex.x = points[i].x;
        mMyShaderAttributes[i].vertex.y = points[i].y;
        mMyShaderAttributes[i].vertex.z = points[i].z;
	}
	
//	Polygon4D poly;
//	poly.numVertices = 6;
	for (int i = 0; i < 6; i++) {
		
		Coordinates3D hsv, rgb;
		hsv.x = i*60;
		hsv.y = 0.9;
		hsv.z = 0.5;
		rgb = hslToRgb(hsv);
		
        
        mMyShaderAttributes[i].color.r = rgb.x;
        mMyShaderAttributes[i].color.g = rgb.y;
        mMyShaderAttributes[i].color.b = rgb.z;
        mMyShaderAttributes[i].color.a = 0;
	}
    RenderPipeline mRenderPipeline;
    
	
	FrameBuffer renderBuffer;
	
	asTexImage2d(&renderBuffer, FBT_RGBA, 64, 64);
    mRenderPipeline.resize(128, 128);
	
	ColorRGBA clearColor = {0,0,0,0};//{127,127,127,0};
	ColorRGBA whiteColor = {255,255,255,0};//{127,127,127,0};
	
	double angle = 0;
	bool autorotate = true;
	int ch;
	bool showPoints = true;
	while (keepRunning == true) {
		erase();
//		renderBuffer.clear(clearColor);
        mRenderPipeline.reset();
	
		Coordinates3D axis = {0,0,1};
		Mat4D rotation = rotationFromAngleAndUnitAxis(angle, axis);
        Coordinates3D axis2 = {0,1,0};
        Mat4D rotation2 = rotationFromAngleAndUnitAxis(sin(angle/2.1)*M_PI/2, axis2);
        Mat4D scale = scaleMatrix(1,12.0/28.0,1);
		Mat4D translation = translationMatrix(64, 64, 128);
		Mat4D model = matrixMultiply(rotation2, rotation);
        model = matrixMultiply(translation, model);
        model = matrixMultiply(scale, model);
        Mat4D perspective = projectionMatrixPerspective(M_PI_2, 12.0/28.0, 100, 0.01);
//		Coordinates4D p[3];
		for(int i = 0; i < 6; i++) {
//			fragments[i].location3D = matrixVectorMultiply(model, points[i]);
//			poly.vertices[i] = matrixVectorMultiply(model, points[i]);
            mMyShaderAttributes[i].vertex = matrixVectorMultiply(model, points[i]);
//            mMyShaderAttributes[i].vertex = matrixVectorMultiply(perspective, mMyShaderAttributes[i].vertex);
		}

        mRenderPipeline.setFragmentShader(shaderBasic);
        
        int hexagonIndes[][3] = {
            {0, 1, 2},
            {2, 3, 4},
            {5, 0, 2},
            {2, 4, 5}
        };
        
//        mRenderPipeline.trianglesFill(mMyShaderAttributes, hexagonIndes, 4);
        for(int i = 0; i < 4; i++) {
            mRenderPipeline.triangleFill(&mMyShaderAttributes[hexagonIndes[i][0]], &mMyShaderAttributes[hexagonIndes[i][1]], &mMyShaderAttributes[hexagonIndes[i][2]], NULL);
        }
		
        mRenderPipeline.renderBufferToTerminal();
		
		
		
		usleep(1000000.0/60.0);
		
		if (autorotate) {
			angle += 0.01;
		}
		
		if ((ch = getch()) == 0x1B) {	// Escape
			keepRunning = false;
		} else if ( ch == KEY_LEFT) {
			angle += 0.05;
		} else if ( ch == KEY_RIGHT) {
			angle -= 0.05;
		} else if (ch == ' ') {
			autorotate = !autorotate;
		}

	}
	
	
	return 0;
}
