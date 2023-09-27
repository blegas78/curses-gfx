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
	init_pair(2, COLOR_YELLOW, COLOR_BLACK);
	init_pair(3, COLOR_GREEN, COLOR_BLACK);
	init_pair(4, COLOR_CYAN, COLOR_BLACK);
	init_pair(5, COLOR_BLUE, COLOR_BLACK);
	init_pair(6, COLOR_MAGENTA, COLOR_BLACK);
	init_pair(7, COLOR_WHITE, COLOR_BLACK);


//	init_pair(5, COLOR_BLACK, COLOR_RED );
//	init_pair(6, COLOR_BLACK, COLOR_GREEN );
//	init_pair(7, COLOR_BLACK, COLOR_CYAN );
//	init_pair(8, COLOR_WHITE, COLOR_BLUE );

	curs_set(0);	// no cursor

//	atexit(destroy);
}

void cleanupConsole() {
	clear();
	endwin();

	std::cout << "Console has been cleaned!" << std::endl;
}




void renderBufferToTerminal(FrameBuffer* fbo) {
	
	Coordinates2D pixel;
	Coordinates3D color;
	int offset, index;
	
	for (int y = 0; y < fbo[0].rows; y++) {
		offset = y * fbo[0].cols;
		for (int x = 0; x < fbo[0].cols; x++) {
			index = x + offset;
			
			pixel.x = x;
			pixel.y = y;
			
			color.x = (double)((ColorRGBA*)fbo[0].data)[index].r / 255.0;
			color.y = (double)((ColorRGBA*)fbo[0].data)[index].g / 255.0;
			color.z = (double)((ColorRGBA*)fbo[0].data)[index].b / 255.0;
//#ifndef FB_SUPPORT
			if (((ColorRGBA*)fbo[0].data)[index].a == 0) {
				setRGB(pixel, color);
			} else {
				Coordinates3D clippedRGB = clipRGB(color);
				Coordinates3D hsl = rgbToHsv(clippedRGB);
				
				int hueIndex = floor(hsl.x + 1);
				
				if (hsl.y < 0.33) {
					hueIndex = 7;
				}
				
				attron(COLOR_PAIR(hueIndex));
				set(pixel, ((ColorRGBA*)fbo[0].data)[index].a);
				attroff(COLOR_PAIR(hueIndex));
			}
#ifdef FB_SUPPORT
			int offsetfb = (y * (fb_bytes_per_length) + x)*fb_bytes;
//			uint16_t finalcolor = 0;
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].r & 0xF8) << (11-3);
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].g & 0xFC) << (5-2);
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].b & 0xF8) >> (3);
//			*(uint16_t*)&fbdata[0 + offsetfb] = finalcolor;
			
			fbdata[2 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].r;
			fbdata[1 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].g;
			fbdata[0 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].b;
			
#endif
		}
	}
}



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
    
    ColorRGBA sample( const double& x, const double& y) {
        int xPart = mod(x*(double)width,width);
        int yPart = mod(y*(double)height,height);
        return data[ xPart + width*yPart];
    }
};



Texture testTexture(4,4);
void shaderBasic(const FragmentInfo& fInfo) {
    Coordinates3D* colorRGB = (Coordinates3D*)fInfo.data;
    MyShaderAttributes* mMyShaderAttributes = (MyShaderAttributes*)fInfo.interpolated;
    //setRGB(fInfo.pixel, *colorRGB);
    
    fInfo.colorOutput->r = mMyShaderAttributes->color.r;
    fInfo.colorOutput->g = mMyShaderAttributes->color.g;
    fInfo.colorOutput->b = mMyShaderAttributes->color.b;
    fInfo.colorOutput->a = 0;
    
    *fInfo.colorOutput = testTexture.sample(mMyShaderAttributes->textureCoord.x, mMyShaderAttributes->textureCoord.y)*0.5 + (*fInfo.colorOutput)*0.5;
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

//    ColorRGBA color = {0,0,255,0};
//    testTexture.set(1.0,1.0, color);
//    testTexture.set(0.0,0.0, color);
//    testTexture.set(0.75,0.75, color);
    
    
    setupTerminal();
	
//	FragmentInfo fragments[3];

	Coordinates4D points[6];
    
    MyShaderAttributes mMyShaderAttributes[6];
	
//	points[0].x = 3*3 - 32;
//	points[0].y = 3*3 - 32;
//
//	points[2].x = 12.*3 - 32;
//	points[2].y = 15.*3 - 32;
//
//	points[1].x = 1.*3 - 32;
//	points[1].y = 15.*3 - 32;
//
//
//	points[3].x = 15.*3 - 32;
//	points[3].y = 1.*3 - 32;
	
	for (int i = 0; i < 6; i++) {
		points[i].x = 63*cos(3.14159/3.0 * (double)i);
		points[i].y = -63*sin(3.14159/3.0 * (double)i);
        
        mMyShaderAttributes[i].textureCoord.x = 10+1.*cos(3.14159/3.0 * (double)i);
        mMyShaderAttributes[i].textureCoord.y = 10+1.*sin(3.14159/3.0 * (double)i);
		
        points[i].z = 1.0;
		
        mMyShaderAttributes[i].vertex.x = points[i].x;
        mMyShaderAttributes[i].vertex.y = points[i].y;
        mMyShaderAttributes[i].vertex.z = points[i].z;
	}
	
//    mMyShaderAttributes[0].textureCoord.x = 1.5;
//    mMyShaderAttributes[0].textureCoord.y = 0;
//
//    mMyShaderAttributes[1].textureCoord.x = 1;
//    mMyShaderAttributes[1].textureCoord.y = 1;
////    mMyShaderAttributes[1].vertex.z = 100;
//
//    mMyShaderAttributes[2].textureCoord.x = 0;
//    mMyShaderAttributes[2].textureCoord.y = 1;
//	for (int i = 0; i < 3; i++) {
//		fragments[i].color.r = 0;
//		fragments[i].color.g = 0;
//		fragments[i].color.b = 0;
//		fragments[i].color.a = 0;
//	}
//	fragments[0].color.r = 255;
////	fragments[0].color.g = 255;
//
//	fragments[1].color.g = 255;
////	fragments[1].color.b = 255;
//
//	fragments[2].color.b = 255;
////	fragments[2].color.r = 255;
//
////	fragments[2].location3D.y = fragments[1].location3D.y;
	
	Polygon4D poly;
	poly.numVertices = 6;
	for (int i = 0; i < poly.numVertices; i++) {
		poly.colors[i].r = 0;
		poly.colors[i].g = 0;
		poly.colors[i].b = 0;
		poly.colors[i].a = 0;
		
		Coordinates3D hsv, rgb;
		hsv.x = i*60;
		hsv.y = 0.9;
		hsv.z = 0.5;
		rgb = hslToRgb(hsv);
		
		poly.colors[i].r = rgb.x;
		poly.colors[i].g = rgb.y;
		poly.colors[i].b = rgb.z;
		poly.colors[i].a = 0;
        
        mMyShaderAttributes[i].color.r = rgb.x;
        mMyShaderAttributes[i].color.g = rgb.y;
        mMyShaderAttributes[i].color.b = rgb.z;
        mMyShaderAttributes[i].color.a = 0;
	}
//	poly.colors[0].r = 255;
//	poly.colors[2].b = 255;
//	poly.colors[1].r = 255;
//	poly.colors[1].b = 255;
//	poly.colors[1].g = 255;
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
		renderBuffer.clear(&clearColor);
        mRenderPipeline.reset();
	
		Coordinates3D axis = {0,0,1};
		Mat4D rotation = rotationFromAngleAndUnitAxis(angle, axis);
        Mat4D scale = scaleMatrix(1,12.0/28.0,1);
		Mat4D translation = translationMatrix(64, 64, 0);
		Mat4D model = matrixMultiply(translation, rotation);
        model = matrixMultiply(scale, model);
//		Coordinates4D p[3];
		for(int i = 0; i < poly.numVertices; i++) {
			points[i].w = 1;
//			fragments[i].location3D = matrixVectorMultiply(model, points[i]);
			poly.vertices[i] = matrixVectorMultiply(model, points[i]);
            mMyShaderAttributes[i].vertex = poly.vertices[i];
		}

        mRenderPipeline.setFragmentShader(shaderBasic);
//		triangleFill(fragments, &renderBuffer);
//		drawPolygonWithTriangles(poly, &renderBuffer);
//        mRenderPipeline.drawPolygonWithTriangles(poly, poly, NULL);
        
//        mRenderPipeline.triangleFill(&mMyShaderAttributes[0], &mMyShaderAttributes[1], &mMyShaderAttributes[2] );
//        mRenderPipeline.triangleFill(&mMyShaderAttributes[2], &mMyShaderAttributes[3], &mMyShaderAttributes[4] );
//        mRenderPipeline.triangleFill(&mMyShaderAttributes[5], &mMyShaderAttributes[0], &mMyShaderAttributes[2] );
//        mRenderPipeline.triangleFill(&mMyShaderAttributes[2], &mMyShaderAttributes[4], &mMyShaderAttributes[5] );
        int hexagonIndes[][3] = {
            {0, 1, 2},
            {2, 3, 4},
            {5, 0, 2},
            {2, 4, 5}
        };
        
        mRenderPipeline.trianglesFill(mMyShaderAttributes, hexagonIndes, 4);
		
//		if(showPoints)
//		for (int i = 0; i < poly.numVertices; i++) {
//			((ColorRGBA*)renderBuffer.data)[(int)poly.vertices[i].y * renderBuffer.cols + (int)poly.vertices[i].x] = whiteColor;
//		}
//		showPoints = !showPoints;
//		renderBufferToTerminal(&renderBuffer);
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
	
	
	
	cleanupConsole();

	
	
	
	
	return 0;
}
