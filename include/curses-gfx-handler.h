#ifndef CURSES_GFX_HANDLER_H
#define CURSES_GFX_HANDLER_H

#include <pthread.h>
#ifdef FB_SUPPORT
#include <linux/fb.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <sys/mman.h>
#endif

#include "curses-gfx.h"
#include "curses-gfx-3d.h"



typedef struct _Renderable {
	Polygon4D* polygons;
	int count;
	Mat4D& modelView;
	Mat4D& projection;
	Mat4D& viewPort;
	void* userData;
	DepthBuffer* depthBuffer;
	void (*fragmentShader)(const FragmentInfo&);
	int &line;
} Renderable;






class CursesGfxHandler {
private:
	pthread_t* threads;
	pthread_mutex_t mutex;
	pthread_mutex_t condMutex;
	pthread_cond_t cond;
	
	int totalThreads;
	bool busyThreads;
	
public:
	CursesGfxHandler();
	~CursesGfxHandler();
	
	void setTotalThreads(int numthreads);
	
};


class RenderPipeline {
private:
	FrameBuffer* fbo;
	FrameBuffer* depthBuffer;
	
	double depthClearColor;
	
//	DepthBuffer tempD;
	
	void (*fragmentShader)(const FragmentInfo&);
	
	void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, int &line);
	
	void setWithDepthBuffer( Coordinates2D pt, char c, double depth);
	void setWithDepthBuffer( Coordinates4D pt, char c, double depth);
	
	void drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData);
	
	void setFloatDotWithDepthBuffer( double x, double y, double depth);
	void drawDotFloat(double x, double y);
	
	void setWithShader( Coordinates2D& pixel, double invDepth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData);
	
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
	
	void setupLinuxFb();
#endif
public:
	RenderPipeline();
	~RenderPipeline();
	
	void resize(int width, int height);
	void reset(); // resets depth buffer
	
	void setFragmentShader(void (*fragmentShader)(const FragmentInfo&));
	
	void rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	void rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	
	void setRenderBuffer(int x, int y, ColorRGBA& color);
	
	void renderBufferToTerminal();
	void depthBufferToTerminal();
};

#endif
