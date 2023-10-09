#ifndef CURSES_GFX_3D_H
#define CURSES_GFX_3D_H

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <mutex>

#include "curses-gfx.h"

// Neat way to view template/macro expansions: https://cppinsights.io
#define REGISTER_VERTEX_LAYOUT(Class) \
auto regist(Class& value) \
{ \
    using Type = Class; \
    return std::make_tuple(
 
#define MEMBER(Name) &value.Name
 
#define END_VERTEX_LAYOUT(Class) \
    ); \
};


// There's a great visual rendering pipeline visual here (helped me with depth buffer): https://fabiensanglard.net/polygon_codec/


//typedef struct _Coordinates2Df {
//    double x;
//    double y;
//
//    _Coordinates2Df() {}
//    _Coordinates2Df( const double& x, const double& y) : x(x), y(y) {}
//    _Coordinates2Df( const _Coordinates2Df& p) : x(p.x), y(p.y) {}
//    _Coordinates2Df operator + (const _Coordinates2Df &p) const { return _Coordinates2Df(p.x+x, p.y+y); }
//    _Coordinates2Df operator * (const double &d) const { return _Coordinates2Df(d*x, d*y); }
//} Coordinates2Df;
//
//typedef struct _Coordinates3D {
//    double x;
//    double y;
//    double z;
//
//    _Coordinates3D() {}
//    _Coordinates3D( const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}
//    _Coordinates3D( const _Coordinates3D& p) : x(p.x), y(p.y), z(p.z) {}
//    _Coordinates3D operator + (const double &d) const { return _Coordinates3D(d+x, d+y, d+z); }
//    _Coordinates3D operator + (const _Coordinates3D &p) const { return _Coordinates3D(p.x+x, p.y+y, p.z+z); }
//    _Coordinates3D operator * (const double &d) const { return _Coordinates3D(d*x, d*y, d*z); }
//    _Coordinates3D operator * (const _Coordinates3D &p) const { return _Coordinates3D(p.x*x, p.y*y, p.z*z); }
//} Coordinates3D;
//
//typedef struct _ColorRGB {
//	uint8_t r;
//	uint8_t g;
//	uint8_t b;
//} ColorRGB;
//
//typedef struct _ColorRGBA {
//        uint8_t r;
//        uint8_t g;
//        uint8_t b;
//        uint8_t a;
//
//    _ColorRGBA() {}
//    _ColorRGBA( const uint8_t& r, const uint8_t& g, const uint8_t& b, const uint8_t& a) : r(r), g(g), b(b), a(a) {}
//    _ColorRGBA( const _ColorRGBA& p) : r(p.r), g(p.g), b(p.b), a(p.a) {}
//    _ColorRGBA operator + (const _ColorRGBA &p) const { return _ColorRGBA(p.r+r, p.g+g, p.b+b, p.a+a); }
//    _ColorRGBA operator * (const double &d) const { return _ColorRGBA(round(d*(double)r),round(d*(double)g),round(d*(double)b),round(d*(double)a)); }
//} ColorRGBA;
//
//
//typedef struct _Coordinates4D {
//	double x;
//	double y;
//	double z;
//	double w;
//
//    _Coordinates4D() {}
//    _Coordinates4D( const double& x, const double& y, const double& z, const double& w) : x(x), y(y), z(z), w(w) {}
//    _Coordinates4D( const _Coordinates4D& p) : x(p.x), y(p.y), z(p.z), w(p.w) {}
//    _Coordinates4D operator + (const _Coordinates4D &p) const { return _Coordinates4D(p.x+x, p.y+y, p.z+z, p.w+w); }
//    _Coordinates4D operator * (const double &d) const { return _Coordinates4D(d*x, d*y, d*z, d*w); }
//} Coordinates4D;
//
//typedef struct _Mat3D {
//	double d[3][3];
//} Mat3D;
//
//typedef struct _Mat4D {
//	double d[4][4];
//} Mat4D;
//
//typedef struct _Polygon4D {
//	int numVertices;
//	Coordinates4D vertices[10];
//	Coordinates3D normals[10];
//	ColorRGBA colors[10];
//} Polygon4D;


enum LineClipStatus {
	LINE_CLIP_REJECT = 0x00,
	LINE_CLIP_ACCEPT = 0x01,
	LINE_CLIP_A = 0x02,
	LINE_CLIP_B = 0x04,
	LINE_CLIP_BOTH = 0x06
};

enum FrameBufferType {
	FBT_RGB,	// 3 bytes
	FBT_RGBA,	// 4 bytes
	FBT_DEPTH	// 4 bytes, doubles
};

typedef struct _FrameBuffer {
	FrameBufferType type;
	
	int rows;
	int cols;
	uint8_t* data;
    
    std::mutex mutex;
	
	void clear(void* color) {
		switch (type) {
			case FBT_RGB:
				std::fill_n((ColorRGB*)&data[0], cols*rows, *(ColorRGB*)color);
				break;
				
			case FBT_RGBA:
				std::fill_n((ColorRGBA*)&data[0], cols*rows, *(ColorRGBA*)color);
				break;
				
			case FBT_DEPTH:
				std::fill_n((double*)&data[0], cols*rows, *(double*)color);
				break;
		}
	}
	
	_FrameBuffer(): rows(0), cols(0), data(NULL) {};
	~_FrameBuffer() {
		if (data) {
			free(data);
			data = NULL;
		}
	};
} FrameBuffer;

void asTexImage2d(FrameBuffer* fbo, FrameBufferType type, int width, int height);

typedef struct _DepthBuffer {
	double* d;
	int width;
	int height;
	
	_DepthBuffer(): width(1), height(1) {
			d = (double*)malloc(sizeof(double));
		}
	
	void setSize(int width, int height) {
		if (width*height != this->width*this->height) {
			d = (double*)realloc(d, sizeof(double)*width*height);
		}
		this->width = width;
		this->height = height;
	}
	
	void reset() {
//		memset(d, 0, sizeof(double)*width*height);
//        std::fill_n(&d[0], width*height, 0);
		std::fill_n(&d[0], width*height, std::numeric_limits<double>::max());
	}
	
	~_DepthBuffer() {
		free(d);
	}
} DepthBuffer;

typedef struct _FragmentInfo {
	Coordinates2D pixel;
	Coordinates4D location3D;
	Coordinates3D normal;
	ColorRGBA color;
//	Mat4D modeView;
	void* data;
    void* interpolated;
	ColorRGBA* colorOutput;
} FragmentInfo;



// Can any software renderer not use this gem from Quake?
float Q_rsqrt( float number );

// Shader stuff
void defaultFragment(const FragmentInfo&);

// Depth buffer stuffs
void setFloatDotWithDepthBuffer( double x, double y, double depth, DepthBuffer* depthBuffer);
void setWithDepthBuffer( Coordinates2D pt, char c, double depth, DepthBuffer* depthBuffer);
void setWithDepthBuffer( Coordinates4D pt, char c, double depth, DepthBuffer* depthBuffer);
void setWithShader( Coordinates2D& pixel, double depth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&));
void lineWithDepthBuffer( Coordinates4D a3, Coordinates4D b3, DepthBuffer* depthBuffer);
void triangleWithDepthBuffer( Coordinates4D a3, Coordinates4D b3, Coordinates4D c3, DepthBuffer* depthBuffer, char fill);
void drawPolygon( Polygon4D& poly, DepthBuffer* depthBuffer, char fill, int &debugLine);
void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line);

// Construction:
Mat4D makeWindowTransform(int screenSizeX, int screenSizeY, double characterAspect);
Mat4D rotationFromAngleAndUnitAxis( double radians, Coordinates3D axis);
Mat4D translationMatrix( double x, double y, double z);
Mat4D scaleMatrix( double x, double y, double z );
Mat4D invert3x3Slow( const Mat4D& mat );
Mat4D transpose( const Mat4D& mat );
Mat4D projectionMatrixOrtho(double width, double height, double zfar, double znear);
Mat4D projectionMatrixPerspective(double fov, double aspect, double zfar, double znear);
void fillPolygonNormals(Polygon4D* polygons, int count);	// 

// Conversions:
Coordinates2D onlyXY(Coordinates3D& input);
Coordinates2D onlyXY(Coordinates4D& input);

// Manipulation:
Mat4D transpose( Mat4D& input );
Coordinates3D normalizeVector(const Coordinates3D& input);
Coordinates3D normalizeVector(Coordinates4D& input);	// only in 3D
Coordinates3D normalizeVectorFast(Coordinates3D& input);
Coordinates3D normalizeVectorFast(Coordinates4D& input);	// only in 3D

// Operations
int mod(int a, int b);
Mat3D matrixMultiply(Mat3D& a, Mat3D& b);
Mat4D matrixMultiply(Mat4D& a, Mat4D& b);
Coordinates3D matrixVectorMultiple(Mat3D& rotation, const Coordinates3D& vec);
Coordinates4D matrixVectorMultiply(Mat4D& rotation, const Coordinates4D& vec);
Coordinates3D matrixVectorMultiply(Mat4D& rotation, const Coordinates3D& vec); // operated in 3d
Coordinates3D crossProduct(Coordinates3D a, Coordinates3D b);
Coordinates4D crossProduct(Coordinates4D a, Coordinates4D b);
double dotProduct(const Coordinates4D& a, const Coordinates4D& b); // operates in 3d
double dotProduct(const Coordinates3D& a, const Coordinates3D& b); // operates in 3d
double dotProduct(const Coordinates4D& a, const Coordinates3D& b); // operates in 3d
double dotProduct(const Coordinates3D& a, const Coordinates4D& b); // operates in 3d
//double vectorSquared(Coordinates4D& a);	// operates in 3D
Coordinates4D vectorAdd(Coordinates4D a, Coordinates4D b);
Coordinates4D vectorSubtract(Coordinates4D a, Coordinates4D b);
Coordinates3D vectorAdd(Coordinates3D a, Coordinates3D b);
Coordinates3D vectorSubtract(Coordinates3D a, Coordinates3D b);
Coordinates3D vectorScale(Coordinates3D a, double scale);


Coordinates3D interpolate(Coordinates3D& a, Coordinates3D& b, double factor);
Coordinates4D interpolate(Coordinates4D& a, Coordinates4D& b, double factor);
ColorRGBA interpolate(ColorRGBA& a, ColorRGBA& b, double factor);
//Coordinates4D perspectiveInterpolate(Coordinates4D& a, Coordinates4D& b, double aDepth, double bDepth, double factor);
Coordinates4D perspectiveInterpolate(Coordinates4D& a, Coordinates4D& b, double factor);
Coordinates3D perspectiveInterpolate(Coordinates3D& a, Coordinates3D& b, double aDepth, double bDepth, double correctDepth, double factor);
Coordinates4D perspectiveInterpolateInv(Coordinates4D& a, Coordinates4D& b, double aInvDepth, double bInvDepth, double correctDepth, double factor);
Coordinates3D perspectiveInterpolateInv(Coordinates3D& a, Coordinates3D& b, double aInvDepth, double bInvDepth, double correctDepth, double factor);

// Rendering
Coordinates3D clipRGB(Coordinates3D rgb);
Coordinates3D rgbToHsv( const Coordinates3D& rgb);
Coordinates3D hsvToRgb(const Coordinates3D& hsvIn);
Coordinates3D hslToRgb(Coordinates3D hsl);
void setRGB( const Coordinates2D& pt, const Coordinates3D& rgb);
void setFrameBufferRGBA(int x, int y, FrameBuffer* fbo, const ColorRGBA& value);
void rasterize(Coordinates4D* vertices, int edgeIndices[][2], int numEdges, Mat4D& windowTransform, DepthBuffer* depthBuffer);
//void triangleFill(Coordinates4D& p1, Coordinates4D& p2, Coordinates4D& p3, FrameBuffer* output);
//void triangleFill(Coordinates4D* vertices, ColorRGBA* colors, FrameBuffer* fbo);
bool compareCoordinatesFrag(FragmentInfo i1, FragmentInfo i2);
void triangleFill(FragmentInfo* fragments, FrameBuffer* fbo);
void drawPolygonWithTriangles( Polygon4D& poly, FrameBuffer* fbo);
void rasterizeTriangle(Coordinates4D* vertices, int edgeIndices[][3], int numEdges, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill=' ', int line=0);
void rasterizeQuads(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill, int &line);
void rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line);

void rasterizePolygon(Polygon4D* polygons, int count, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill, int &line);
void rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line);

char getp( Coordinates4D* pts, double err);
void drawHorizonalLineWithDepthBuffer(int x1, int x2, int y, char c, double depth1, double depth2, DepthBuffer* depthBuffer);
void drawHorizonalLineWithShader(int x1, int x2, int y, double depth1, double depth2, Coordinates4D& point1, Coordinates4D& point2,  Coordinates3D& normal1, Coordinates3D& normal2, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&));
void drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&));

// Specific algorithms:
LineClipStatus liangBarskyHomogeneous(Coordinates4D &point1, Coordinates4D &point2);
//int clipTriangle(Coordinates4D input[3], Coordinates4D* output, int& debugLine);
Polygon4D clipTriangle(Coordinates4D input[3], int& line);
int clipPolygon(Polygon4D& input, Polygon4D* output, int& debugLine);

#endif
