#ifndef CURSES_GFX_TYPES_H
#define CURSES_GFX_TYPES_H

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <mutex>

#include "curses-gfx.h"


typedef struct _Coordinates2Df {
    double x;
    double y;
    
    _Coordinates2Df() {}
    _Coordinates2Df( const double& x, const double& y) : x(x), y(y) {}
    _Coordinates2Df( const _Coordinates2Df& p) : x(p.x), y(p.y) {}
    _Coordinates2Df operator + (const _Coordinates2Df &p) const { return _Coordinates2Df(p.x+x, p.y+y); }
    _Coordinates2Df operator * (const double &d) const { return _Coordinates2Df(d*x, d*y); }
} Coordinates2Df;

typedef struct _Coordinates3D {
    double x;
    double y;
    double z;
    
    _Coordinates3D() {}
    _Coordinates3D( const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}
    _Coordinates3D( const _Coordinates3D& p) : x(p.x), y(p.y), z(p.z) {}
    _Coordinates3D operator + (const double &d) const { return _Coordinates3D(d+x, d+y, d+z); }
    _Coordinates3D operator + (const _Coordinates3D &p) const { return _Coordinates3D(p.x+x, p.y+y, p.z+z); }
    _Coordinates3D operator * (const double &d) const { return _Coordinates3D(d*x, d*y, d*z); }
    _Coordinates3D operator * (const _Coordinates3D &p) const { return _Coordinates3D(p.x*x, p.y*y, p.z*z); }
} Coordinates3D;

typedef struct _ColorRGB {
	uint8_t r;
	uint8_t g;
	uint8_t b;
} ColorRGB;

typedef struct _ColorRGBA {
        uint8_t r;
        uint8_t g;
        uint8_t b;
        uint8_t a;
    
    _ColorRGBA() {}
    _ColorRGBA( const uint8_t& r, const uint8_t& g, const uint8_t& b, const uint8_t& a) : r(r), g(g), b(b), a(a) {}
    _ColorRGBA( const _ColorRGBA& p) : r(p.r), g(p.g), b(p.b), a(p.a) {}
    _ColorRGBA operator + (const _ColorRGBA &p) const { return _ColorRGBA(p.r+r, p.g+g, p.b+b, p.a+a); }
    _ColorRGBA operator * (const double &d) const { return _ColorRGBA(round(d*(double)r),round(d*(double)g),round(d*(double)b),round(d*(double)a)); }
} ColorRGBA;


typedef struct _Coordinates4D {
	double x;
	double y;
	double z;
	double w;
    
    _Coordinates4D() {}
    _Coordinates4D( const double& x, const double& y, const double& z, const double& w) : x(x), y(y), z(z), w(w) {}
    _Coordinates4D( const _Coordinates4D& p) : x(p.x), y(p.y), z(p.z), w(p.w) {}
    _Coordinates4D operator + (const _Coordinates4D &p) const { return _Coordinates4D(p.x+x, p.y+y, p.z+z, p.w+w); }
    _Coordinates4D operator * (const double &d) const { return _Coordinates4D(d*x, d*y, d*z, d*w); }
} Coordinates4D;

typedef struct _Mat3D {
	double d[3][3];
} Mat3D;

typedef struct _Mat4D {
	double d[4][4];
} Mat4D;

typedef struct _Polygon4D {
	int numVertices;
	Coordinates4D vertices[10];
	Coordinates3D normals[10];
	ColorRGBA colors[10];
} Polygon4D;



#endif
