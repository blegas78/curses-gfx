#ifndef CURSES_GFX_TYPES_H
#define CURSES_GFX_TYPES_H

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <mutex>

//#include "curses-gfx.h"
/* TODO: I'd like to make the folowing template but there are some architectural challenges to implement it and be ultra fast:
template <typename T, unsigned int length> class GfxVector {
    T d[length];
    GfxVector() {}
//    _Coordinates3D( const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}
    GfxVector( const GfxVector<T,length>& p) { std::memcpy(d, p.d, length*sizeof(T)); }
    GfxVector operator + (const double &d) const { return GfxVector<T,length>(d+x, d+y, d+z); }
    GfxVector operator + (const GfxVector<T,length> &p) const { return _Coordinates3D(p.x+x, p.y+y, p.z+z); }
    GfxVector operator * (const double &d) const { return GfxVector<T,length>(d*x, d*y, d*z); }
    GfxVector operator * (const GfxVector<T,length> &p) const { return GfxVector<T,length>(p.x*x, p.y*y, p.z*z); }
};
*/

typedef struct _Coordinates2D {
    int x;
    int y;
    
    _Coordinates2D() {}
    _Coordinates2D( const int& x, const int& y) : x(x), y(y) {}
    _Coordinates2D( const _Coordinates2D& p) : x(p.x), y(p.y) {}
    _Coordinates2D operator + (const _Coordinates2D &p) const { return _Coordinates2D(p.x+x, p.y+y); }
    _Coordinates2D operator * (const double &d) const { return _Coordinates2D(d*(double)x, d*(double)y); }
} Coordinates2D;

typedef struct _Coordinates2Df {
    double x;
    double y;
    
    _Coordinates2Df() {}
    _Coordinates2Df( const double& x, const double& y) : x(x), y(y) {}
    _Coordinates2Df( const _Coordinates2Df& p) : x(p.x), y(p.y) {}
    _Coordinates2Df operator + (const _Coordinates2Df &p) const { return _Coordinates2Df(p.x+x, p.y+y); }
    _Coordinates2Df operator * (const double &d) const { return _Coordinates2Df(d*x, d*y); }
} Coordinates2Df;

class Coordinates4D;

typedef struct Coordinates3D {
    double x;
    double y;
    double z;
    
    Coordinates3D() {}
    Coordinates3D( const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}
    Coordinates3D( const Coordinates3D& p) : x(p.x), y(p.y), z(p.z) {}
    Coordinates3D operator + (const double &d) const { return Coordinates3D(d+x, d+y, d+z); }
    Coordinates3D operator + (const Coordinates3D &p) const { return Coordinates3D(p.x+x, p.y+y, p.z+z); }
    Coordinates3D operator * (const double &d) const { return Coordinates3D(d*x, d*y, d*z); }
    Coordinates3D operator * (const Coordinates3D &p) const { return Coordinates3D(p.x*x, p.y*y, p.z*z); }
    
    
    Coordinates3D( const Coordinates4D& p);
} Coordinates3D;

class Coordinates4D {
public:
    double x;
    double y;
    double z;
    double w;
    
    Coordinates4D() {}
    Coordinates4D( const double& x, const double& y, const double& z, const double& w) : x(x), y(y), z(z), w(w) {}
    Coordinates4D( const Coordinates4D& p) : x(p.x), y(p.y), z(p.z), w(p.w) {}
    Coordinates4D operator + (const Coordinates4D &p) const { return Coordinates4D(p.x+x, p.y+y, p.z+z, p.w+w); }
    Coordinates4D operator * (const double &d) const { return Coordinates4D(d*x, d*y, d*z, d*w); }
//    _Coordinates4D operator = (const Coordinates3D &d) const {return _Coordinates4D(d.x, d.y, d.z, 1);}
    
    Coordinates4D& operator += (const Coordinates4D &p) { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
    
    Coordinates4D( const Coordinates3D& p, const double& w);
};

typedef struct _ColorRGB {
	uint8_t r;
	uint8_t g;
	uint8_t b;
} ColorRGB;

double clamp(const double& input, const double& min, const double& max);

typedef struct _ColorRGBA {
        uint8_t r;
        uint8_t g;
        uint8_t b;
        uint8_t a;
    
    _ColorRGBA() {}
    _ColorRGBA( const uint8_t& r, const uint8_t& g, const uint8_t& b, const uint8_t& a) : r(r), g(g), b(b), a(a) {}
    _ColorRGBA( const _ColorRGBA& p) : r(p.r), g(p.g), b(p.b), a(p.a) {}
    _ColorRGBA operator + (const _ColorRGBA &p) const { return _ColorRGBA(p.r+r, p.g+g, p.b+b, p.a+a); }
    _ColorRGBA operator * (const double &d) const { return _ColorRGBA(clamp(round(d*(double)r),0,255),clamp(round(d*(double)g),0,255),clamp(round(d*(double)b),0,255),clamp(round(d*(double)a),0,255)); }
//    _ColorRGBA& operator = (const _ColorRGBA &d) {this->r = d.r; this->g = d.g; this->b = d.b; this->a = d.a; return *this;}
//    _ColorRGBA& operator = (const Coordinates4D &d) {this->r = d.x; this->g = d.y; this->b = d.z; this->a = d.w; return *this;}
//    _ColorRGBA& operator = (const Coordinates3D &d) {this->r = d.x; this->g = d.y; this->b = 127; this->a = 0; return *this;}
    
    
    _ColorRGBA(const Coordinates4D &p) : r(clamp(p.x*255,0,255)), g(clamp(p.y*255,0,255)), b(clamp(p.z*255,0,255)), a(clamp(p.w*255,0,255)) {}

} ColorRGBA;




typedef struct _Mat3D {
	double d[3][3];
} Mat3D;

typedef struct _Mat4D {
	double d[4][4];
    
    _Mat4D() {}
//    _Coordinates4D( const double* p[4]) { for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) d[i][j] = p[i][j]; }
//    _Mat4D operator * (const double &d) const { return _Mat4D(d*x, d*y, d*z, d*w); }
} Mat4D;

typedef struct _Polygon4D {
	int numVertices;
	Coordinates4D vertices[10];
	Coordinates3D normals[10];
	ColorRGBA colors[10];
} Polygon4D;



#endif
