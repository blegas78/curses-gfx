#include "curses-gfx-3d.h"
#include "curses-gfx.h"

#include <cmath>
#include <cassert>
#include <unistd.h>
#include <stdint.h>

#include <ncurses.h>
#include <vector>
#include <algorithm>

float Q_rsqrt( float number )
{
//	long i;
//	float x2, y;
//	const float threehalfs = 1.5F;
//
//	x2 = number * 0.5F;
//	y  = number;
////	i  = * ( long * ) &y;                       // evil floating point bit level hacking
//	memcpy(&i, &y, sizeof(i));                       // evil floating point bit level hacking
////	i  = 0x5f3759df - ( (i >> 1) & 0x7fffffffffffffff );               // what the fuck?
//	i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
//	y  = * ( float * ) &i;
//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
////	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
//
//	return y;
	
	// Per wikipedia, there are aliasing issues in modern C with the above casting so using unions seems to work:
	union {
			float    f;
			uint32_t i;
		} conv = { .f = number };
		conv.i  = 0x5f3759df - (conv.i >> 1);
		conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
		return conv.f;
}

Coordinates2D onlyXY(Coordinates3D& input) {
	Coordinates2D result;
	result.x = input.x + 0.5;
	result.y = input.y + 0.5;
	return result;
}

Coordinates2D onlyXY(Coordinates4D& input) {
	Coordinates2D result;
	result.x = input.x + 0.5;
	result.y = input.y + 0.5;
	return result;
}


Mat4D makeWindowTransform(int screenSizeX, int screenSizeY, double characterAspect) {
	
//	int screenSizeX, screenSizeY;
//	getmaxyx(stdscr, screenSizeY, screenSizeX);
	
	double screenAspect = (double)screenSizeX/(double)screenSizeY / characterAspect;
	
	// Window
//	Mat4D windowScale = scaleMatrix((double)(screenSizeX-1)/2, (double)(screenSizeY-1)/2, 1);	// Full accurate
	//	Mat4D translationScreen = translationMatrix((double)screenSizeX/2 - 0.5, (double)screenSizeY/2 - 0.5, 0);
//	Mat4D windowScale = scaleMatrix((double)(screenSizeX+2)/2, -(double)(screenSizeY+2)/2, 1);	// note the negative in Y, since y is positive down in a terminal
    Mat4D windowScale = scaleMatrix((double)(screenSizeX-2)/2, -(double)(screenSizeY-2)/2, 1);    // note the negative in Y, since y is positive down in a terminal
	Mat4D translationScreen = translationMatrix((double)screenSizeX/2, (double)screenSizeY/2, 0);
	return matrixMultiply(translationScreen, windowScale);
}

void setFloatDotWithDepthBuffer( double x, double y, double depth, DepthBuffer* depthBuffer) {
	if (depthBuffer != NULL) {
		if ( (int)x < 0 || (int)y < 0 || depthBuffer->width <= (int)x || depthBuffer->height <= (int)y) {
			return;
		}
		if (depth < depthBuffer->d[(int)x + depthBuffer->width*(int)y]) {
//			set(pt, c);
			drawDotFloat(x, y);
			depthBuffer->d[(int)x + depthBuffer->width*(int)y] = depth;
		}
	} else {
//		set(pt, c);
		drawDotFloat(x, y);
	}
}

void setWithDepthBuffer( Coordinates2D pt, char c, double depth, DepthBuffer* depthBuffer) {
	if (depthBuffer != NULL) {
		if ( pt.x < 0 || pt.y < 0 || depthBuffer->width <= pt.x || depthBuffer->height <= pt.y) {
			return;
		}
		if (depth < depthBuffer->d[pt.x + depthBuffer->width*pt.y]) {
			set(pt, c);
			depthBuffer->d[pt.x + depthBuffer->width*pt.y] = depth;
			return;
			
			
			if (depth > 0.4) {
				set(pt, '`');
			} else if (depth > 0.35) {
				set(pt, '.');
			} else if (depth > 0.3) {
				set(pt, ',');
			} else if (depth > 0.25) {
				set(pt, '=');
			} else if (depth > 0.2) {
				set(pt, '+');
			} else if (depth > 0.15) {
				set(pt, '*');
			} else if (depth > 0.10) {
				set(pt, '#');
			} else if (depth > 0.05) {
				set(pt, '%');
			} else {
				set(pt, '@');
				
			}
		}
	} else {
		set(pt, c);
	}
}

void defaultFragment(const FragmentInfo& fInfo) {
//	set(fInfo.pixel, '.');
#ifdef ASCII_EXTENDED
	//mvprintw(fInfo.pixel.y, fInfo.pixel.x, "░"); // light "▒" ▓ █
#else
	//set(fInfo.pixel, '+');
#endif
	fInfo.colorOutput->r = 255;
	fInfo.colorOutput->g = 255;
	fInfo.colorOutput->b = 0;
	fInfo.colorOutput->a = '+';
}

void setWithDepthBuffer( Coordinates4D pt, char c, double depth, DepthBuffer* depthBuffer) {
	setWithDepthBuffer(onlyXY(pt), c, depth, depthBuffer);
}

void setWithShader( Coordinates2D& pixel, double depth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&)) {
	if ( pixel.x < 0 || pixel.y < 0 || depthBuffer->width <= pixel.x || depthBuffer->height <= pixel.y) {
		return;
	}
	int depthIndex = pixel.x + depthBuffer->width*pixel.y;
	if (depth < depthBuffer->d[depthIndex]) {
		depthBuffer->d[depthIndex] = depth;
		
		FragmentInfo fInfo;
		fInfo.pixel = pixel;
		fInfo.location3D = pt3D;
		fInfo.normal = normal;
        fInfo.data = userData;
//        ColorRGBA colorOutput;
//        fInfo.colorOutput = &colorOutput;
		fragmentShader(fInfo);
//        Coordinates3D color = {(double)colorOutput.r, (double)colorOutput.g, (double)colorOutput.b};
        
        //setRGB(pixel, color);
		return;
//		const double maxDepth = 10.0;
//		const double minDepth = 3.0;
//		depth = -pt3D.z;
//		if (depth > maxDepth) {
//			set(pixel, ' ');
//		} else if (depth > maxDepth*0.9) {
//			set(pixel, '`');
//		} else if (depth > maxDepth*0.8) {
//			set(pixel, '.');
//		} else if (depth > maxDepth*0.7) {
//			set(pixel, ',');
//		} else if (depth > maxDepth*0.6) {
//			set(pixel, '=');
//		} else 250if (depth > maxDepth*0.5) {
//			set(pixel, '+');
//		} else if (depth > maxDepth*0.4) {
//			set(pixel, '*');
//		} else if (depth > maxDepth*0.3) {
//			set(pixel, '#');
//		} else if (depth > maxDepth*0.2) {
//			set(pixel, '%');
//		} else {
//			set(pixel, '@');
//		}
	}
}

// This is basically https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
// The silly thing here is that for determining individual characters,
//  floating operations are done which circumvents the integer-only intention
//  in the original algorithm
// This was super helpful: https://zingl.github.io/Bresenham.pdf
void lineWithDepthBuffer( Coordinates4D a3, Coordinates4D b3, DepthBuffer* depthBuffer) {
	Coordinates2D a = onlyXY(a3);
	Coordinates2D b = onlyXY(b3);
	setWithDepthBuffer( a, 'o', a3.z, depthBuffer);
	setWithDepthBuffer( b, 'o', b3.z, depthBuffer);
	//	set( a, 'O');
	//	set( b, 'O');
	
	int dx = abs(b.x - a.x);
	int sx = a.x < b.x ? 1 : -1;
	int dy = -abs(b.y - a.y);
	int sy = a.y < b.y ? 1 : -1;
	int err = dx + dy;
	int e2;
	
	int errs[3];
	double errsNormalized[3];
	Coordinates2D pts[3];
	
	double alpha;
	double dxn, dyn;
	double depths[3];
	
	// this is a pretty shoehorned, inefficient method of depth tracking
	//	double depth = a3.y;
	double lineMagnitudeSq = sqrt(dx*dx + dy*dy);
	
	int count = 0;
	while (true)
	{
		//		double depth = a3.z - alpha*(a3.z - b3.z);
		if (count++ > 2) {
			setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
			
		}
		//		set( pts[1], getp( pts, errsNormalized[1] + 0.5));
		
		pts[0] = pts[1];
		pts[1] = pts[2];
		pts[2] = a;
		
		errs[0] = errs[1];
		errs[1] = errs[2];
		errs[2] = err;
		
		errsNormalized[0] = errsNormalized[1];
		errsNormalized[1] = errsNormalized[2];
		errsNormalized[2] = (double)(errs[1])/(dx - dy);
		if (sy == -1 ) {
			errsNormalized[2] = 0.0 - errsNormalized[2];
		}
		
		dxn = a.x-b.x;
		dyn = a.y-b.y;
		alpha = ((double)sqrt(dxn*dxn + dyn*dyn))/lineMagnitudeSq;
		depths[0] = depths[1];
		depths[1] = depths[2];
		depths[2] = 1.0/(1.0/b3.z - alpha*(1.0/b3.z - 1.0/a3.z));
		
		if (a.x == b.x && a.y == b.y) {
			break;
		}
		
		e2 = 2*err;
		if (e2 >= dy) {
			err += dy;
			a.x += sx;
		} else {
			errs[2] -= dy;
		}
		if (e2 <= dx) {
			err += dx;
			a.y += sy;
		} else {
			errs[2] -= dx;
		}
	}
	
	// add the final point
	setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
	//	set( pts[1], getp( pts, errsNormalized[1] + 0.5));
	
}

char getp( Coordinates4D* pts, double err)
{
	if (abs(pts[0].y - pts[2].y) < 2)	// checks the slope
	{
		// 28 pixel charaxter height
		// _ -> 23
		// - -> 16 -> separation for above = (23-16)/2/28
		// " -> 9
		// ` -> 8
		
		if (err > 1.0-((9.0+8.0)/(2.0*28.0))) { //  _-"`
			return '`';
		}
		if (err > 1.0-((16.0+9.0)/(2.0*28.0))) { //  _-"`
			return '"';
		}
		if (err > 1.0-((23.0+16.0)/(2.0*28.0)))
		{
			return '-';
		}
		return '_';
	}
	
	if (abs(pts[0].x - pts[2].x) < 2 &&
		(pts[0].x >= pts[2].x || pts[1].x != pts[2].x) &&
		(pts[0].x <= pts[2].x || pts[1].x != pts[0].x))
	{
		return '|';
	}
	
	int mX = pts[0].y < pts[2].y ? pts[0].x : pts[2].x;
	return mX < pts[1].x ? '\\' : '/';
}

void drawHorizonalLineWithDepthBuffer(int x1, int x2, int y, char c, double depth1, double depth2, DepthBuffer* depthBuffer) {
	Coordinates2D point;
	int diff = x2-x1;
	if (abs(diff) < 2) {
		return;
	}
	int increment = x2 > x1 ? 1 : -1;
	double depth;
	point.y = y;
	for (int i = x1+increment; i != x2; i += increment) {
		point.x = i;
		
		depth = depth1 - ((double)(i - x1)/(double)diff) * (depth1 - depth2);
		
		setWithDepthBuffer(point, c, depth, depthBuffer);
	}
	
}

void drawHorizonalLineWithShader(int x1, int x2, int y, double depth1, double depth2, Coordinates4D& point1, Coordinates4D& point2,  Coordinates3D& normal1, Coordinates3D& normal2, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&)) {
	Coordinates2D pixel;
	int diff = x2-x1;
	if (abs(diff) < 2) {
		return;
	}
	int increment = x2 > x1 ? 1 : -1;
	double factor;
	
	Coordinates4D point;
	Coordinates3D normal;
	double depth;
	
	pixel.y = y;
	for (int i = x1+increment; i != x2; i += increment) {
		pixel.x = i;
		
		factor = ((double)(i - x1)/(double)diff);
		
//		point = vectorSubtract(point2, point1);
//		point.x *= factor;
//		point.y *= factor;
//		point.z *= factor;
//		point.w *= factor;
//		point = vectorAdd(point, point1);
//		point = interpolate(point1, point2, factor);
//		point.z = 1.0/(1.0/point1.z + factor*(1.0/point2.z - 1.0/point1.z));
		
//		point = perspectiveInterpolate(point1, point2, 0, 0, 0, factor);
		point = perspectiveInterpolate(point1, point2, factor);
		
		depth = 1.0/(1.0/depth1 + factor*(1.0/depth2 - 1.0/depth1));
		
//		normal = vectorSubtract(normal2, normal1);
//		normal.x *= factor;
//		normal.y *= factor;
//		normal.z *= factor;
//		normal = vectorAdd(normal, normal1);
		normal = interpolate(normal1, normal2, factor);
		normal = normalizeVector(normal);
		
		setWithShader(pixel, depth, point, normal, userData, depthBuffer, fragmentShader);
		
	}
	
}

void triangleWithDepthBuffer( Coordinates4D a3, Coordinates4D b3, Coordinates4D c3, DepthBuffer* depthBuffer, char fill) {
		
//	int line = 0;
//	char string[100];
	
	Coordinates4D ac = a3;
	Coordinates4D bc = b3;
	Coordinates4D cc = c3;
	
	// Find highest vertex
	if (a3.y > b3.y) {
		if (a3.y > c3.y) {
			ac = a3;
			if (b3.y > c3.y) {
				bc = b3;
				cc = c3;
			} else {
				bc = c3;
				cc = b3;
			}
		} else {
			ac = c3;
			bc = a3;
			cc = b3;
		}
	} else {
		if (b3.y > c3.y) {
			ac = b3;
			if (a3.y > c3.y) {
				bc = a3;
				cc = c3;
			} else {
				bc = c3;
				cc = a3;
			}
		} else {
			ac = c3;
			bc = b3;
			cc = a3;
		}
	}
	
	setWithDepthBuffer( ac, 'o', ac.z, depthBuffer);
	setWithDepthBuffer( bc, 'o', bc.z, depthBuffer);
	setWithDepthBuffer( cc, 'o', cc.z, depthBuffer);
	
	Coordinates2D a = onlyXY(ac);
	Coordinates2D b = onlyXY(bc);
	Coordinates2D c = onlyXY(cc);
	
//	mvprintw(line++, 0, "A: %d %d", a.x, a.y);
//	mvprintw(line++, 0, "B: %d %d", b.x, b.y);
//	mvprintw(line++, 0, "C: %d %d", c.x, c.y);
	
	int dx = abs(b.x - a.x);
	int sx = a.x < b.x ? 1 : -1;
	int dy = -abs(b.y - a.y);
	int sy = a.y < b.y ? 1 : -1;
	int err = dx + dy;
	int e2;
	
	Coordinates2D a2 = a;
	int dx2 = abs(c.x - a.x);
	int sx2 = a.x < c.x ? 1 : -1;
	int dy2 = -abs(c.y - a.y);
	int sy2 = a.y < c.y ? 1 : -1;
	int err2 = dx2 + dy2;
	int e22;
	
	int errs[3];
	double errsNormalized[3];
	Coordinates2D pts[3];
	
	double alpha;
	double dxn, dyn;
	double depths[3];
	
	int errs2[3];
	double errsNormalized2[3];
	Coordinates2D pts2[3];
	
	double alpha2;
	double dxn2, dyn2;
	double depths2[3];
	
	// this is a pretty shoehorned, inefficient method of depth tracking
	//	double depth = a3.y;
	double lineMagnitudeSq = sqrt(dx*dx + dy*dy);
	double lineMagnitudeSq2 = sqrt(dx2*dx2 + dy2*dy2);
	
	
	bool skipMainLine = false;
	bool skipOtherLine = true;
	int lineStartX[3];// = a.x - sx;
	int lineEndX[3];
	int lineY[3];// = a.y;
	double lineDepthStart[3];
	double lineDepthEnd[3];
	
	
	bool linePrimed = false;
	lineStartX[0] = a.x;// + (a2.x < a.x ? -1 : 1);
	lineEndX[0] = lineStartX[0];
	lineY[0] = a.y;
	lineStartX[1] = lineStartX[0];
	lineEndX[1] = lineStartX[0];
	lineY[1] = a.y;
	lineStartX[2] = lineStartX[0];
	lineEndX[2] = lineStartX[0];
	lineY[2] = a.y;
	//	err2 = errs2[2];
	int t = 0;
	int t2 = 0;
	static double depth = -10;
	depth *= -1;

	
	bool ignoreFirstLine = true;
	int lineCount = 0;
	while(true)
	{
//				refresh();
//				usleep(1000000);
		//		double depth = a3.z - alpha*(a3.z - b3.z);
		
		//		set( pts[1], getp( pts, errsNormalized[1] + 0.5));
		if (!skipMainLine) {
			
			if (t++ > 2) {
				setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
				//				setWithDepthBuffer( pts[2], '.', depths[1], depthBuffer);
			}
			pts[0] = pts[1];
			pts[1] = pts[2];
			pts[2] = a;
			
			errs[0] = errs[1];
			errs[1] = errs[2];
			errs[2] = err;
			
			errsNormalized[0] = errsNormalized[1];
			errsNormalized[1] = errsNormalized[2];
			errsNormalized[2] = (double)(errs[1])/(dx - dy);
			if (sy == -1 ) {
				errsNormalized[2] = 0.0 - errsNormalized[2];
			}
			
			dxn = a.x-b.x;
			dyn = a.y-b.y;
			alpha = ((double)sqrt(dxn*dxn + dyn*dyn))/lineMagnitudeSq;
			depths[0] = depths[1];
			depths[1] = depths[2];
//			depths[2] = bc.z - alpha*(bc.z - ac.z);
			depths[2] = 1.0/(1.0/bc.z - alpha*(1.0/bc.z - 1.0/ac.z));
			lineDepthStart[2] = depths[2];
		}
		
		
		if (a.x == b.x && a.y == b.y) {
			break;
		}
		
		if (!skipMainLine) {
			
			
			e2 = 2*err;
			if (e2 >= dy) {
				err += dy;
				a.x += sx;
			} else {
				errs[2] -= dy;
			}
			if (e2 <= dx) {
				err += dx;
				a.y += sy;
				lineStartX[0] = lineStartX[1];
				lineStartX[1] = lineStartX[2];
				lineY[0] = lineY[1];
				lineY[1] = lineY[2];
				lineEndX[0] = lineEndX[1];
				lineEndX[1] = lineEndX[2];
				
				lineDepthStart[0] = lineDepthStart[1];
				lineDepthStart[1] = lineDepthStart[2];
				lineDepthEnd[0] = lineDepthEnd[1];
				lineDepthEnd[1] = lineDepthEnd[2];
//				if (sx > 0 && sx2 > 0) {
//					lineStartX[2] = a.x + (a2.x < a.x ? 0 : 1);
//				} else {
				lineStartX[2] = a.x;// + (a2.x < a.x ? -1 : 1);
					
//				}
					lineY[2] = a.y;
				if (lineCount++ > 1) {
					
//					mvprintw(line++, 0, "line from %d-%d @ row %d", lineStartX[0], lineEndX[0], lineY[0]);
					drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
				}
				
				
				skipOtherLine = false;
				
				//				sprintf(string, "AB incremented, lineStartX[1] = %d", lineStartX[1]);
				//				mvaddstr(line++, 0, string);
			} else {
				if (abs(a2.x - a.x) < abs(a2.x - lineStartX[2])) {
//					mvaddstr(line++, 0, string);
					lineStartX[2] = a.x;// + (a2.x < a.x ? -1 : 1);
//					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
				}
				//				}
				errs[2] -= dx;
				skipOtherLine = true;
			}
		}
		
		if (!skipOtherLine) {
			if (t2++ > 2) {
				setWithDepthBuffer( pts2[1], getp( pts2, errsNormalized2[1] + 0.5), depths2[1], depthBuffer);
				//				setWithDepthBuffer( pts2[2], '`', depths2[1], depthBuffer);
			}
			pts2[0] = pts2[1];
			pts2[1] = pts2[2];
			pts2[2] = a2;
			errs2[0] = errs2[1];
			errs2[1] = errs2[2];
			errs2[2] = err2;
			
			errsNormalized2[0] = errsNormalized2[1];
			errsNormalized2[1] = errsNormalized2[2];
			errsNormalized2[2] = (double)(errs2[1])/(dx2 - dy2);
			if (sy2 == -1 ) {
				errsNormalized2[2] = 0.0 - errsNormalized2[2];
			}
			dxn2 = a2.x-c.x;
			dyn2 = a2.y-c.y;
			alpha2 = ((double)sqrt(dxn2*dxn2 + dyn2*dyn2))/lineMagnitudeSq2;
			depths2[0] = depths2[1];
			depths2[1] = depths2[2];
//			depths2[2] = cc.z - alpha2*(cc.z - ac.z);
			depths2[2] = 1.0/(1.0/cc.z - alpha2*(1.0/cc.z - 1.0/ac.z));
			lineDepthEnd[2] = depths2[2];
		}
		
		
		if (!skipOtherLine) {
			e22 = 2*err2;
			if (e22 >= dy2) {
				err2 += dy2;
				
				a2.x += sx2;
			} else {
				errs2[2] -= dy2;
			}
			
			if (e22 <= dx2) {
				err2 += dx2;
				a2.y += sy2;
				//				lineEndX[0] = lineEndX[1];
				lineEndX[2] = a2.x;// + (a2.x < lineStartX[1] ? -1 : 0);
				skipMainLine = false;
				
				if (a.x == b.x && a.y == b.y) {
					skipMainLine = true;
					skipOtherLine = false;
				}
				
				
				
			} else {
				errs2[2] -= dx2;
				skipMainLine = true;
				
				if (!(a.x == b.x && a.y == b.y)) {
					if (abs(a2.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
						lineEndX[1] = a2.x;
					}
				}
			}
		}
		
		
	}
	
	// add the final points
	setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
	//setWithDepthBuffer( pts2[1], getp( pts2, errsNormalized2[1] + 0.5), depths2[1], depthBuffer);
	
	pts[0] = pts[1];
	pts[1] = pts[2];
	pts[2] = b;
	depths[1] = depths[2];
	errsNormalized[1] = errsNormalized[2];
	setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
	//	setWithDepthBuffer( pts2[2], getp( pts2, errsNormalized2[2] + 0.5), depths2[2], depthBuffer);
	
	
//	lineEndX[2] = lineEndX[1];
//	lineEndX[2] = a2.x + (a2.x > a.x ? 100 : -100);
	
	lineStartX[0] = lineStartX[1];
	lineStartX[1] = lineStartX[2];
	
	lineY[0] = lineY[1];
	lineY[1] = lineY[2];
	lineEndX[0] = lineEndX[1];
	lineEndX[1] = lineEndX[2];
	
	lineDepthStart[0] = lineDepthStart[1];
	lineDepthStart[1] = lineDepthStart[2];
	lineDepthEnd[0] = lineDepthEnd[1];
	lineDepthEnd[1] = lineDepthEnd[2];

	if (lineY[0] != (int)ac.y) {
//		mvprintw(line++, 0, "printing midline because %d != %d", lineY[0], (int)ac.y);
		drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
	}
	
	double oldAcZ = ac.z;
	ac = bc;
	bc = cc;
	
	
	a = b;
	b = c;
	dx = abs(b.x - a.x);
	sx = a.x < b.x ? 1 : -1;
	dy = -abs(b.y - a.y);
	sy = a.y < b.y ? 1 : -1;
	err = dx + dy;
	//	e2;
	t = 0;
	skipOtherLine = true;
	skipMainLine = false;
	
	lineMagnitudeSq = sqrt(dx*dx + dy*dy);
	
	ignoreFirstLine = true;
	
	
	while(true)
	{
		if (a.y != a2.y) {
			skipMainLine = true;
			skipOtherLine = false;
		}
		
		//		set( pts[1], getp( pts, errsNormalized[1] + 0.5));
		if (!skipMainLine) {
			
			if (t++ > 2) {
				skipOtherLine = false;
				setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
			}
			pts[0] = pts[1];
			pts[1] = pts[2];
			pts[2] = a;
			
			errs[0] = errs[1];
			errs[1] = errs[2];
			errs[2] = err;
			
			errsNormalized[0] = errsNormalized[1];
			errsNormalized[1] = errsNormalized[2];
			errsNormalized[2] = (double)(errs[1])/(dx - dy);
			if (sy == -1 ) {
				errsNormalized[2] = 0.0 - errsNormalized[2];
			}
			
			dxn = a.x-b.x;
			dyn = a.y-b.y;
			alpha = ((double)sqrt(dxn*dxn + dyn*dyn))/lineMagnitudeSq;
			depths[0] = depths[1];
			depths[1] = depths[2];
//			depths[2] = bc.z - alpha*(bc.z - ac.z);
			depths[2] = 1.0/(1.0/bc.z - alpha*(1.0/bc.z - 1.0/ac.z));
			lineDepthStart[2] = depths[2];
		}
		
		
		if (a.x == b.x && a.y == b.y) {
			//			break;
			skipMainLine = true;
			skipOtherLine = false;
		}
		
		
		
		if (!skipMainLine) {
			e2 = 2*err;
			if (e2 >= dy) {
				err += dy;
				a.x += sx;
			} else {
				errs[2] -= dy;
			}
			if (e2 <= dx) {
				err += dx;
				a.y += sy;
				lineStartX[0] = lineStartX[1];
				lineStartX[1] = lineStartX[2];
				lineY[0] = lineY[1];
				lineY[1] = lineY[2];
				lineEndX[0] = lineEndX[1];
				lineEndX[1] = lineEndX[2];
				
				lineDepthStart[0] = lineDepthStart[1];
				lineDepthStart[1] = lineDepthStart[2];
				lineDepthEnd[0] = lineDepthEnd[1];
				lineDepthEnd[1] = lineDepthEnd[2];
				
				lineStartX[2] = a.x;// + (a2.x < a.x ? -1 : 1);
				//				lineStartX[1] = a.x + (b.x > a.x ? -1 : 1);
				lineY[2] = a.y;
				
				skipOtherLine = false;
				
				if (ignoreFirstLine) {
					ignoreFirstLine = false;
				} else {
					
					drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
//				drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], 'x', lineDepthStart[0], depth, depthBuffer);
				}
			} else {
				
				if(lineStartX[2] == lineEndX[2]) {
					//					sprintf(string, "lineStartX[2] == lineEndX");
//					mvprintw(line++, 0, "lineStartX[2] == lineEndX[2] @ y=%d", lineY[2]);
				} else
//					if (abs(a2.x - (a.x + (a2.x < a.x ? -1 : 1))) < abs(a2.x - lineStartX[2])) {
					if (abs(a2.x - a.x ) < abs(a2.x - lineStartX[2])) {
						//					mvaddstr(line++, 0, string);
						lineStartX[2] = a.x;// + (a2.x < a.x ? -1 : 1);
//					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
				}
				//				}
				errs[2] -= dx;
				if (a.y == a2.y) {
					skipOtherLine = true;
					
				}
			}
		}
		
		
		if (!skipOtherLine) {
			if (t2++ > 2) {
				setWithDepthBuffer( pts2[1], getp( pts2, errsNormalized2[1] + 0.5), depths2[1], depthBuffer);
			}
			pts2[0] = pts2[1];
			pts2[1] = pts2[2];
			pts2[2] = a2;
			errs2[0] = errs2[1];
			errs2[1] = errs2[2];
			errs2[2] = err2;
			
			errsNormalized2[0] = errsNormalized2[1];
			errsNormalized2[1] = errsNormalized2[2];
			errsNormalized2[2] = (double)(errs2[1])/(dx2 - dy2);
			if (sy2 == -1 ) {
				errsNormalized2[2] = 0.0 - errsNormalized2[2];
			}
			dxn2 = a2.x-c.x;
			dyn2 = a2.y-c.y;
			alpha2 = ((double)sqrt(dxn2*dxn2 + dyn2*dyn2))/lineMagnitudeSq2;
			depths2[0] = depths2[1];
			depths2[1] = depths2[2];
			depths2[2] = 1.0/(1.0/cc.z - alpha2*(1.0/cc.z - 1.0/oldAcZ ));
			
			lineDepthEnd[2] = depths2[2];

			if (a2.x == c.x && a2.y == c.y) {
				//			break;
				skipOtherLine = true;
				skipMainLine = false;
				if (a.x == b.x && a.y == b.y) {
					break;
				}
			}
		}
		if (!skipOtherLine) {
			e22 = 2*err2;
			if (e22 >= dy2) {
				err2 += dy2;
				a2.x += sx2;
			} else {
				errs2[2] -= dy2;
			}
			
			if (e22 <= dx2) {
				err2 += dx2;
				a2.y += sy2;
				
				lineEndX[2] = a2.x;
				skipMainLine = false;
			} else {
				errs2[2] -= dx2;
				
				skipMainLine = true;
				if (!(a.x == b.x && a.y == b.y)) {
					if (abs(a2.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
						lineEndX[1] = a2.x;
					}
				}
			}
		}
		
	}
	
	// add the final points for each line
	setWithDepthBuffer( pts[1], getp( pts, errsNormalized[1] + 0.5), depths[1], depthBuffer);
	setWithDepthBuffer( pts2[1], getp( pts2, errsNormalized2[1] + 0.5), depths2[1], depthBuffer);
	
	// draw the final fill line
		drawHorizonalLineWithDepthBuffer(lineStartX[1], lineEndX[1], lineY[1], fill, lineDepthStart[1], lineDepthEnd[1], depthBuffer);
		
}

Mat4D transpose( Mat4D& input ) {
	Mat4D result;
	
	for(int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.d[i][j] = input.d[j][i];
		}
	}
	
	return result;
}

Mat4D rotationFromAngleAndUnitAxis( double radians, Coordinates3D axis) {
	Mat4D result;
	double ct = cos(radians);
	double st = sin(radians);
	
	double omct = 1 - ct;
	
	result.d[0][0] = ct + axis.x*axis.x*omct;
	result.d[0][1] = axis.x*axis.y*omct - axis.z*st;
	result.d[0][2] = axis.x*axis.z*omct + axis.y*st;
	result.d[0][3] = 0;
	
	result.d[1][0] = axis.x*axis.y*omct + axis.z*st;
	result.d[1][1] = ct + axis.y*axis.y*omct;
	result.d[1][2] = axis.y*axis.z*omct - axis.x*st;
	result.d[1][3] = 0;
	
	result.d[2][0] = axis.x*axis.z*omct - axis.y*st;
	result.d[2][1] = axis.y*axis.z*omct + axis.x*st;
	result.d[2][2] = ct + axis.z*axis.z*omct;
	result.d[2][3] = 0;
	
	result.d[3][0] = 0;
	result.d[3][1] = 0;
	result.d[3][2] = 0;
	result.d[3][3] = 1;
	
	return result;
}

Mat4D translationMatrix( double x, double y, double z) {
	Mat4D result;
	memset(&result, 0, sizeof(result));
	
	
	result.d[0][0] = 1;
	result.d[1][1] = 1;
	result.d[2][2] = 1;
	result.d[3][3] = 1;
	
	result.d[0][3] = x;
	result.d[1][3] = y;
	result.d[2][3] = z;
	
	return result;
}

Mat4D scaleMatrix( double x, double y, double z ) {
	Mat4D result;
	memset(&result, 0, sizeof(result));
	
	result.d[0][0] = x;
	result.d[1][1] = y;
	result.d[2][2] = z;
	result.d[3][3] = 1;
	
	return result;
}

Mat3D matrixMultiply(Mat3D& a, Mat3D& b) {
	Mat3D result;
	memset(&result, 0, sizeof(result));
	
	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			double sum = 0;
			for (int k = 0; k < 3; k++) {
				sum += a.d[col][k] * b.d[k][row];
			}
			result.d[col][row] = sum;
		}
	}
	
	return result;
}

Mat4D matrixMultiply(Mat4D& a, Mat4D& b) {
	Mat4D result;
	memset(&result, 0, sizeof(result));
	
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			double sum = 0;
			for (int k = 0; k < 4; k++) {
				sum += a.d[col][k] * b.d[k][row];
			}
			result.d[col][row] = sum;
		}
	}
	
	return result;
}

Coordinates3D matrixVectorMultiple(Mat3D& rotation, Coordinates3D& vec) {
	Coordinates3D result = {0, 0 , 0};
	
	result.x  = rotation.d[0][0] * vec.x;
	result.x += rotation.d[0][1] * vec.y;
	result.x += rotation.d[0][2] * vec.z;
	
	result.y  = rotation.d[1][0] * vec.x;
	result.y += rotation.d[1][1] * vec.y;
	result.y += rotation.d[1][2] * vec.z;
	
	result.z  = rotation.d[2][0] * vec.x;
	result.z += rotation.d[2][1] * vec.y;
	result.z += rotation.d[2][2] * vec.z;
	
	return result;
}

Coordinates4D matrixVectorMultiply(Mat4D& rotation, Coordinates4D& vec) {
	Coordinates4D result;
	
	result.x  = rotation.d[0][0] * vec.x;
	result.x += rotation.d[0][1] * vec.y;
	result.x += rotation.d[0][2] * vec.z;
	result.x += rotation.d[0][3] * vec.w;
	
	result.y  = rotation.d[1][0] * vec.x;
	result.y += rotation.d[1][1] * vec.y;
	result.y += rotation.d[1][2] * vec.z;
	result.y += rotation.d[1][3] * vec.w;
	
	result.z  = rotation.d[2][0] * vec.x;
	result.z += rotation.d[2][1] * vec.y;
	result.z += rotation.d[2][2] * vec.z;
	result.z += rotation.d[2][3] * vec.w;
	
	result.w  = rotation.d[3][0] * vec.x;
	result.w += rotation.d[3][1] * vec.y;
	result.w += rotation.d[3][2] * vec.z;
	result.w += rotation.d[3][3] * vec.w;
	
	return result;
}

Coordinates3D matrixVectorMultiply(Mat4D& rotation, Coordinates3D& vec) {
	Coordinates3D result;
	
	result.x  = rotation.d[0][0] * vec.x;
	result.x += rotation.d[0][1] * vec.y;
	result.x += rotation.d[0][2] * vec.z;
//	result.x += rotation.d[0][3] * vec.w;
	
	result.y  = rotation.d[1][0] * vec.x;
	result.y += rotation.d[1][1] * vec.y;
	result.y += rotation.d[1][2] * vec.z;
//	result.y += rotation.d[1][3] * vec.w;
	
	result.z  = rotation.d[2][0] * vec.x;
	result.z += rotation.d[2][1] * vec.y;
	result.z += rotation.d[2][2] * vec.z;
//	result.z += rotation.d[2][3] * vec.w;
	
//	result.w  = rotation.d[3][0] * vec.x;
//	result.w += rotation.d[3][1] * vec.y;
//	result.w += rotation.d[3][2] * vec.z;
//	result.w += rotation.d[3][3] * vec.w;
	
	return result;
}

Coordinates3D crossProduct(Coordinates3D a, Coordinates3D b) {
	Coordinates3D result;
	
	result.x = a.y*b.z - a.z*b.y;
	result.y = a.z*b.x - a.x*b.z;
	result.z = a.x*b.y - a.y*b.x;
	
	return result;
}
Coordinates4D crossProduct(Coordinates4D a, Coordinates4D b) {
	Coordinates4D result;
	
	result.x = a.y*b.z - a.z*b.y;
	result.y = a.z*b.x - a.x*b.z;
	result.z = a.x*b.y - a.y*b.x;
	result.w = 1;	// shrug
	
	return result;
}

double dotProduct(const Coordinates3D& a, const Coordinates3D& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

double dotProduct(const Coordinates4D& a, const Coordinates4D& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
double dotProduct(const Coordinates4D& a, const Coordinates3D& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
double dotProduct(const Coordinates3D& a, const Coordinates4D& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

//double vectorSquared(Coordinates4D& a) {
//	return <#expression#>;
//}

Coordinates4D vectorAdd(Coordinates4D a, Coordinates4D b) {
	Coordinates4D result;
	
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	result.w = a.w + b.w;
	
	return result;
}

Coordinates3D vectorAdd(Coordinates3D a, Coordinates3D b) {
	Coordinates3D result;
	
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	
	return result;
}

Coordinates4D vectorSubtract(Coordinates4D a, Coordinates4D b) {
	Coordinates4D result;
	
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	result.w = a.w - b.w;
	
	return result;
}

Coordinates3D vectorSubtract(Coordinates3D a, Coordinates3D b) {
	Coordinates3D result;
	
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	
	return result;
}

Coordinates3D vectorScale(Coordinates3D a, double scale) {
	a.x *= scale;
	a.y *= scale;
	a.z *= scale;
	return a;
}

Coordinates3D interpolate(Coordinates3D& a, Coordinates3D& b, double factor) {
	Coordinates3D result = vectorSubtract(b, a);
	result.x *= factor;
	result.y *= factor;
	result.z *= factor;
	return vectorAdd(result, a);
}

Coordinates4D interpolate(Coordinates4D& a, Coordinates4D& b, double factor) {
	Coordinates4D result = vectorSubtract(b, a);
	result.x *= factor;
	result.y *= factor;
	result.z *= factor;
	result.w *= factor;
	return vectorAdd(result, a);
}


// https://www.comp.nus.edu.sg/~lowkl/publications/lowk_persp_interp_techrep.pdf
double perspectiveInterpolate(double& a, double& b, double& aDepth, double& bDepth, double& correctDepth, double& factor) {
	return correctDepth*(a/aDepth + factor*(b/bDepth - a/aDepth));
}
double perspectiveInterpolateInv(double& a, double& b, double& aInvDepth, double& bInvDepth, double& correctDepth, double& factor) {
	return correctDepth*(a*aInvDepth + factor*(b*bInvDepth - a*aInvDepth));
	//return correctDepth*(a*aInvDepth* (1 - factor) + factor*(b*bInvDepth));
}

Coordinates4D perspectiveInterpolate(Coordinates4D& a, Coordinates4D& b, double factor) {
	Coordinates4D result;
	
//	result.z = 1.0/(1.0/a.z + factor*(1.0/b.z - 1.0/a.z));	// This is Zt, where Z1 is with "a" and Z2 is with "b"
//	result.x = perspectiveInterpolate(a.x, b.x, a.z, b.z, result.z, factor);
//	result.y = perspectiveInterpolate(a.y, b.y, a.z, b.z, result.z, factor);
//	result.w = perspectiveInterpolate(a.w, b.w, a.z, b.z, result.z, factor);
	
	double aInvDepth = 1.0/a.z;
	double bInvDepth = 1.0/b.z;
	result.z = 1.0/(aInvDepth + factor*(bInvDepth - aInvDepth));	// This is Zt, where Z1 is with "a" and Z2 is with "b"
	result.x = perspectiveInterpolateInv(a.x, b.x, aInvDepth, bInvDepth, result.z, factor);
	result.y = perspectiveInterpolateInv(a.y, b.y, aInvDepth, bInvDepth, result.z, factor);
	result.w = perspectiveInterpolateInv(a.w, b.w, aInvDepth, bInvDepth, result.z, factor);
	
	
	return result;
}

Coordinates3D perspectiveInterpolate(Coordinates3D& a, Coordinates3D& b, double aDepth, double bDepth, double correctDepth, double factor) {
	Coordinates3D result;
	
//	result.z = perspectiveInterpolate(a.z, b.z, aDepth, bDepth, correctDepth, factor);
//	result.x = perspectiveInterpolate(a.x, b.x, aDepth, bDepth, correctDepth, factor);
//	result.y = perspectiveInterpolate(a.y, b.y, aDepth, bDepth, correctDepth, factor);
	//result.w = perspectiveInterpolate(a.w, b.w, a.z, b.z, result.z, factor);
	
	double aInvDepth = 1.0/aDepth;
	double bInvDepth = 1.0/bDepth;
	result.z = perspectiveInterpolateInv(a.z, b.z, aInvDepth, bInvDepth, correctDepth, factor);
	result.x = perspectiveInterpolateInv(a.x, b.x, aInvDepth, bInvDepth, correctDepth, factor);
	result.y = perspectiveInterpolateInv(a.y, b.y, aInvDepth, bInvDepth, correctDepth, factor);
	
	return result;
}

Coordinates4D perspectiveInterpolateInv(Coordinates4D& a, Coordinates4D& b, double aInvDepth, double bInvDepth, double correctDepth, double factor) {
	Coordinates4D result;

//	 aInvDepth = 1.0/a.z;
//	 bInvDepth = 1.0/b.z;
	result.z = 1.0/(aInvDepth + factor*(bInvDepth - aInvDepth));	// This is Zt, where Z1 is with "a" and Z2 is with "b"
	result.x = perspectiveInterpolateInv(a.x, b.x, aInvDepth, bInvDepth, result.z, factor);
	result.y = perspectiveInterpolateInv(a.y, b.y, aInvDepth, bInvDepth, result.z, factor);
	result.w = perspectiveInterpolateInv(a.w, b.w, aInvDepth, bInvDepth, result.z, factor);
	
	
	return result;
}

Coordinates3D perspectiveInterpolateInv(Coordinates3D& a, Coordinates3D& b, double aInvDepth, double bInvDepth, double correctDepth, double factor) {
	Coordinates3D result;

	result.z = perspectiveInterpolateInv(a.z, b.z, aInvDepth, bInvDepth, correctDepth, factor);
	result.x = perspectiveInterpolateInv(a.x, b.x, aInvDepth, bInvDepth, correctDepth, factor);
	result.y = perspectiveInterpolateInv(a.y, b.y, aInvDepth, bInvDepth, correctDepth, factor);
	
	return result;
}

Coordinates3D normalizeVector(Coordinates3D& input) {
	Coordinates3D result;
	double invMag = 1.0/sqrt(input.x*input.x + input.y*input.y + input.z*input.z);

	result.x = input.x*invMag;
	result.y = input.y*invMag;
	result.z = input.z*invMag;
	
	return result;
}

Coordinates3D normalizeVector(Coordinates4D& input) {
	Coordinates3D result;
	double invMag = 1.0/sqrt(input.x*input.x + input.y*input.y + input.z*input.z);

	result.x = input.x*invMag;
	result.y = input.y*invMag;
	result.z = input.z*invMag;
	
	return result;
}

Coordinates3D normalizeVectorFast(Coordinates3D& input) {
	Coordinates3D result;
	
	double magInv = Q_rsqrt(input.x*input.x + input.y*input.y + input.z*input.z);

	result.x = input.x * magInv;
	result.y = input.y * magInv;
	result.z = input.z * magInv;
	
	return result;
}

Coordinates3D normalizeVectorFast(Coordinates4D& input) {
	Coordinates3D result;
	
	double magInv = Q_rsqrt(input.x*input.x + input.y*input.y + input.z*input.z);
	
	result.x = input.x * magInv;
	result.y = input.y * magInv;
	result.z = input.z * magInv;
	
	return result;
}


Mat4D projectionMatrixOrtho(double width, double height, double zfar, double znear) {
	Mat4D result;
	memset(&result, 0, sizeof(result));
	
	result.d[0][0] = 1/width;
	result.d[1][1] = 1/height;
	result.d[2][2] = -2/(zfar-znear);
	result.d[2][3] = -(zfar+znear)/(zfar-znear);
//    result.d[3][2] = -(zfar+znear)/(zfar-znear);

	
	result.d[3][3] = 1;
	
	return result;
}

Mat4D projectionMatrixPerspective(double fov, double aspect, double zfar, double znear) {
	Mat4D result;
	memset(&result, 0, sizeof(result));
	
	result.d[1][1] = 1.0/tan(fov/2);
	result.d[0][0] = result.d[1][1] / aspect;
	
	result.d[2][2] = -(zfar+znear)/(zfar-znear);
	//	result.d[2][2] = -(zfar)/(zfar-znear);
	result.d[2][3] = -2*(zfar*znear)/(zfar-znear);
	// z' = -(zfar+znear)/(zfar-znear)*z + -2(zfar*znear)/(zfar-znear)

	result.d[3][2] = -1;
    
//    result.d[2][2] = -(zfar)/(zfar-znear);
//    result.d[2][3] = -(zfar*znear)/(zfar-znear);
//    result.d[3][2] =-1;
    
	result.d[3][3] = 0;
	
	// w' = -z
	
	// x = 1/(0,0)
	// y = 1/(1,1)
	// z = (z' - (2,3)) / (2,2)
	
	//	result = transpose(result);
	
	return result;
}


void fillPolygonNormals(Polygon4D* polygons, int count) {
	for (int i = 0; i < count; i++) {
		Polygon4D* polygon = &polygons[i];
		double crossProductMagSq = 0;
		Coordinates4D crossProd;
		int v1 = 0, v2 = 1, v3 = 2;
		while (crossProductMagSq == 0 && v1 < polygon->numVertices) {
			Coordinates4D vector1 = vectorSubtract(polygon->vertices[v2], polygon->vertices[v1]);
			Coordinates4D vector2 = vectorSubtract(polygon->vertices[v3], polygon->vertices[v1]);
			
			crossProd = crossProduct(vector1, vector2);
			crossProductMagSq = dotProduct(crossProd, crossProd);
			
			v1++;
			v2 = (v1+1) % polygon->numVertices;
			v3 = (v1+2) % polygon->numVertices;
		}
		if (crossProductMagSq != 0) {
//			polygon->normals[0] = normalizeVector(crossProd);
			polygon->normals[0] = normalizeVectorFast(crossProd);
			for (int v = 1; v < polygon->numVertices; v++) {
				polygon->normals[v] = polygon->normals[0];
			}
		}
	}
}


void rasterize(Coordinates4D* vertices, int edgeIndices[][2], int numEdges, Mat4D& windowTransform, DepthBuffer* depthBuffer) {
	// Frustum culling:
	
	uint8_t validEdge[numEdges];
	memset(validEdge, 1, sizeof(validEdge));
	Coordinates4D imagePoints[numEdges][2];
	for (int i = 0; i < numEdges; i++) {
		imagePoints[i][0] = vertices[edgeIndices[i][0]];
		imagePoints[i][1] = vertices[edgeIndices[i][1]];
		
		// Clip lines
		if ((validEdge[i] &= (liangBarskyHomogeneous(imagePoints[i][0], imagePoints[i][1]) != LINE_CLIP_REJECT)) == false) {
			// Reject, nothing to do, not even a transform
			continue;
		}
		
		// Convert homogeneous
		for (int v = 0; v < 2; v++) {
			if( imagePoints[i][v].w != 1) {
				imagePoints[i][v].x = imagePoints[i][v].x / imagePoints[i][v].w;
				imagePoints[i][v].y = imagePoints[i][v].y / imagePoints[i][v].w;
				imagePoints[i][v].z = imagePoints[i][v].z / imagePoints[i][v].w;
				imagePoints[i][v].w = 1;
			}
		}
		
		
		// Move viewport from -1,1 (all dimensions) to the full screen:
		for (int v = 0; v < 2; v++) {
			imagePoints[i][v] = matrixVectorMultiply(windowTransform, imagePoints[i][v]);
		}
	}
	
	// Rasterize
	for (int i = 0; i < numEdges; i++) {
		if(validEdge[i]) {
			
#ifdef PRINT_DEPTH
			printf("%f\t%f\n", imagePoints[i][0].z, imagePoints[i][1].z);
#endif
			lineWithDepthBuffer( imagePoints[i][0], imagePoints[i][1], depthBuffer);
		}
	}
}

Polygon4D clipTriangle(Coordinates4D input[3], int& line) {
	Polygon4D in, out;
	in.vertices[0] = input[0];
	in.vertices[1] = input[1];
	in.vertices[2] = input[2];
	in.numVertices = 3;
	clipPolygon(in, &out, line);
//	for (int i = 0; i < out.numVertices; i++) {
//		output[i] = out.vertices[i];
//	}
	return out;
	
//	LineClipStatus stat[3];
//
//	Coordinates4D edges[2][3];
//	edges[0][0] = input[0]; edges[1][0] = input[1];
//	edges[0][1] = input[1]; edges[1][1] = input[2];
//	edges[0][2] = input[2]; edges[1][2] = input[0];
//
//	//		// Clip lines
//	stat[0] = liangBarskyHomogeneous(edges[0][0], edges[1][0]);
//	stat[1] = liangBarskyHomogeneous(edges[0][1], edges[1][1]);
//	stat[2] = liangBarskyHomogeneous(edges[0][2], edges[1][2]);
//
//	output[0][0] = input[0];	// TBD move this to trivial accept
//	output[0][1] = input[1];
//	output[0][2] = input[2];
//
//	// Scenario 0: all lines rejected
//	if (stat[0] == LINE_CLIP_REJECT && stat[1] == LINE_CLIP_REJECT && stat[2] == LINE_CLIP_REJECT ) {
////		mvprintw(line++ + 10, 0, "triangle rejected");
//		return 0;
//	}
//
//	// Scenario 1: all edged accepted
//	if (stat[0] == LINE_CLIP_ACCEPT && stat[1] == LINE_CLIP_ACCEPT && stat[2] == LINE_CLIP_ACCEPT) {
////		mvprintw(line++ + 10, 0, "triangle accepted");
//
//		return 1;
//	}
//
//	// Scenario 2: one point clipped
//	if ((stat[0] | stat[1] | stat[2]) & LINE_CLIP_ACCEPT) {
////		mvprintw(line++ + 10, 0, "triangle single point clip");
//		Coordinates4D clippedEdge[2][2];
//		int acceptedEdge;
//
//		for (int i = 0, j = 0; i < 3; i++) {
//			if (stat[i] & LINE_CLIP_ACCEPT ) {
//				acceptedEdge = i;
//				continue;
//			}
//			if (stat[i] & LINE_CLIP_A) {
//				clippedEdge[0][j] = edges[0][i];
//				clippedEdge[1][j++] = edges[1][i];
//			} else {
//				clippedEdge[0][j] = edges[1][i];
//				clippedEdge[1][j++] = edges[0][i];
//			}
//
//		}
//
//
//		output[0][0] = edges[0][acceptedEdge];
//		output[0][1] = edges[1][acceptedEdge];
//		output[0][2] = clippedEdge[0][0];
//
//		output[1][0] = clippedEdge[0][0];
//		output[1][1] = edges[1][acceptedEdge];
//		output[1][2] = clippedEdge[0][1];
//
//		if (clippedEdge[0][0].x == clippedEdge[0][1].x ||
//			clippedEdge[0][0].y == clippedEdge[0][1].y) {
////			mvprintw(line++ + 10, 0, " - same clip plane");
//		} else {
////			output[2][0] = clippedEdge[0][1];
////			output[2][1] = clippedEdge[0][0];
//
////			output[2][2] = clippedEdge[0][1].;
//			return 2;
//		}
//
//		return 2;
//	}
//
////	mvprintw(line++ + 10, 0, "triangle clipped - unhandled");
//
//	return 1;
}

void rasterizeTriangle(Coordinates4D* vertices, int edgeIndices[][3], int numEdges, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill, int line) {
	// Frustum culling:
//	int line = 0;
//	bool validEdge[3];
//	Coordinates4D normalCheck;
//	uint8_t validTriangle[numEdges];
//	memset(validTriangle, 1, sizeof(validTriangle));
//	Coordinates4D imagePoints[numEdges][3];
	
//	Coordinates4D clippedTriangles[3][3];
	
	Polygon4D triangleAsPolygon;
	for (int i = 0; i < numEdges; i++) {
		
		
		
//		imagePoints[i][0] = vertices[edgeIndices[i][0]];
//		imagePoints[i][1] = vertices[edgeIndices[i][1]];
//		imagePoints[i][2] = vertices[edgeIndices[i][2]];
		
		triangleAsPolygon.vertices[0] = vertices[edgeIndices[i][0]];
		triangleAsPolygon.vertices[1] = vertices[edgeIndices[i][1]];
		triangleAsPolygon.vertices[2] = vertices[edgeIndices[i][2]];
		triangleAsPolygon.numVertices = 3;
		
		
		rasterizePolygon(&triangleAsPolygon, 1, windowTransform, depthBuffer, fill, line);
		
	}
////		int newTriangleCount = clipTriangle(imagePoints[i], clippedTriangles, line);
//		clippedTriangle = clipTriangle(imagePoints[i], line);
//
//		for (int t = 0; t < clippedTriangle.numVertices; t++) {
//
//
//		// Convert homogeneous
////		for (int v = 0; v < 3; v++) {
////			if( clippedTriangle. .w != 1) {
////				clippedTriangles[t][v].x = clippedTriangles[t][v].x / clippedTriangles[t][v].w;
////				clippedTriangles[t][v].y = -clippedTriangles[t][v].y / clippedTriangles[t][v].w;
////				clippedTriangles[t][v].z = clippedTriangles[t][v].z / clippedTriangles[t][v].w;
////				clippedTriangles[t][v].w = 1;
////			}
////		}
//
//
//		normalCheck = crossProduct(vectorSubtract( clippedTriangles[t][1], clippedTriangles[t][0]),
//								   vectorSubtract( clippedTriangles[t][2], clippedTriangles[t][1]) );
//		if (normalCheck.z > 0) {
////			mvprintw(line++ + 10, 0, "triangle culled - normal");
//			continue;
//		}
//
//
//		// Move viewport from -1,1 (all dimensions) to the full screen:
//		for (int v = 0; v < 3; v++) {
//			clippedTriangles[t][v] = matrixVectorMultiply(windowTransform, clippedTriangles[t][v]);
//		}
////	}
////
////	// Rasterize
////	for (int i = 0; i < numEdges; i++) {
////		if(validTriangle[i]) {
//
//#ifdef PRINT_DEPTH
//			printf("%f\t%f\n", clippedTriangles[t][0].z, clippedTriangles[t][1].z);
//#endif
//			//			lineWithDepthBuffer( clippedTriangles[t][0], clippedTriangles[t][1], depthBuffer);
//			//			lineWithDepthBuffer( clippedTriangles[t][1], clippedTriangles[t][2], depthBuffer);
//			//			lineWithDepthBuffer( clippedTriangles[t][2], clippedTriangles[t][0], depthBuffer);
//			triangleWithDepthBuffer( clippedTriangles[t][0], clippedTriangles[t][1], clippedTriangles[t][2], depthBuffer, fill);
////		}
//		}
//	}
}

void rasterizeQuads(Coordinates4D* vertices, int edgeIndices[][4], int count, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill, int &line) {
	Polygon4D quadAsPolygon;
	for (int i = 0; i < count; i++) {
		
		quadAsPolygon.vertices[0] = vertices[edgeIndices[i][0]];
		quadAsPolygon.vertices[1] = vertices[edgeIndices[i][1]];
		quadAsPolygon.vertices[2] = vertices[edgeIndices[i][2]];
		quadAsPolygon.vertices[3] = vertices[edgeIndices[i][3]];
		quadAsPolygon.numVertices = 4;
		
		rasterizePolygon(&quadAsPolygon, 1, windowTransform, depthBuffer, fill, line);
		
	}
}

void rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line) {
	Polygon4D quadAsPolygon;
	for (int i = 0; i < count; i++) {
		
		quadAsPolygon.vertices[0] = vertices[quadIndices[i][0]];
		quadAsPolygon.vertices[1] = vertices[quadIndices[i][1]];
		quadAsPolygon.vertices[2] = vertices[quadIndices[i][2]];
		quadAsPolygon.vertices[3] = vertices[quadIndices[i][3]];
		quadAsPolygon.numVertices = 4;
		
		fillPolygonNormals(&quadAsPolygon, 1);
		
//		rasterizePolygon(&quadAsPolygon, 1, windowTransform, depthBuffer, fill, line);
		rasterizePolygonsShader(&quadAsPolygon, 1, modelView, projection, viewport, userData, depthBuffer, fragmentShader, line);
		
	}
}

void rasterizePolygon(Polygon4D* polygons, int count, Mat4D& windowTransform, DepthBuffer* depthBuffer, char fill, int &line) {
	Polygon4D polygon;
	for (int i = 0; i < count; i++) {
		polygon = polygons[i];
		if (polygons[i].numVertices < 3) {
			continue;
		}
		
		double crossZ = 0;
		int v = 0, v2 = 1, v3 = 2;
		while (crossZ == 0 && v3 < polygon.numVertices) {
			double ax = polygon.vertices[v].x/polygon.vertices[v].w;
			double ay = -polygon.vertices[v].y/polygon.vertices[v].w;
			double bx = polygon.vertices[v2].x/polygon.vertices[v2].w;
			double by = -polygon.vertices[v2].y/polygon.vertices[v2].w;
			double cx = polygon.vertices[v3].x/polygon.vertices[v3].w;
			double cy = -polygon.vertices[v3].y/polygon.vertices[v3].w;
			v++;
			v2++;
			
			crossZ  = (bx-ax)*(cy-by) - (by-ay)*(cx-bx);
		}
		if(crossZ > 0) {
			continue;
		}
		
//		polygon = polygons[i];
		// Homogenous clipping:
		int valid = clipPolygon(polygons[i], &polygon, line);
		
		
		// Convert homogeneous
		for (int v = 0; v < polygon.numVertices; v++) {
			if( polygon.vertices[v].w != 1) {
				polygon.vertices[v].x = polygon.vertices[v].x / polygon.vertices[v].w;
				polygon.vertices[v].y = polygon.vertices[v].y / polygon.vertices[v].w;
				polygon.vertices[v].z = polygon.vertices[v].z / polygon.vertices[v].w;
				polygon.vertices[v].w = 1;
			}
		}
		
		// Normal check
		
		// Viewport
		for (int v = 0; v < polygon.numVertices; v++) {
			polygon.vertices[v] = matrixVectorMultiply(windowTransform, polygon.vertices[v]);
		}
		
		// Draw:
//		mvprintw(line++, 0, "Calling drawPolygon() with %d vertices", polygon.numVertices);
//		refresh();
		drawPolygon( polygon, depthBuffer, fill, line);
//		mvprintw(line++, 0, "Calling drawPolygon() Done");
//		refresh();
		
	}
}

void rasterizePolygonsShader(Polygon4D* polygons, int count,  Mat4D& modelView, Mat4D& projection, Mat4D& viewPort, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line) {
	Polygon4D polygon;
	Polygon4D polygonProjected;
	
	// x = 1/(0,0)
	// y = 1/(1,1)
	// z = (z' - (2,3)) / (2,2)
	Mat4D inverseProjection;
	memset(&inverseProjection, 0x00, sizeof(inverseProjection));
	inverseProjection.d[0][0] = 1.0/projection.d[0][0];
	inverseProjection.d[1][1] = 1.0/projection.d[1][1];
	
	
	// Original:
//	result.d[2][2] = -(zfar+znear)/(zfar-znear);
//	//	result.d[2][2] = -(zfar)/(zfar-znear);
//	result.d[2][3] = -2*(zfar*znear)/(zfar-znear);
//	// z' = -(zfar+znear)/(zfar-znear)*z + -2(zfar*znear)/(zfar-znear)
//
//	result.d[3][2] = -1;
//	result.d[3][3] = 0;
	
	
	// determinant = 1.0/(0 - (2,3)*(3,2)) = -1/((3,2)*(2,3))
//	inverseProjection.d[2][2] = projection.d[3][3]/(-bc); // always 0
	
	inverseProjection.d[2][3] = -1.0; // -(-(3,2))/((3,2)*(2,3)) = 1/(2,3) = -1
	
	inverseProjection.d[3][2] = 1.0/projection.d[2][3];  // -(-(3,2))/((3,2)*(2,3)) = 1/(2,3)

	inverseProjection.d[3][3] = projection.d[2][2]/projection.d[2][3];  // -((2,2))/((3,2)*(2,3)) = -(2,2)/((2,3)*-1) = (2,2)/(2,3)
	
	Polygon4D polygonRestored;
	
	for (int i = 0; i < count; i++) {
		if (polygons[i].numVertices < 3) {
			continue;
		}
		
		polygon = polygons[i];
		
		// ModelView
		for (int v = 0; v < polygon.numVertices; v++) {
			polygon.vertices[v] = matrixVectorMultiply(modelView, polygon.vertices[v]);
			polygon.normals[v] = matrixVectorMultiply(modelView, polygon.normals[v]);
		}
		
		
		
		// Normal clipping (w==1 in this case)
		Coordinates4D surfaceNormal, a, b;
		double dotCheck = 0;
		int v = 0, v2 = 1, v3 = 2;
		while (dotCheck == 0 && v3 < polygon.numVertices) {
			a = vectorSubtract(polygon.vertices[v2], polygon.vertices[v]);
			b = vectorSubtract(polygon.vertices[v3], polygon.vertices[v2]);
			
			surfaceNormal = crossProduct(a, b);
			
			dotCheck = dotProduct(surfaceNormal, polygon.vertices[v2]);
			///
			v++;
			v2++;
			v3++;
			
//			crossZ  = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
		}
		if(dotCheck > 0) {
//			mvprintw(line++, 0, "Normal culled: dotCheck = %f", dotCheck);
			continue;
		}
		
		polygonProjected = polygon;
		
		// Projection
		for (int v = 0; v < polygonProjected.numVertices; v++) {
			polygonProjected.vertices[v] = matrixVectorMultiply(projection, polygon.vertices[v]);
//			polygon.normals[v] = matrixVectorMultiply(projection, polygon.normals[v]);
		}
		
		
		
		// Homogenous clipping:
		int valid = clipPolygon(polygonProjected, &polygonProjected, line);
		
		// now get back the original corrdinates or separate handling
		polygonRestored.numVertices = polygonProjected.numVertices;
		for (int v = 0; v < polygonRestored.numVertices; v++) {
			polygonRestored.vertices[v] = matrixVectorMultiply(inverseProjection, polygonProjected.vertices[v]);
			polygonRestored.normals[v] = polygonProjected.normals[v];
		}
		
		// // Convert homogeneous
		for (int v = 0; v < polygonProjected.numVertices; v++) {
			if( polygonProjected.vertices[v].w != 1) {
				polygonProjected.vertices[v].x = polygonProjected.vertices[v].x / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].y = polygonProjected.vertices[v].y / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].z = polygonProjected.vertices[v].z / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].w = 1;
			}
		}
		
		// Viewport
		for (int v = 0; v < polygonProjected.numVertices; v++) {
			polygonProjected.vertices[v] = matrixVectorMultiply(viewPort, polygonProjected.vertices[v]);
		}
		
		// Draw:
//		mvprintw(line++, 0, "Calling drawPolygon() with %d vertices", polygon.numVertices);
//		refresh();
//		drawPolygon( polygon, depthBuffer, fill, line);
		drawPolygonShader( polygonProjected, polygonRestored, userData, depthBuffer, fragmentShader, line);
//		mvprintw(line++, 0, "Calling drawPolygon() Done");
//		refresh();
		
	}
}


// this is from: https://nanopdf.com/download/recall-liang-barsky-3d-clipping_pdf
LineClipStatus liangBarskyHomogeneous(Coordinates4D &point1, Coordinates4D &point2) {
	
	double tIn = 0.0, tOut = 1.0, tHit;
	double aBC[6], cBC[6];
	int aOutCode=0, cOutCode=0;
	
	int i;
	
	aBC[0] = point1.w + point1.x;
	aBC[1] = point1.w - point1.x;
	aBC[2] = point1.w + point1.y;
	aBC[3] = point1.w - point1.y;
	aBC[4] = point1.w + point1.z;
	aBC[5] = point1.w - point1.z;
	
	cBC[0] = point2.w + point2.x;
	cBC[1] = point2.w - point2.x;
	cBC[2] = point2.w + point2.y;
	cBC[3] = point2.w - point2.y;
	cBC[4] = point2.w + point2.z;
	cBC[5] = point2.w - point2.z;
	
	for (i = 0; i < 6; i++) {
		if (aBC[i] < 0) {
			aOutCode |= 0x01 << i;
		}
		if (cBC[i] < 0) {
			cOutCode |= 0x01 << i;
		}
	}
#ifdef PRINT_DEPTH
	printf("aOutCode = %d\tcOutCode = %d, point1.w = %f, point2.w = %f\n", aOutCode, cOutCode, point1.w	, point2.w);
#endif
	
	if((aOutCode & cOutCode)!=0)//trivialreject
		return LINE_CLIP_REJECT;
	if((aOutCode | cOutCode)==0)//trivialaccept
		return LINE_CLIP_ACCEPT;
	
	for(i=0;i<6;i++) { //clipagainsteachplane
		if(cBC[i]<0) {//Cisoutsidewalli(exitsotOut)
			tHit = aBC[i]/(aBC[i] - cBC[i]);
			tOut = tOut < tHit ? tOut : tHit;
		} if(aBC[i]<0) {
			tHit = aBC[i]/(aBC[i] - cBC[i]);
			tIn = tIn > tHit ? tIn : tHit;
		}
		if (tIn > tOut) {
			return LINE_CLIP_REJECT;
		}
	}
	
	LineClipStatus result = LINE_CLIP_ACCEPT;
	Coordinates4D tmp;
	if (aOutCode != 0) {
		tmp.x = point1.x + tIn*(point2.x - point1.x);
		tmp.y = point1.y + tIn*(point2.y - point1.y);
		tmp.z = point1.z + tIn*(point2.z - point1.z);
		tmp.w = point1.w + tIn*(point2.w - point1.w);
		result = LINE_CLIP_A;
	}
	if (cOutCode != 0) {
		point2.x = point1.x + tOut*(point2.x - point1.x);
		point2.y = point1.y + tOut*(point2.y - point1.y);
		point2.z = point1.z + tOut*(point2.z - point1.z);
		point2.w = point1.w + tOut*(point2.w - point1.w);
		result = result == LINE_CLIP_A ? LINE_CLIP_BOTH : LINE_CLIP_B;
	}
	
	if (aOutCode != 0) {
		point1 = tmp;
	}
	return result;
}


int mod(int a, int b)
{
	assert(b != 0);
	int r = a % b;
	return r < 0 ? r + b : r;
}


void drawPolygon( Polygon4D& poly, DepthBuffer* depthBuffer, char fill, int &line) {
	if (poly.numVertices < 3) {
//		mvprintw(line++, 0, "Not enough Polygon vertices: %d", poly.numVertices);
		return;
	}
	int minIndex = 0, maxIndex = 0;
	double minY = poly.vertices[0].y;
	double maxY = poly.vertices[0].y;
	for (int i = 1; i < poly.numVertices; i++) {
		if (poly.vertices[i].y < minY) {
			minY = poly.vertices[i].y;
			minIndex = i;
		}
		if (poly.vertices[i].y > maxY) {
			maxY = poly.vertices[i].y;
			maxIndex = i;
		}
	}
	
	
	char c = 'A';
	for ( int i = minIndex; ; ) {
//		mvprintw(line++, 0, "%d %c %f,%f %s", i, c, poly.vertices[i].x, poly.vertices[i].y, i == minIndex ? "Min" : (i == maxIndex ? "Max" : ""));
//		set(onlyXY(poly.vertices[i]), c++);
		setWithDepthBuffer(onlyXY(poly.vertices[i]), 'o', poly.vertices[i].z, depthBuffer);
//		refresh();
		i = mod(i + 1, poly.numVertices);
		if (i == minIndex)
			break;
	}
	
	
	int indexRight = minIndex;
	int indexLeft = minIndex;
	
	Coordinates2D ar = onlyXY(poly.vertices[minIndex]);
	Coordinates2D br = onlyXY(poly.vertices[mod(minIndex+1, poly.numVertices)]);
	double arDepth = poly.vertices[minIndex].z;
	double brDepth = poly.vertices[mod(minIndex+1, poly.numVertices)].z;
	
	Coordinates2D al = onlyXY(poly.vertices[minIndex]);
	Coordinates2D bl = onlyXY(poly.vertices[mod(minIndex-1, poly.numVertices)]);
	double alDepth = poly.vertices[minIndex].z;
	double blDepth = poly.vertices[mod(minIndex-1, poly.numVertices)].z;
	
	// The right line:
	int dxr = abs(br.x - ar.x);
	int sxr = ar.x < br.x ? 1 : -1;
	int dyr = -abs(br.y - ar.y);
	int syr = ar.y < br.y ? 1 : -1;
	int errr = dxr + dyr;
	int e2r;
	double lineMagnitudeSqr = sqrt(dxr*dxr + dyr*dyr);
	
	int errsr[3];
	double errsNormalizedr[3];
	Coordinates2D ptsr[3];
	
	double alphar;
	double dxnr, dynr;
	double depthsr[3];
	
	// The left line:
	int dxl = abs(bl.x - al.x);
	int sxl = al.x < bl.x ? 1 : -1;
	int dyl = -abs(bl.y - al.y);
	int syl = al.y < bl.y ? 1 : -1;
	int errl = dxl + dyl;
	int e2l;
	double lineMagnitudeSql = sqrt(dxl*dxl + dyl*dyl);
	
	int errsl[3];
	double errsNormalizedl[3];
	Coordinates2D ptsl[3];
	
	double alphal;
	double dxnl, dynl;
	double depthsl[3];
	
	
	
	int lineStartX[3] = {0,0,0};// = a.x - sx;
	int lineEndX[3] = {0,0,0};
	int lineY[3] = {0,0,0};// = a.y;
	double lineDepthStart[3] = {0,0,0};
	double lineDepthEnd[3] = {0,0,0};
	
	bool skipRightLine = false;
	bool skipLeftLine = true;
	
	bool rightComplete = false;
	bool leftComplete = false;
	
	int lineCount = 0;
	int tr = 0;
	int tl = 0;
	
	
//	while (indexRight != maxIndex || indexLeft != maxIndex) {
	while (!leftComplete || !rightComplete) {
//		refresh();
//		usleep(100000);
//
		if (al.y >= ar.y && !rightComplete) {
			skipLeftLine = true;
			skipRightLine = false;
		} else {
			skipLeftLine = false;
			skipRightLine = true;

		}
		
//		if (rightComplete) {
//			skipRightLine = true;
//		}
//		if (leftComplete) {
//			skipLeftLine = true;
//		}
		
		if (!skipRightLine && !rightComplete) {
			
			if (tr++ > 2) {
				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
				//				setWithDepthBuffer( pts[2], '.', depths[1], depthBuffer);
			}
			ptsr[0] = ptsr[1];
			ptsr[1] = ptsr[2];
			ptsr[2] = ar;
			
			errsr[0] = errsr[1];
			errsr[1] = errsr[2];
			errsr[2] = errr;
			
			errsNormalizedr[0] = errsNormalizedr[1];
			errsNormalizedr[1] = errsNormalizedr[2];
			errsNormalizedr[2] = (double)(errsr[1])/(dxr - dyr);
			if (syr == -1 ) {
				errsNormalizedr[2] = 0.0 - errsNormalizedr[2];
			}
			
			dxnr = ar.x-br.x;
			dynr = ar.y-br.y;
			alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr))/lineMagnitudeSqr;
			depthsr[0] = depthsr[1];
			depthsr[1] = depthsr[2];
//			depthsr[2] = brDepth - alphar*(brDepth - arDepth);
			depthsr[2] = 1.0/(1.0/brDepth - alphar*(1.0/brDepth - 1.0/arDepth));
			lineDepthStart[2] = depthsr[2];
//		}
		
		
		if (ar.x == br.x && ar.y == br.y && indexRight != maxIndex) {
			// new line time
			
			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
			
			
			indexRight = mod(indexRight + 1, poly.numVertices);
			int nextIndex = mod(indexRight + 1, poly.numVertices);
			
			if (indexRight == maxIndex) {
//				mvprintw(line++, 0, "Right line complete!");
//				refresh();
//				skipLeftLine = false;
//				skipRightLine = true;
				rightComplete = true;
				continue;
			}
			
//			mvprintw(line++, 0, "New Right line: %d - %d", indexRight, nextIndex);
			
			ar = onlyXY(poly.vertices[indexRight]);
			br = onlyXY(poly.vertices[mod(indexRight+1, poly.numVertices)]);
			arDepth = poly.vertices[indexRight].z;
			brDepth = poly.vertices[mod(indexRight+1, poly.numVertices)].z;
			
			dxr = abs(br.x - ar.x);
			sxr = ar.x < br.x ? 1 : -1;
			dyr = -abs(br.y - ar.y);
			syr = ar.y < br.y ? 1 : -1;
			errr = dxr + dyr;
//			e2r;
			lineMagnitudeSqr = sqrt(dxr*dxr + dyr*dyr);
			
			tr = 0;
			
//			break;
//			skipRightLine = false;
//			skipLeftLine = false;
			continue;
		}
		
//		if (!skipRightLine) {
			
			
			e2r = 2*errr;
			if (e2r >= dyr) {
				errr += dyr;
				ar.x += sxr;
			} else {
				errsr[2] -= dyr;
			}
			if (e2r <= dxr) {
				errr += dxr;
				ar.y += syr;
				lineStartX[0] = lineStartX[1];
				lineStartX[1] = lineStartX[2];
				lineY[0] = lineY[1];
				lineY[1] = lineY[2];
				lineEndX[0] = lineEndX[1];
				lineEndX[1] = lineEndX[2];
				
				lineDepthStart[0] = lineDepthStart[1];
				lineDepthStart[1] = lineDepthStart[2];
				lineDepthEnd[0] = lineDepthEnd[1];
				lineDepthEnd[1] = lineDepthEnd[2];
//				if (sx > 0 && sx2 > 0) {
//					lineStartX[2] = a.x + (a2.x < a.x ? 0 : 1);
//				} else {
				lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
					
//				}
					lineY[2] = ar.y;
				if (lineCount++ > 1) {
					
//					mvprintw(line++, 0, "line from %d-%d @ row %d", lineStartX[0], lineEndX[0], lineY[0]);
					if (fill != 0x00) {
						
						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
					}
				}
				
				
//				skipLeftLine = false;
				
				//				sprintf(string, "AB incremented, lineStartX[1] = %d", lineStartX[1]);
				//				mvaddstr(line++, 0, string);
			} else {
				if (abs(al.x - ar.x) < abs(al.x - lineStartX[2])) {
//					mvaddstr(line++, 0, string);
					lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
//					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
				}
				//				}
				errsr[2] -= dxr;
//				skipLeftLine = true;
			}
		}
		
		if (!skipLeftLine) {
			if (tl++ > 2) {
				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
				//				setWithDepthBuffer( pts2[2], '`', depths2[1], depthBuffer);
				
				
			}
			ptsl[0] = ptsl[1];
			ptsl[1] = ptsl[2];
			ptsl[2] = al;
			errsl[0] = errsl[1];
			errsl[1] = errsl[2];
			errsl[2] = errl;
			
			errsNormalizedl[0] = errsNormalizedl[1];
			errsNormalizedl[1] = errsNormalizedl[2];
			errsNormalizedl[2] = (double)(errsl[1])/(dxl - dyl);
			if (syl == -1 ) {
				errsNormalizedl[2] = 0.0 - errsNormalizedl[2];
			}
			dxnl = al.x-bl.x;
			dynl = al.y-bl.y;
			alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/lineMagnitudeSql;
			depthsl[0] = depthsl[1];
			depthsl[1] = depthsl[2];
			depthsl[2] = 1.0/(1.0/blDepth - alphal*(1.0/blDepth - 1.0/alDepth));
			lineDepthEnd[2] = depthsl[2];
			
//			if(al.x == bl.x && al.y == bl.y) {
//				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
//			}
//		}
		
		if (al.x == bl.x && al.y == bl.y && indexLeft != maxIndex) {
			// new line time
//			mvprintw(line++, 0, "New Left line!");
//
			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
			
			indexLeft = mod(indexLeft - 1, poly.numVertices);
			int nextIndex = mod(indexLeft - 1, poly.numVertices);
			
			if (indexLeft == maxIndex) {
//				mvprintw(line++, 0, "Left line complete!");
//				refresh();
//				skipLeftLine = true;
//				skipRightLine = false;
				leftComplete = true;
				continue;
			}
//			mvprintw(line++, 0, "New Left line: %d - %d", indexLeft, nextIndex);
			
			al = onlyXY(poly.vertices[indexLeft]);
			bl = onlyXY(poly.vertices[mod(indexLeft-1, poly.numVertices)]);
			alDepth = poly.vertices[indexLeft].z;
			blDepth = poly.vertices[mod(indexLeft-1, poly.numVertices)].z;
			
			dxl = abs(bl.x - al.x);
			sxl = al.x < bl.x ? 1 : -1;
			dyl = -abs(bl.y - al.y);
			syl = al.y < bl.y ? 1 : -1;
			errl = dxl + dyl;
//			e2l;
			lineMagnitudeSql = sqrt(dxl*dxl + dyl*dyl);
			
//			skipLeftLine = false;
//			skipRightLine = true;
			
			tl = 0;
			
			continue;
		}
		
//		if (!skipLeftLine) {
			e2l = 2*errl;
			if (e2l >= dyl) {
				errl += dyl;
				
				al.x += sxl;
			} else {
				errsl[2] -= dyl;
			}
			
			if (e2l <= dxl) {
				errl += dxl;
				al.y += syl;
				lineEndX[2] = al.x;
//				skipRightLine = false;
	
				
				
				
			} else {
				errsl[2] -= dxl;
//				skipRightLine = true;
				
				if (!(ar.x == br.x && ar.y == br.y)) {
					if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
						lineEndX[1] = al.x;
					}
				}
			}
		}
		
		
	}
	if (tr > 2) {
		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
	}
	if (tl > 2) {
		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
	}
	// draw the final fill line
	
	if (fill != 0x00 && lineCount > 1) {
//		mvprintw(line++, 0, "final line: %d,%d", lineStartX[1],lineEndX[1]);
		drawHorizonalLineWithDepthBuffer(lineStartX[1], lineEndX[1], lineY[1], fill, lineDepthStart[1], lineDepthEnd[1], depthBuffer);
	}
		
}


//void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, Mat4D modelView,  void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line) {
//	if (poly.numVertices < 3) {
////		mvprintw(line++, 0, "Not enough Polygon vertices: %d", poly.numVertices);
//		return;
//	}
//	int minIndex = 0, maxIndex = 0;
//	double minY = poly.vertices[0].y;
//	double maxY = poly.vertices[0].y;
//	for (int i = 1; i < poly.numVertices; i++) {
//		if (poly.vertices[i].y < minY) {
//			minY = poly.vertices[i].y;
//			minIndex = i;
//		}
//		if (poly.vertices[i].y > maxY) {
//			maxY = poly.vertices[i].y;
//			maxIndex = i;
//		}
//	}
//
//
////	char c = 'A';
//	for ( int i = minIndex; ; ) {
////		mvprintw(line++, 0, "%d %c %f,%f %s", i, c, poly.vertices[i].x, poly.vertices[i].y, i == minIndex ? "Min" : (i == maxIndex ? "Max" : ""));
////		set(onlyXY(poly.vertices[i]), c++);
//		setWithDepthBuffer(onlyXY(poly.vertices[i]), 'o', poly.vertices[i].z, depthBuffer);
////		refresh();
//		i = mod(i + 1, poly.numVertices);
//		if (i == minIndex)
//			break;
//	}
//
//
//	int indexRight = minIndex;
//	int indexLeft = minIndex;
//
//	Coordinates2D ar = onlyXY(poly.vertices[minIndex]);
//	Coordinates2D br = onlyXY(poly.vertices[mod(minIndex+1, poly.numVertices)]);
//	Coordinates3D arNormal = restored.normals[minIndex];
//	Coordinates4D ar3dPoint = restored.vertices[minIndex];
//	Coordinates3D brNormal = restored.normals[mod(minIndex+1, poly.numVertices)];
//	Coordinates4D br3dPoint = restored.vertices[mod(minIndex+1, poly.numVertices)];
//	double arDepth = poly.vertices[minIndex].z;
//	double brDepth = poly.vertices[mod(minIndex+1, poly.numVertices)].z;
//
//	Coordinates2D al = onlyXY(poly.vertices[minIndex]);
//	Coordinates2D bl = onlyXY(poly.vertices[mod(minIndex-1, poly.numVertices)]);
//	Coordinates3D alNormal = restored.normals[minIndex];
//	Coordinates4D al3dPoint = restored.vertices[minIndex];
//	Coordinates3D blNormal = restored.normals[mod(minIndex-1, poly.numVertices)];
//	Coordinates4D bl3dPoint = restored.vertices[mod(minIndex-1, poly.numVertices)];
//	double alDepth = poly.vertices[minIndex].z;
//	double blDepth = poly.vertices[mod(minIndex-1, poly.numVertices)].z;
//
//	// The right line:
//	int dxr = abs(br.x - ar.x);
//	int sxr = ar.x < br.x ? 1 : -1;
//	int dyr = -abs(br.y - ar.y);
//	int syr = ar.y < br.y ? 1 : -1;
//	int errr = dxr + dyr;
//	int e2r;
//	double lineMagnitudeSqr = sqrt(dxr*dxr + dyr*dyr);
//
//	int errsr[3];
//	double errsNormalizedr[3];
//	Coordinates2D ptsr[3];
//
//	double alphar;
//	double dxnr, dynr;
//	double depthsr[3];
//
//	// The left line:
//	int dxl = abs(bl.x - al.x);
//	int sxl = al.x < bl.x ? 1 : -1;
//	int dyl = -abs(bl.y - al.y);
//	int syl = al.y < bl.y ? 1 : -1;
//	int errl = dxl + dyl;
//	int e2l;
//	double lineMagnitudeSql = sqrt(dxl*dxl + dyl*dyl);
//
//	int errsl[3];
//	double errsNormalizedl[3];
//	Coordinates2D ptsl[3];
//
//	double alphal;
//	double dxnl, dynl;
//	double depthsl[3];
//
//
//
//	int lineStartX[3] = {0,0,0};// = a.x - sx;
//	int lineEndX[3] = {0,0,0};
//	int lineY[3] = {0,0,0};// = a.y;
//	double lineDepthStart[3] = {0,0,0};
//	double lineDepthEnd[3] = {0,0,0};
//
//	Coordinates4D point3dr[3];
//	Coordinates4D point3dl[3];
//	Coordinates3D normalr[3];
//	Coordinates3D normall[3];
//
//	bool skipRightLine = false;
//	bool skipLeftLine = true;
//
//	bool rightComplete = false;
//	bool leftComplete = false;
//
//	int lineCount = 0;
//	int tr = 0;
//	int tl = 0;
//
//
////	while (indexRight != maxIndex || indexLeft != maxIndex) {
//	while (!leftComplete || !rightComplete) {
////		refresh();
////		usleep(500000);
////
//		if (al.y >= ar.y && !rightComplete) {
//			skipLeftLine = true;
//			skipRightLine = false;
//		} else {
//			skipLeftLine = false;
//			skipRightLine = true;
//
//		}
//
////		if (rightComplete) {
////			skipRightLine = true;
////		}
////		if (leftComplete) {
////			skipLeftLine = true;
////		}
//
//		if (!skipRightLine && !rightComplete) {
//
//			if (tr++ > 2) {
//				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
//				//				setWithDepthBuffer( pts[2], '.', depths[1], depthBuffer);
//			}
//			ptsr[0] = ptsr[1];
//			ptsr[1] = ptsr[2];
//			ptsr[2] = ar;
//
//			errsr[0] = errsr[1];
//			errsr[1] = errsr[2];
//			errsr[2] = errr;
//
//			errsNormalizedr[0] = errsNormalizedr[1];
//			errsNormalizedr[1] = errsNormalizedr[2];
//			errsNormalizedr[2] = (double)(errsr[1])/(dxr - dyr);
//			if (syr == -1 ) {
//				errsNormalizedr[2] = 0.0 - errsNormalizedr[2];
//			}
//
//			dxnr = ar.x-br.x;
//			dynr = ar.y-br.y;
//			alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr))/lineMagnitudeSqr;
//			depthsr[0] = depthsr[1];
//			depthsr[1] = depthsr[2];
////			depthsr[2] = brDepth - alphar*(brDepth - arDepth);
//			depthsr[2] = 1.0/(1.0/brDepth + alphar*(1.0/arDepth - 1.0/brDepth));
//			lineDepthStart[2] = depthsr[2];
//
//
//			normalr[0] = normalr[1];
//			normalr[1] = normalr[2];
//			normalr[2] = interpolate(brNormal, arNormal, alphar);
//
//			point3dr[0] = point3dr[1];
//			point3dr[1] = point3dr[2];
////			point3dr[2] = interpolate(br3dPoint, ar3dPoint, alphar);
//			point3dr[2] = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, depthsr[2], alphar);
////			point3dr[2].z = 1.0/(1.0/br3dPoint.z + alphar*(1.0/ar3dPoint.z - 1.0/br3dPoint.z));
////		}
//
//
//		if (ar.x == br.x && ar.y == br.y && indexRight != maxIndex) {
//			// new line time
//
//			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
//
//
//			indexRight = mod(indexRight + 1, poly.numVertices);
//			int nextIndex = mod(indexRight + 1, poly.numVertices);
//
//			if (indexRight == maxIndex) {
////				mvprintw(line++, 0, "Right line complete!");
////				refresh();
////				skipLeftLine = false;
////				skipRightLine = true;
//				rightComplete = true;
//				continue;
//			}
//
////			mvprintw(line++, 0, "New Right line: %d - %d", indexRight, nextIndex);
//
//			ar = onlyXY(poly.vertices[indexRight]);
//			br = onlyXY(poly.vertices[nextIndex]);
//			arNormal = restored.normals[indexRight];
//			ar3dPoint = restored.vertices[indexRight];
//			brNormal = restored.normals[nextIndex];
//			br3dPoint = restored.vertices[nextIndex];
//
//			arDepth = poly.vertices[indexRight].z;
//			brDepth = poly.vertices[nextIndex].z;
//
//			dxr = abs(br.x - ar.x);
//			sxr = ar.x < br.x ? 1 : -1;
//			dyr = -abs(br.y - ar.y);
//			syr = ar.y < br.y ? 1 : -1;
//			errr = dxr + dyr;
////			e2r;
//			lineMagnitudeSqr = sqrt(dxr*dxr + dyr*dyr);
//
//			tr = 0;
//
////			break;
////			skipRightLine = false;
////			skipLeftLine = false;
//			continue;
//		}
//
////		if (!skipRightLine) {
//
//
//			e2r = 2*errr;
//			if (e2r >= dyr) {
//				errr += dyr;
//				ar.x += sxr;
//			} else {
//				errsr[2] -= dyr;
//			}
//			if (e2r <= dxr) {
//				errr += dxr;
//				ar.y += syr;
//
//				lineStartX[0] = lineStartX[1];
//				lineStartX[1] = lineStartX[2];
//				lineY[0] = lineY[1];
//				lineY[1] = lineY[2];
//				lineEndX[0] = lineEndX[1];
//				lineEndX[1] = lineEndX[2];
//
//				lineDepthStart[0] = lineDepthStart[1];
//				lineDepthStart[1] = lineDepthStart[2];
//				lineDepthEnd[0] = lineDepthEnd[1];
//				lineDepthEnd[1] = lineDepthEnd[2];
////				if (sx > 0 && sx2 > 0) {
////					lineStartX[2] = a.x + (a2.x < a.x ? 0 : 1);
////				} else {
//				lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
//
////				}
//					lineY[2] = ar.y;
////				if (lineCount++ > 0) {
////
//////					if (lineCount == 2) {
////						mvprintw(line++, 1, "line from %d-%d @ row %d", lineStartX[1], lineEndX[1], lineY[1]);
//////					}
//////					if (fill != 0x00) {
////					static bool shoulDplot = true;
////					if (shoulDplot) {
////
////						shoulDplot = false;
//////						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
////					drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
////					} else {
////						shoulDplot = true;
////					}
//////					}
////				}
//
//
////				skipLeftLine = false;
//
//				//				sprintf(string, "AB incremented, lineStartX[1] = %d", lineStartX[1]);
//				//				mvaddstr(line++, 0, string);
//			} else {
//				if (abs(al.x - ar.x) < abs(al.x - lineStartX[2])) {
////					mvaddstr(line++, 0, string);
//					lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
////					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
//				}
//				//				}
//				errsr[2] -= dxr;
////				skipLeftLine = true;
//			}
//		}
//
//		if (!skipLeftLine) {
//			if (tl++ > 2) {
//				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
//				//				setWithDepthBuffer( pts2[2], '`', depths2[1], depthBuffer);
//
//
//			}
//			ptsl[0] = ptsl[1];
//			ptsl[1] = ptsl[2];
//			ptsl[2] = al;
//			errsl[0] = errsl[1];
//			errsl[1] = errsl[2];
//			errsl[2] = errl;
//
//			errsNormalizedl[0] = errsNormalizedl[1];
//			errsNormalizedl[1] = errsNormalizedl[2];
//			errsNormalizedl[2] = (double)(errsl[1])/(dxl - dyl);
//			if (syl == -1 ) {
//				errsNormalizedl[2] = 0.0 - errsNormalizedl[2];
//			}
//			dxnl = al.x-bl.x;
//			dynl = al.y-bl.y;
//			alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/lineMagnitudeSql;
//			depthsl[0] = depthsl[1];
//			depthsl[1] = depthsl[2];
////			depthsl[2] = blDepth - alphal*(blDepth - alDepth);
//			depthsl[2] = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
////			lineDepthEnd[2] = depthsl[2];
//
//
////			normall[0] = normall[1];
////			normall[1] = normall[2];
////			normall[2] = interpolate(blNormal, alNormal, alphal);
//
////			point3dl[0] = point3dl[1];
////			point3dl[1] = point3dl[2];
//////			point3dl[2] = interpolate(bl3dPoint, al3dPoint, alphal);
//////
//////			point3dl[2].z = 1.0/(1.0/bl3dPoint.z + alphal*(1.0/al3dPoint.z - 1.0/bl3dPoint.z));
////
////
////			point3dl[2] = perspectiveInterpolate(bl3dPoint, al3dPoint, blDepth, alDepth, depthsl[2], alphal);
////
////			if(al.x == bl.x && al.y == bl.y) {
////				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
////			}
////		}
//
//		if (al.x == bl.x && al.y == bl.y && indexLeft != maxIndex) {
//			// new line time
////			mvprintw(line++, 0, "New Left line!");
////
//			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
//
//			indexLeft = mod(indexLeft - 1, poly.numVertices);
//			int nextIndex = mod(indexLeft - 1, poly.numVertices);
//
//			if (indexLeft == maxIndex) {
////				mvprintw(line++, 0, "Left line complete!");
////				refresh();
////				skipLeftLine = true;
////				skipRightLine = false;
//				leftComplete = true;
//				continue;
//			}
////			mvprintw(line++, 0, "New Left line: %d - %d", indexLeft, nextIndex);
//
////			int nextIndex = mod(indexLeft-1, poly.numVertices);
//			al = onlyXY(poly.vertices[indexLeft]);
//			bl = onlyXY(poly.vertices[nextIndex]);
//
//			alNormal = restored.normals[indexLeft];
//			al3dPoint = restored.vertices[indexLeft];
//			blNormal = restored.normals[nextIndex];
//			bl3dPoint = restored.vertices[nextIndex];
//
//			alDepth = poly.vertices[indexLeft].z;
//			blDepth = poly.vertices[nextIndex].z;
//
//			dxl = abs(bl.x - al.x);
//			sxl = al.x < bl.x ? 1 : -1;
//			dyl = -abs(bl.y - al.y);
//			syl = al.y < bl.y ? 1 : -1;
//			errl = dxl + dyl;
////			e2l;
//			lineMagnitudeSql = sqrt(dxl*dxl + dyl*dyl);
//
////			skipLeftLine = false;
////			skipRightLine = true;
//
//			tl = 0;
//
//			continue;
//		}
//
////		if (!skipLeftLine) {
//			e2l = 2*errl;
//			if (e2l >= dyl) {
//				errl += dyl;
//
//				al.x += sxl;
//			} else {
//				errsl[2] -= dyl;
//			}
//
//			if (e2l <= dxl) {
//				errl += dxl;
//				al.y += syl;
//				lineEndX[2] = al.x;
////				skipRightLine = false;
//
//
//				if (lineCount++ > 0) {
//
//					dxnl = lineEndX[1]-bl.x;
//					dynl = (al.y-1)-bl.y;
//					alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/lineMagnitudeSql;
////					depthsl[0] = depthsl[1];
////					depthsl[1] = depthsl[2];
////		//			depthsl[2] = blDepth - alphal*(blDepth - alDepth);
////					depthsl[2] = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
////					lineDepthEnd[2] = depthsl[2];
//					lineDepthEnd[1] = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
//
////					normall[0] = normall[1];
////					normall[1] = normall[2];
//					normall[1] = interpolate(blNormal, alNormal, alphal);
//
//					point3dl[1] = perspectiveInterpolate(bl3dPoint, al3dPoint, blDepth, alDepth, lineDepthEnd[1], alphal);
//
////						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
//					//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//					drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], al.y-1, lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//
//				}
//
//
//
//			} else {
//				errsl[2] -= dxl;
////				skipRightLine = true;
//
//				//if (!(ar.x == br.x && ar.y == br.y)) {
//					if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
////						mvprintw(line++, 1, "Correcting lineEndX[1]");
//						lineEndX[1] = al.x;
//					}
//				//}
//			}
//		}
//
//
//	}
//	if (tr > 2) {
//		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
//	}
//	if (tl > 2) {
//		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
//	}
//	// draw the final fill line
//
//	if ( lineCount > 0 ) {
////		mvprintw(line++, 0, "final line: %d,%d", lineStartX[1],lineEndX[1]);
////		drawHorizonalLineWithDepthBuffer(lineStartX[1], lineEndX[1], lineY[1], fill, lineDepthStart[1], lineDepthEnd[1], depthBuffer);
//
////		mvprintw(line++, 1, "Final line from %d-%d @ row %d", lineStartX[1], lineEndX[1], lineY[1]);
//		//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//	}
//
//}

//typedef union _Vector4D {
//	Coordinates4D c;
//	double d[4];
//} Vector4D;

#define CLIP_PLANE_W (0.00001)

int clipW(Polygon4D input, Polygon4D* output, int& line) {
	output->numVertices = 0;
	
	double factor;
	
	char previousDot;
	char currentDot;
	
	Coordinates4D *current = &input.vertices[0];
	Coordinates4D *prior = &input.vertices[input.numVertices-1];
	previousDot = prior->w < CLIP_PLANE_W ? -1 : 1;
	for (; current != &input.vertices[input.numVertices]; ) {
		currentDot = (current->w < CLIP_PLANE_W) ? -1 : 1;
//			factor = current->w/current->x;
		// duplicate and shift
		
		if (previousDot * currentDot < 0) {
			mvprintw(line++, 0, " - - Clip against w=%f", CLIP_PLANE_W);
			factor = (CLIP_PLANE_W - prior->w) / (prior->w - current->w);
			Coordinates4D diff = vectorSubtract(*current, *prior);
			diff.x *= factor;
			diff.y *= factor;
			diff.z *= factor;
			diff.w *= factor;
			diff = vectorAdd(*prior, diff);
			
			output->vertices[output->numVertices] = diff;
			output->numVertices++;
		}
		
		if (currentDot > 0) {
//				mvprintw(line++, 0, " - - Vertex Insertion");
			output->vertices[output->numVertices] = *current;
			output->numVertices++;
		}
		
		previousDot = currentDot;
		prior = current;
		current++;
	}
	return 1;
}

int clipPlane(Polygon4D input, Polygon4D* output, int axis, int plane, int& line) {
	output->numVertices = 0;
	
	double factor;
	
//	char previousDot;
//	char currentDot;
//	int previousDot;
//	int currentDot;
	bool previousDot;
	bool currentDot;
	
	Coordinates4D *current = &input.vertices[0];
	Coordinates3D *currentNormal = &input.normals[0];
	ColorRGBA *currentColor = &input.colors[0];
	Coordinates4D *prior = &input.vertices[input.numVertices-1];
	previousDot = ((double*)prior)[axis] <= prior->w ;
	for (; current != &input.vertices[input.numVertices]; ) {
		currentDot = ((double*)current)[axis] <= current->w;
//			factor = current->w/current->x;
		// duplicate and shift
		
		if (previousDot != currentDot ) {
//			mvprintw(line++, 0, " - - Clip against plane %d", axis);
			factor = (prior->w - ((double*)prior)[axis])/((prior->w - ((double*)prior)[axis]) - (current->w - ((double*)current)[axis])) ;
			Coordinates4D diff = vectorSubtract(*current, *prior);
			diff.x *= factor;
			diff.y *= factor;
			diff.z *= factor;
			diff.w *= factor;
			diff = vectorAdd(*prior, diff);
			
			output->vertices[output->numVertices] = diff;
			output->normals[output->numVertices] = *currentNormal; // this should be interpolated as well
			output->colors[output->numVertices] = *currentColor; // this should be interpolated as well
			output->numVertices++;
		}
		
		if (currentDot ) {
//				mvprintw(line++, 0, " - - Vertex Insertion");
			output->vertices[output->numVertices] = *current;
			output->normals[output->numVertices] = *currentNormal;
			output->colors[output->numVertices] = *currentColor;
			output->numVertices++;
		}
		
		previousDot = currentDot;
		prior = current;
		current++;
		currentNormal++;
		currentColor++;
	}
	
	input = *output;
	output->numVertices = 0;
	
	current = &input.vertices[0];
	currentNormal = &input.normals[0];
	currentColor = &input.colors[0];
	prior = &input.vertices[input.numVertices-1];
	previousDot = -((double*)prior)[axis] <= prior->w;
	for (; current != &input.vertices[input.numVertices]; ) {
		currentDot = -((double*)current)[axis] <= current->w;
//			factor = current->w/current->x;
		// duplicate and shift
		
		if (previousDot != currentDot ) {
//			mvprintw(line++, 0, " - - Clip against w=x");
//			mvprintw(line++, 0, " - - Clip against plane -%d", axis);
			factor = (prior->w + ((double*)prior)[axis])/((prior->w + ((double*)prior)[axis]) - (current->w + ((double*)current)[axis])) ;
			Coordinates4D diff = vectorSubtract(*current, *prior);
			diff.x *= factor;
			diff.y *= factor;
			diff.z *= factor;
			diff.w *= factor;
			diff = vectorAdd(*prior, diff);
			
			output->vertices[output->numVertices] = diff;
			output->normals[output->numVertices] = *currentNormal; // this should be interpolated as well
			output->colors[output->numVertices] = *currentColor; // this should be interpolated as well
			output->numVertices++;
		}
		
		if (currentDot) {
//				mvprintw(line++, 0, " - - Vertex Insertion");
			output->vertices[output->numVertices] = *current;
			output->normals[output->numVertices] = *currentNormal;
			output->colors[output->numVertices] = *currentColor;
			output->numVertices++;
		}
		
		previousDot = currentDot;
		prior = current;
		current++;
		currentNormal++;
		currentColor++;
	}
	return 1;
}

int clipPolygon(Polygon4D& input, Polygon4D* output, int& line) {
	
//	mvprintw(line++, 0, "Polygon %d");
//	for (int i = 0; i < input.numVertices; i++) {
//		if (input.vertices[i].w <= 0) {
//			mvprintw(line++, 0, "Polygon W clip needed: %f", input.vertices[i].w);
//			return 0;
//		}
//	}
	
	// The following functions are from https://fabiensanglard.net/polygon_codec/
	
//	clipW(input, output, line);	// I don't think this is necessary (well it is, but it's an ultra rare corner case?)
	clipPlane(input, output, 0, 1, line);
	clipPlane(*output, output, 1, 1, line);
	clipPlane(*output, output, 2, 1, line);

	
	
	return 1;

}


void asTexImage2d(FrameBuffer* fbo, FrameBufferType type, int width, int height) {
	fbo->rows = height;
	fbo->cols = width;
	fbo->type = type;
    
	int channels = 0;
	switch (type) {
		case FBT_RGB:
			channels = 3;
			break;
		case FBT_RGBA:
			channels = 4;
			break;
		case FBT_DEPTH:
			channels = 1 * sizeof(double);
			break;
	}
	if(fbo->data) {
		free(fbo->data);
	}
	fbo->data = (uint8_t*)malloc(width*height*channels);
}



void setFrameBufferRGBA(int x, int y, FrameBuffer* fbo, const ColorRGBA& value) {
	if (x < 0 || y < 0 || x >= fbo->cols || y >= fbo->rows) {
		return;
	}
	
	((ColorRGBA*)fbo->data)[y*fbo->cols + x] = value;
	
}


ColorRGBA interpolate(ColorRGBA& a, ColorRGBA& b, double factor) {
	ColorRGBA result;
	result.r = a.r + factor*(b.r - a.r);
	result.g = a.g + factor*(b.g - a.g);
	result.b = a.b + factor*(b.b - a.b);
	result.a = a.a + factor*(b.a - a.a);
	
	return result;
}

void drawHorizonalLineRGBA(double x1, double x2, int y, ColorRGBA& color1, ColorRGBA& color2, FrameBuffer* fbo) {
//	Coordinates2D pixel;
//	int diff = x2-x1;
//	if (abs(diff) < 2) {
//		return;
//	}
	if ((int)x1 == (int)x2) {
		
		setFrameBufferRGBA(x1, y, fbo, interpolate(color1, color2, 0.5));
		return;
	}
	
	int increment = x2 > x1 ? 1.0 : -1.0;
//	double factor;
	
//	Coordinates4D point;
//	Coordinates3D normal;
//	double depth;
	
//	pixel.y = y;
//	for (int i = x1+increment; i != x2; i += increment) {
	int i;
	int end = x2+increment;
	
	// i goes from x1 to x2
	// c goes from 0 to 1, or from i-x1 to
	double cF = increment/((double)x2 - (double)x1);
	ColorRGBA color;
	double factor = ( x1-floor(x1))*cF;
	
	for (i = x1; i != end; i += increment, factor += cF) {
//		pixel.x = i;
		
//		color = interpolate(color1, color2, ((double)i-x1)*cF);
		color = interpolate(color1, color2, factor);
//		color = interpolate(color1, color2, (double)(xDouble-x1)*cF);
		
		setFrameBufferRGBA(i, y, fbo, color);
		
	}
//	setFrameBufferRGBA(i, y, fbo, color);
	
}

void sortedIndicesByY(Coordinates4D* points, int* indices, int count, int level=0) {
	
	if (count == level) {
		return;
	}
	indices[level] = level;
	double minY = points[level].y;
	
	for(int i = level+1; i < count; i++) {
		if (points[i].y < minY) {
			minY = points[i].y;
			indices[i] = i;
		}
	}
	
	sortedIndicesByY(points, indices, count, level+1);
	
}


bool compareCoordinates(Coordinates4D i1, Coordinates4D i2)
{
	return (i1.y < i2.y);
}

bool compareCoordinatesFrag(FragmentInfo i1, FragmentInfo i2)
{
	return (i1.location3D.y < i2.location3D.y);
}


//void triangleFill(Coordinates4D& p1, Coordinates4D& p2, Coordinates4D& p3, FrameBuffer* fbo) {
//void triangleFill(Coordinates4D* vertices, ColorRGBA* colors, FrameBuffer* fbo) {
void triangleFill(FragmentInfo* fragments, FrameBuffer* fbo) {
	ColorRGBA color1, color2;
	
	FragmentInfo fragment;
	
	std::vector<FragmentInfo> f(fragments, fragments+3);
	
//	for (int i = 0; i < 3; i++) {
//		fragment.location3D = vertices[i];
//		fragment.color = colors[i];
//		f.push_back(fragment);
//	}
	std::sort(f.begin(), f.end(), compareCoordinatesFrag);

//	std::vector<Coordinates4D> v(vertices, vertices+3);
//	std::sort(v.begin(), v.end(), compareCoordinates);

	
//	drawHorizonalLineRGBA(v[0].x, v[0].x, v[0].y, color, fbo);
//	setFrameBufferRGBA(v[0].x, v[0].y, fbo, color);
//	setFrameBufferRGBA(v[1].x, v[1].y, fbo, color);
//	setFrameBufferRGBA(v[2].x, v[2].y, fbo, color);
	
	double slope12 = (f[1].location3D.x - f[0].location3D.x)/(f[1].location3D.y-f[0].location3D.y);
	if((int)f[1].location3D.y == (int)f[0].location3D.y) {
		slope12 = 0;
	}
	double slope13 = (f[2].location3D.x - f[0].location3D.x)/(f[2].location3D.y-f[0].location3D.y);
	
	double y = f[0].location3D.y;
	double yp = 0;// = y - v[0].y; <- add v[0].y to get y
	double x1, x2;
	
//	color.r = 0;
	
	// y goes from v[0].y to v[2].y;
	// yp goes from 0 to v[2].y - v[0].y
	// cf2 goes from 0 to 1, or along yp/(v[2].y - v[0].y)
	double colorFactor = 1.0/(f[2].location3D.y - f[0].location3D.y);
	
	// cf goes from 0 to 1, or along yp/(v[1].y - v[0].y)
	double colorFactor2 = 1.0/(f[1].location3D.y - f[0].location3D.y);
	
	
	for( ; (yp+f[0].location3D.y) < f[1].location3D.y; yp += 1.0) {
		
//		x1 = slope12*(y-v[0].y) + v[0].x;
//		x2 = slope13*(y-v[0].y) + v[0].x;
		x1 = slope12*(yp) + f[0].location3D.x;
		x2 = slope13*(yp) + f[0].location3D.x;
//		drawHorizonalLineRGBA(x1, x2, y, color, fbo);
		
//		color.r = color[0].r + colorFactor*(color[2].r - color[0].r);
		color2 = interpolate(f[0].color, f[2].color, colorFactor*yp);
		color1 = interpolate(f[0].color, f[1].color, colorFactor2*yp);
		
		drawHorizonalLineRGBA(x1, x2, yp+f[0].location3D.y, color1, color2, fbo);
	}
	yp =f[1].location3D.y - f[0].location3D.y;
	
	double slope23 = (f[2].location3D.x - f[1].location3D.x)/(f[2].location3D.y-f[1].location3D.y);
	if((int)f[2].location3D.y == (int)f[1].location3D.y) {
		slope23 = 0;
	}
	
	// now yp goes from v[1].y-v[0].y to v[2].y - v[0].y
	// cf2 goes from 0 to 1, or along (yp-(v[1].y-v[0].y))/((v[2].y - v[0].y) - (v[1].y-v[0].y))
//	colorFactor2 = 1.0/(f[2].location3D.y - (f[1].location3D.y-f[0].location3D.y));
	colorFactor2 = 1.0/((f[2].location3D.y-f[0].location3D.y) - (f[1].location3D.y-f[0].location3D.y));
	
	double end = f[2].location3D.y;
	for ( ; (yp+f[0].location3D.y) < end; yp += 1.0) {
		
//		x1 = slope23*(yp+v[0].y-v[1].y) + v[1].x;
//		x2 = slope13*(yp+v[0].y-v[0].y) + v[0].x;
		x1 = slope23*(yp+f[0].location3D.y-f[1].location3D.y) + f[1].location3D.x;
		x2 = slope13*(yp) + f[0].location3D.x;

		
		color2 = interpolate(f[0].color, f[2].color, colorFactor*yp);
		color1 = interpolate(f[1].color, f[2].color, colorFactor2*(yp-(f[1].location3D.y-f[0].location3D.y)));
		drawHorizonalLineRGBA(x1, x2, yp+f[0].location3D.y, color1, color2, fbo);
	}
	
}


//void drawPolygonWithTriangles( Polygon4D& poly, DepthBuffer* depthBuffer, char fill) {
void drawPolygonWithTriangles( Polygon4D& poly, FrameBuffer* fbo) {
//void drawPolygonWithTriangles( Polygon4D& poly, FrameBuffer* fbo) {
	if (poly.numVertices < 3) {
//		mvprintw(line++, 0, "Not enough Polygon vertices: %d", poly.numVertices);
		return;
	}
//	int minIndex = 0, maxIndex = 0;
//	double minY = poly.vertices[0].y;
//	double maxY = poly.vertices[0].y;
//	for (int i = 1; i < poly.numVertices; i++) {
//		if (poly.vertices[i].y < minY) {
//			minY = poly.vertices[i].y;
//			minIndex = i;
//		}
//		if (poly.vertices[i].y > maxY) {
//			maxY = poly.vertices[i].y;
//			maxIndex = i;
//		}
//	}
	
	FragmentInfo fragments[3];
	fragments[0].location3D = poly.vertices[0];
	fragments[0].color = poly.colors[0];
	for (int i = 2; i < poly.numVertices; i++) {
		fragments[1].location3D = poly.vertices[i-1];
		fragments[2].location3D = poly.vertices[i];
		fragments[1].color = poly.colors[i-1];
		fragments[2].color = poly.colors[i];
		triangleFill(fragments, fbo);
	}
	
	
}
