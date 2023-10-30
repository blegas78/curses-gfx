#include "curses-gfx-3d.h"
#include "curses-gfx.h"

#include <cmath>
//#include <unistd.h>
//#include <stdint.h>

#include <ncurses.h>


Coordinates3D rgbToHsv( const Coordinates3D& rgb) {
	Coordinates3D result;
	
	double cMin = 1e100, delta;
	
	result.z = -1e100;
	
    cMin = fmin(fmin(rgb.x, rgb.y) ,rgb.z);
    
//	if ( rgb.x < cMin ) {
//		cMin = rgb.x;
//	}
//	if ( rgb.y < cMin ) {
//		cMin = rgb.y;
//	}
//	if ( rgb.z < cMin ) {
//		cMin = rgb.z;
//	}
//
//	if ( rgb.x > result.z ) {
//		result.z = rgb.x;
//	}
//	if ( rgb.y > result.z) {
//		result.z = rgb.y;
//	}
//	if ( rgb.z > result.z ) {
//		result.z = rgb.z;
//	}
    result.z = fmax(fmax(rgb.x, rgb.y), rgb.z);
	
	delta = result.z - cMin;
	
	
	if(delta == 0) {
		result.x = 0;
	} else if (result.z == rgb.x) {
		result.x = ((rgb.y - rgb.z)/delta + (rgb.y < rgb.z ? 6 : 0));
	} else if (result.z == rgb.y) {
		result.x = ((rgb.z - rgb.x)/delta + 2);
	} else if (result.z == rgb.z) {
		result.x = ((rgb.x - rgb.y)/delta + 4);
	}
	if (result.z == 0) {
		result.y = 0;
	} else {
		result.y = delta/result.z;
	}
	
    result.x *= 60.0;
//	result.z = cMax;
	
	return result;
}

Coordinates3D hsvToRgb(const Coordinates3D& hsvIn) {
    double r = 0, g = 0, b = 0;
//    Coordinates3D hsv = hsvIn;
    if (hsvIn.y == 0)
    {
        r = hsvIn.z*255;
//        g = hsv.z;
//        b = hsv.z;
        return {r, r, r};
    }
//    else
//    {
//        int i;
//        double f, p, q, t;
//
//        hsv.x = fmod(hsv.x, 360);
////        if (hsv.x == 360)
////            hsv.x = 0;
////        else
//        hsv.x = hsv.x / 60;
//
//        i = floor(hsv.x);// (int)trunc(hsv.x);
//        f = hsv.x - i;
//
//        p = hsv.z * (1.0 - hsv.y);
//        q = hsv.z * (1.0 - (hsv.y * f));
//        t = hsv.z * (1.0 - (hsv.y * (1.0 - f)));
//
//        switch (i)
//        {
//        case 0:
//            r = hsv.z;
//            g = t;
//            b = p;
//            break;
//
//        case 1:
//            r = q;
//            g = hsv.z;
//            b = p;
//            break;
//
//        case 2:
//            r = p;
//            g = hsv.z;
//            b = t;
//            break;
//
//        case 3:
//            r = p;
//            g = q;
//            b = hsv.z;
//            break;
//
//        case 4:
//            r = t;
//            g = p;
//            b = hsv.z;
//            break;
//
//        default:
//            r = hsv.z;
//            g = p;
//            b = q;
//            break;
//        }
//
//    }

    Coordinates3D rgb;
//    rgb.x = r * 255;
//    rgb.y = g * 255;
//    rgb.z = b * 255;
    
    double C = hsvIn.y*hsvIn.z;
    double X = C*(1 - fabs(fmod(hsvIn.x/60.0, 2) - 1));
    double m = hsvIn.z - C;

    if(hsvIn.x < 60) {
        r = C; g = X; b = 0;
    } else if(hsvIn.x < 120) {
        r = X; g = C; b = 0;
    } else if(hsvIn.x < 180) {
        r = 0; g = C; b = X;
    } else if(hsvIn.x < 240) {
        r = 0; g = X; b = C;
    } else if(hsvIn.x < 300) {
        r = X; g = 0; b = C;
    } else { //} if(hsvIn.x < 360) {
        r = C; g = 0; b = X;
    }

    rgb.x = (r+m) * 255;
    rgb.y = (g+m) * 255;
    rgb.z = (b+m) * 255;
    
    return rgb;
}

// from https://www.programmingalgorithms.com/algorithm/hsl-to-rgb/c/
float hueToRgb(float v1, float v2, float vH)
{
	if (vH < 0)
		vH += 1;

	if (vH > 1)
		vH -= 1;

	if ((6 * vH) < 1)
		return (v1 + (v2 - v1) * 6 * vH);

	if ((2 * vH) < 1)
		return v2;

	if ((3 * vH) < 2)
		return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);

	return v1;
}

Coordinates3D hslToRgb(Coordinates3D hsl) {
	Coordinates3D rgb;

	if (hsl.y == 0)
	{
		rgb.x = rgb.y = rgb.z = (unsigned char)(hsl.z * 255);
	}
	else
	{
		float v1, v2;
		float hue = (float)hsl.x / 360;

		v2 = (hsl.z < 0.5) ? (hsl.z * (1 + hsl.y)) : ((hsl.z + hsl.y) - (hsl.z * hsl.y));
		v1 = 2 * hsl.z - v2;

		rgb.x = (unsigned char)(255 * hueToRgb(v1, v2, hue + (1.0f / 3)));
		rgb.y = (unsigned char)(255 * hueToRgb(v1, v2, hue));
		rgb.z = (unsigned char)(255 * hueToRgb(v1, v2, hue - (1.0f / 3)));
	}

	return rgb;
}

Coordinates3D rgbToHsl( const Coordinates3D& rgb) {
	Coordinates3D result;
	
	double cMax = -1e100, cMin = 1e100, delta;
	
	if ( rgb.x < cMin ) {
		cMin = rgb.x;
	}
	if ( rgb.y < cMin ) {
		cMin = rgb.y;
	}
	if ( rgb.z < cMin ) {
		cMin = rgb.z;
	}
	
	if ( rgb.x > cMax ) {
		cMax = rgb.x;
	}
	if ( rgb.y > cMax ) {
		cMax = rgb.y;
	}
	if ( rgb.z > cMax ) {
		cMax = rgb.z;
	}
	
	delta = cMax - cMin;
	
	result.z = (cMax + cMin)/2.0;
	
	if(delta == 0) {
		result.x = 0;
		result.y = 0;
	} else if (cMax == rgb.x) {
		result.x = ((rgb.y - rgb.z)/delta + (rgb.y < rgb.z ? 6 : 0));
		result.y = delta/(1 - fabs(2*result.z - 1));
	} else if (cMax == rgb.y) {
		result.x = ((rgb.z - rgb.x)/delta + 2);
		result.y = delta/(1 - fabs(2*result.z - 1));
	} else if (cMax == rgb.z) {
		result.x = ((rgb.x - rgb.y)/delta + 4);
		result.y = delta/(1 - fabs(2*result.z - 1));
	}
	
	return result;
}

Coordinates3D clipRGB(Coordinates3D rgb) {
	
	if (rgb.x > 1) {
		rgb.x = 1;
	} else if (rgb.x < 0) {
		rgb.x = 0;
	}
	if (rgb.y > 1) {
		rgb.y = 1;
	} else if (rgb.y < 0) {
		rgb.y = 0;
	}
	if (rgb.z > 1) {
		rgb.z = 1;
	} else if (rgb.z < 0) {
		rgb.z = 0;
	}
	
	
	return rgb;
}


void setRGB( FrameBuffer* fbo, const Coordinates2D& pixel, const Coordinates3D& rgb) {
	int offset = pixel.y * fbo->cols + pixel.x;
	fbo->data[offset + 0] = rgb.x;
	fbo->data[offset + 1] = rgb.y;
	fbo->data[offset + 2] = rgb.z;
}
