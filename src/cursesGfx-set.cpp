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
	
	if ( rgb.x < cMin ) {
		cMin = rgb.x;
	}
	if ( rgb.y < cMin ) {
		cMin = rgb.y;
	}
	if ( rgb.z < cMin ) {
		cMin = rgb.z;
	}
	
	if ( rgb.x > result.z ) {
		result.z = rgb.x;
	}
	if ( rgb.y > result.z) {
		result.z = rgb.y;
	}
	if ( rgb.z > result.z ) {
		result.z = rgb.z;
	}
	
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
    Coordinates3D hsv = hsvIn;
    if (hsv.y == 0)
    {
        r = hsv.z;
        g = hsv.z;
        b = hsv.z;
    }
    else
    {
        int i;
        double f, p, q, t;

        hsv.x = fmod(hsv.x, 360);
//        if (hsv.x == 360)
//            hsv.x = 0;
//        else
        hsv.x = hsv.x / 60;

        i = floor(hsv.x);// (int)trunc(hsv.x);
        f = hsv.x - i;

        p = hsv.z * (1.0 - hsv.y);
        q = hsv.z * (1.0 - (hsv.y * f));
        t = hsv.z * (1.0 - (hsv.y * (1.0 - f)));

        switch (i)
        {
        case 0:
            r = hsv.z;
            g = t;
            b = p;
            break;

        case 1:
            r = q;
            g = hsv.z;
            b = p;
            break;

        case 2:
            r = p;
            g = hsv.z;
            b = t;
            break;

        case 3:
            r = p;
            g = q;
            b = hsv.z;
            break;

        case 4:
            r = t;
            g = p;
            b = hsv.z;
            break;

        default:
            r = hsv.z;
            g = p;
            b = q;
            break;
        }

    }

    Coordinates3D rgb;
    rgb.x = r * 255;
    rgb.y = g * 255;
    rgb.z = b * 255;

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

//struct AsciiLutInfo {
//    char c;
//    uint8_t type;
//};
//#ifdef __APPLE__
//
//struct AsciiLutInfo asciiLUT[256] = {{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'`', 2},{'.', 2},{'.', 2},{'.', 2},{'.', 2},{'.', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'-', 2},{'\'', 2},{'\'', 2},{'\'', 2},{'\'', 2},{'\'', 2},{':', 2},{':', 2},{':', 2},{',', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'_', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'"', 2},{'~', 2},{'~', 2},{'~', 2},{'~', 2},{';', 2},{';', 2},{';', 2},{';', 2},{';', 2},{';', 2},{'|', 2},{'|', 2},{'|', 2},{'|', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'!', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'+', 2},{'/', 2},{'/', 2},{'\\', 2},{'\\', 2},{'^', 2},{'=', 2},{'i', 2},{'i', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'r', 2},{'l', 2},{'l', 2},{'l', 2},{'l', 2},{'L', 2},{'L', 2},{'c', 2},{'*', 2},{'*', 2},{'*', 2},{'v', 2},{'v', 2},{'v', 2},{'J', 2},{'J', 2},{'J', 2},{'T', 2},{'T', 2},{'7', 2},{'t', 2},{'t', 2},{'f', 2},{'1', 2},{'j', 2},{'j', 2},{'Y', 2},{'s', 2},{'u', 2},{'n', 2},{'(', 2},{'I', 2},{'C', 2},{'C', 2},{'}', 2},{'}', 2},{'{', 2},{'o', 2},{'o', 2},{'o', 2},{'e', 2},{'e', 2},{']', 2},{'y', 2},{'2', 2},{'2', 2},{'2', 2},{'a', 2},{'a', 2},{'S', 2},{'E', 2},{'h', 2},{'k', 2},{'V', 2},{'Z', 2},{'Z', 2},{'Z', 2},{'5', 2},{'5', 2},{'5', 2},{'5', 2},{'U', 2},{'U', 2},{'U', 2},{'p', 2},{'G', 2},{'q', 2},{'w', 2},{'d', 2},{'d', 2},{'A', 2},{'A', 2},{'H', 2},{'H', 2},{'K', 2},{'K', 2},{'K', 2},{'O', 2},{'O', 2},{'O', 2},{'D', 2},{'#', 2},{'9', 2},{'6', 2},{'$', 2},{'$', 2},{'$', 2},{'R', 2},{'R', 2},{'R', 2},{'R', 2},{'B', 2},{'B', 2},{'B', 2},{'g', 2},{'8', 2},{'8', 2},{'8', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'%', 2},{'%', 2},{'%', 2},{'%', 2},{'@', 2},{'0', 2},{'0', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2}};
//
//#else
//
//struct AsciiLutInfo asciiLUT[256] = {{' ', 1},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{' ', 2},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'`', 1},{'.', 1},{'.', 1},{'.', 1},{'-', 1},{'-', 1},{'-', 1},{'\'', 1},{'\'', 1},{':', 1},{',', 1},{'_', 1},{'`', 2},{'`', 2},{'.', 2},{'"', 1},{'"', 1},{';', 1},{';', 1},{';', 1},{'|', 1},{'|', 1},{'!', 1},{'!', 1},{'-', 2},{'+', 1},{'\\', 1},{'^', 1},{'=', 1},{'i', 1},{'\'', 2},{'r', 1},{'r', 1},{'?', 1},{'l', 1},{'L', 1},{'_', 2},{'v', 1},{'J', 1},{'J', 1},{'T', 1},{'f', 1},{'j', 1},{'x', 1},{'z', 1},{'C', 1},{'I', 1},{'o', 1},{'~', 2},{'y', 1},{'2', 1},{'[', 1},{'a', 1},{'3', 1},{'V', 1},{'P', 1},{';', 2},{'5', 1},{'U', 1},{'G', 1},{'m', 1},{'d', 1},{'A', 1},{'K', 1},{'!', 2},{'O', 1},{'D', 1},{'$', 1},{'$', 1},{'R', 1},{'g', 1},{'+', 2},{'8', 1},{'8', 1},{'N', 1},{'/', 2},{'Q', 1},{'&', 1},{'@', 1},{'=', 2},{'i', 2},{'i', 2},{'0', 1},{'r', 2},{'M', 1},{'M', 1},{'W', 1},{'W', 1},{'W', 1},{'W', 1},{'l', 2},{'l', 2},{'L', 2},{'L', 2},{'c', 2},{'*', 2},{'*', 2},{'*', 2},{'v', 2},{'v', 2},{'v', 2},{'J', 2},{'J', 2},{'J', 2},{'T', 2},{'T', 2},{'7', 2},{'t', 2},{'t', 2},{'f', 2},{'1', 2},{'j', 2},{'j', 2},{'Y', 2},{'F', 2},{'u', 2},{'n', 2},{'(', 2},{'I', 2},{'C', 2},{'C', 2},{'}', 2},{'}', 2},{'{', 2},{'o', 2},{'o', 2},{'o', 2},{'e', 2},{'e', 2},{']', 2},{'y', 2},{'2', 2},{'2', 2},{'2', 2},{'a', 2},{'a', 2},{'S', 2},{'E', 2},{'h', 2},{'k', 2},{'V', 2},{'Z', 2},{'Z', 2},{'Z', 2},{'5', 2},{'5', 2},{'5', 2},{'5', 2},{'U', 2},{'U', 2},{'U', 2},{'p', 2},{'G', 2},{'q', 2},{'w', 2},{'d', 2},{'d', 2},{'A', 2},{'A', 2},{'H', 2},{'H', 2},{'K', 2},{'K', 2},{'K', 2},{'O', 2},{'O', 2},{'O', 2},{'D', 2},{'#', 2},{'9', 2},{'6', 2},{'$', 2},{'$', 2},{'$', 2},{'R', 2},{'R', 2},{'R', 2},{'R', 2},{'B', 2},{'B', 2},{'B', 2},{'g', 2},{'8', 2},{'8', 2},{'8', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'N', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'Q', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'&', 2},{'%', 2},{'%', 2},{'%', 2},{'%', 2},{'@', 2},{'0', 2},{'0', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'M', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2},{'W', 2}};
//#endif //__APPLE__
//
//void setRGB( const Coordinates2D& pixel, const Coordinates3D& rgb) {
//
//	Coordinates3D clippedRGB = clipRGB(rgb);
//	Coordinates3D hsl = rgbToHsv(clippedRGB);
//
////	int hueIndex = floor(hsl.x + 1.5);
////    if(hueIndex > 6) {
////        hueIndex = 1;
////    }
////
////	if (hsl.y < 0.33) {
////		hueIndex = 7;
////	}
//    int colorIndex = floor(hsl.x/6.0*31.0);
//    colorIndex |= (((uint8_t)(hsl.y*7)) << 5) & 0xE0;
////    colorIndex |= (((uint8_t)(hsl.y*8)) << 5) & 0xE0;
//
//    attron(COLOR_PAIR(colorIndex));
//    struct AsciiLutInfo aLut = asciiLUT[(int)(hsl.z*255)];
//    if (aLut.type == 1) {
//        attron(WA_DIM);
////        attrset(WA_DIM);
//        set(pixel, aLut.c);
//        attroff(WA_DIM);
////        move(pixel.y, pixel.x);
////        putchar(aLut.c);
////        fprintf(stderr, "%c", aLut.c);
////        move(pixel.y, pixel.x);
////        mvprintw(pixel.x, pixel.y, "\x1b[2;37;40m%c\x1b[0m", aLut.c);
//    } else {
//        set(pixel, aLut.c);
//    }
//
//    attroff(COLOR_PAIR(colorIndex));
//    return;
//    /*
//	attron(COLOR_PAIR(hueIndex));
////    attron(COLOR_PAIR(1));
//
//
//
//	char c = ' ';
//	const double maxL = 1.0;
//#ifdef ASCII_EXTENDED
//	const double numLevels = 20;
//#else
//	const double numLevels = 16;
//#endif
//	double step = numLevels - 1.0;
//	double LperLs = maxL/numLevels;
//
//#ifdef ASCII_EXTENDED
//	if (hsl.z > step-- * LperLs) {
////		c = '█';
//		mvprintw(pixel.y, pixel.x, "▓");	//dark
//		attroff(COLOR_PAIR(hueIndex));
//		return;
//	} else if (hsl.z > step-- * LperLs) {
//		mvprintw(pixel.y, pixel.x, "█");	// full
//		attroff(COLOR_PAIR(hueIndex));
//		return;
//	} else
//#endif
//		if (hsl.z > step-- * LperLs) {
//		c = '$';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '@';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '%';
//
//#ifdef ASCII_EXTENDED
//	} else if (hsl.z > step-- * LperLs) {
//
//
//		mvprintw(pixel.y, pixel.x, "▒"); // medium "▒" ▓ █
//		attroff(COLOR_PAIR(hueIndex));
//		return;
//#endif
////	}  else if (hsl.z > 0.50) {
////		c = '*';
//
//	} else if (hsl.z > step-- * LperLs) {
//		c = '#';
//
//	}	else if (hsl.z > step-- * LperLs) {
//		c = 'o';
//#ifdef ASCII_EXTENDED
//	} else if (hsl.z > step-- * LperLs) {
//
//
//		mvprintw(pixel.y, pixel.x, "░"); // light "▒" ▓ █
//		attroff(COLOR_PAIR(hueIndex));
//		return;
//#endif
//	} else if (hsl.z > step-- * LperLs) {
//		c = '=';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '*';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '+';
//	} else if (hsl.z > step-- * LperLs) {
//		c = ';';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '~';
////	} else if (hsl.z > maxL*6.0/numLevels) {
////		c = '"';
//	} else if (hsl.z > step-- * LperLs) {
//		c = ',';
//	} else if (hsl.z > step-- * LperLs) {
//		c = ':';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '-';
//	} else if (hsl.z > step-- * LperLs) {
//		c = '.';
//	} else if (hsl.z > step-- * LperLs) {
//
//		c = '`';
//	}
//
////	attron(COLOR_PAIR(4));
//	set(pixel, c);
//
//	attroff(COLOR_PAIR(hueIndex));
////    attroff(COLOR_PAIR(1));
//     */
//}


void setRGB( FrameBuffer* fbo, const Coordinates2D& pixel, const Coordinates3D& rgb) {
	int offset = pixel.y * fbo->cols + pixel.x;
	fbo->data[offset + 0] = rgb.x;
	fbo->data[offset + 1] = rgb.y;
	fbo->data[offset + 2] = rgb.z;
}
