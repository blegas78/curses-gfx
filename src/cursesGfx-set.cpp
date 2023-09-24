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
	
	
//	result.z = cMax;
	
	return result;
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


void setRGB( const Coordinates2D& pixel, const Coordinates3D& rgb) {
	
	Coordinates3D clippedRGB = clipRGB(rgb);
	Coordinates3D hsl = rgbToHsv(clippedRGB);
	
	int hueIndex = floor(hsl.x + 1);
	
	if (hsl.y < 0.33) {
		hueIndex = 7;
	}
	
	attron(COLOR_PAIR(hueIndex));
//    attron(COLOR_PAIR(1));
	

	
	char c = ' ';
	const double maxL = 1.0;
#ifdef ASCII_EXTENDED
	const double numLevels = 20;
#else
	const double numLevels = 16;
#endif
	double step = numLevels - 1.0;
	double LperLs = maxL/numLevels;
	
#ifdef ASCII_EXTENDED
	if (hsl.z > step-- * LperLs) {
//		c = '█';
		mvprintw(pixel.y, pixel.x, "▓");	//dark
		attroff(COLOR_PAIR(hueIndex));
		return;
	} else if (hsl.z > step-- * LperLs) {
		mvprintw(pixel.y, pixel.x, "█");	// full
		attroff(COLOR_PAIR(hueIndex));
		return;
	} else
#endif
		if (hsl.z > step-- * LperLs) {
		c = '$';
	} else if (hsl.z > step-- * LperLs) {
		c = '@';
	} else if (hsl.z > step-- * LperLs) {
		c = '%';
	
#ifdef ASCII_EXTENDED
	} else if (hsl.z > step-- * LperLs) {
		
		
		mvprintw(pixel.y, pixel.x, "▒"); // medium "▒" ▓ █
		attroff(COLOR_PAIR(hueIndex));
		return;
#endif
//	}  else if (hsl.z > 0.50) {
//		c = '*';
		
	} else if (hsl.z > step-- * LperLs) {
		c = '#';
		
	}	else if (hsl.z > step-- * LperLs) {
		c = 'o';
#ifdef ASCII_EXTENDED
	} else if (hsl.z > step-- * LperLs) {
		
		
		mvprintw(pixel.y, pixel.x, "░"); // light "▒" ▓ █
		attroff(COLOR_PAIR(hueIndex));
		return;
#endif
	} else if (hsl.z > step-- * LperLs) {
		c = '=';
	} else if (hsl.z > step-- * LperLs) {
		c = '*';
	} else if (hsl.z > step-- * LperLs) {
		c = '+';
	} else if (hsl.z > step-- * LperLs) {
		c = ';';
	} else if (hsl.z > step-- * LperLs) {
		c = '~';
//	} else if (hsl.z > maxL*6.0/numLevels) {
//		c = '"';
	} else if (hsl.z > step-- * LperLs) {
		c = ',';
	} else if (hsl.z > step-- * LperLs) {
		c = ':';
	} else if (hsl.z > step-- * LperLs) {
		c = '-';
	} else if (hsl.z > step-- * LperLs) {
		c = '.';
	} else if (hsl.z > step-- * LperLs) {
		
		c = '`';
	}
	
//	attron(COLOR_PAIR(4));
	set(pixel, c);
	
	attroff(COLOR_PAIR(hueIndex));
//    attroff(COLOR_PAIR(1));
}


void setRGB( FrameBuffer* fbo, const Coordinates2D& pixel, const Coordinates3D& rgb) {
	int offset = pixel.y * fbo->cols + pixel.x;
	fbo->data[offset + 0] = rgb.x;
	fbo->data[offset + 1] = rgb.y;
	fbo->data[offset + 2] = rgb.z;
}
