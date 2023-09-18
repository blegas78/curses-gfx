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
	

	
	char c = ' ';
	const double maxL = 1.0;
#ifdef ASCII_EXTENDED
	const double numLevels = 20;
#else
	const double numLevels = 16;
#endif
	double step = numLevels - 1.0;
	
#ifdef ASCII_EXTENDED
	if (hsl.z > maxL* step--/numLevels) {
//		c = '█';
		mvprintw(pixel.y, pixel.x, "▓");	//dark
		attroff(COLOR_PAIR(hueIndex));
		return;
	} else if (hsl.z > maxL * step--/numLevels) {
		mvprintw(pixel.y, pixel.x, "█");	// full
		attroff(COLOR_PAIR(hueIndex));
		return;
	} else
#endif
		if (hsl.z > maxL*step--/numLevels) {
		c = '$';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '@';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '%';
	
#ifdef ASCII_EXTENDED
	} else if (hsl.z > maxL*step--/numLevels) {
		
		
		mvprintw(pixel.y, pixel.x, "▒"); // medium "▒" ▓ █
		attroff(COLOR_PAIR(hueIndex));
		return;
#endif
//	}  else if (hsl.z > 0.50) {
//		c = '*';
		
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '#';
		
	}	else if (hsl.z > maxL*step--/numLevels) {
		c = 'o';
#ifdef ASCII_EXTENDED
	} else if (hsl.z > maxL*step--/numLevels) {
		
		
		mvprintw(pixel.y, pixel.x, "░"); // light "▒" ▓ █
		attroff(COLOR_PAIR(hueIndex));
		return;
#endif
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '=';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '*';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '+';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = ';';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '~';
//	} else if (hsl.z > maxL*6.0/numLevels) {
//		c = '"';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = ',';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = ':';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '-';
	} else if (hsl.z > maxL*step--/numLevels) {
		c = '.';
	} else if (hsl.z > maxL*step--/numLevels) {
		
		c = '`';
	}
	
//	attron(COLOR_PAIR(4));
	set(pixel, c);
	
	attroff(COLOR_PAIR(hueIndex));
}


void setRGB( FrameBuffer* fbo, const Coordinates2D& pixel, const Coordinates3D& rgb) {
	int offset = pixel.y * fbo->cols + pixel.x;
	fbo->data[offset + 0] = rgb.x;
	fbo->data[offset + 1] = rgb.y;
	fbo->data[offset + 2] = rgb.z;
}
