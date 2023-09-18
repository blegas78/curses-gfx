#include "curses-gfx-handler.h"

#include "curses-gfx-3d.h"

#include <ncurses.h>
RenderPipeline::RenderPipeline() {
	fbo = new FrameBuffer;
	depthBuffer = new FrameBuffer;
	
	asTexImage2d(fbo, FBT_RGBA, 1, 1);
	asTexImage2d(depthBuffer, FBT_DEPTH, 1, 1);
	depthClearColor = 00;
#ifdef FB_SUPPORT
	setupLinuxFb();
#endif
}

#ifdef FB_SUPPORT
void RenderPipeline::setupLinuxFb() {
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
	
	fb_data_size = fb_width * fb_height * fb_bytes*64;

	fbdata = (char*)mmap (0, fb_data_size, PROT_READ | PROT_WRITE, MAP_SHARED, fbfd, (off_t)0);
	
	printf("Framebuffer info: res %dx%d, bpp=%d, Bytes=%d\n", fb_width, fb_height, fb_bpp, fb_bytes);
}
#endif

void saveFrameBufferToFile(const char* filename, FrameBuffer* fbo) {
	FILE *fp = fopen(filename, "wb"); /* b - binary mode */
	switch (fbo->type) {
		case FBT_RGBA:
			fprintf(fp, "P6\n%d %d\n255\n", fbo->cols, fbo->rows);
			for (int i = 0; i < fbo->cols*fbo->rows; i++) {
				fwrite( &((ColorRGBA*)fbo->data)[i], 1, 3, fp);
			}
			break;
			
		case FBT_DEPTH:
			fprintf(fp, "P5\n%d %d\n255\n", fbo->cols, fbo->rows);
			for (int i = 0; i < fbo->cols*fbo->rows; i++) {
				uint8_t data;
				if ( ((double*)fbo->data)[i] == 0 ) {
					data = 0;
				} else {
					data = (((double*)fbo->data)[i]-1)*200000;
				}
//				printf("%f == %d\n", ((double*)fbo->data)[i], data);
				fwrite( &data, 1, 1, fp);
			}
			
		default:
			break;
	}
	
	fclose(fp);
}

RenderPipeline::~RenderPipeline() {
	// Fuck it, save the render and depth buffer to a file:

	saveFrameBufferToFile("render.ppm", &fbo[0]);
	saveFrameBufferToFile("depth.pgm", &depthBuffer[0]);
	
	delete fbo;
	delete depthBuffer;
}


void RenderPipeline::resize(int width, int height) {
	asTexImage2d(fbo, FBT_RGBA, width, height);
	asTexImage2d(depthBuffer, FBT_DEPTH, width, height);
	
//	tempD.setSize(width, height);
}
void RenderPipeline::reset() {
	depthBuffer->clear(&depthClearColor);
	
	ColorRGBA clearColor;
	clearColor.r = 0;
	clearColor.g = 0;
	clearColor.b = 0;
	clearColor.a = 0;
	fbo->clear(&clearColor);
}

void RenderPipeline::setFragmentShader(void (*fragmentShader)(const FragmentInfo&)) {
	this->fragmentShader = fragmentShader;
};

void RenderPipeline::rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line) {
	
	
//	void rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line) {
		Polygon4D quadAsPolygon;
		for (int i = 0; i < count; i++) {
			
			quadAsPolygon.vertices[0] = vertices[quadIndices[i][0]];
			quadAsPolygon.vertices[1] = vertices[quadIndices[i][1]];
			quadAsPolygon.vertices[2] = vertices[quadIndices[i][2]];
			quadAsPolygon.vertices[3] = vertices[quadIndices[i][3]];
			quadAsPolygon.numVertices = 4;
			
			fillPolygonNormals(&quadAsPolygon, 1);
			
			rasterizePolygonsShader(&quadAsPolygon, 1, modelView, projection, viewport, userData, line);
			
		}
//	}
}
void RenderPipeline::rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line) {
	
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
		
		// now get back the original coordinates or separate handling
		polygonRestored.numVertices = polygonProjected.numVertices;
		for (int v = 0; v < polygonRestored.numVertices; v++) {
			polygonRestored.vertices[v] = matrixVectorMultiply(inverseProjection, polygonProjected.vertices[v]);
			polygonRestored.normals[v] = polygonProjected.normals[v];
		}
		
		// // Convert homogeneous
		for (int v = 0; v < polygonProjected.numVertices; v++) {
			if( polygonProjected.vertices[v].w != 1) {
				polygonProjected.vertices[v].x = polygonProjected.vertices[v].x / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].y = -polygonProjected.vertices[v].y / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].z = polygonProjected.vertices[v].z / polygonProjected.vertices[v].w;
				polygonProjected.vertices[v].w = 1;
			}
		}
		
		// Viewport
		for (int v = 0; v < polygonProjected.numVertices; v++) {
			polygonProjected.vertices[v] = matrixVectorMultiply(viewport, polygonProjected.vertices[v]);
		}
		
		// Draw:
//		mvprintw(line++, 0, "Calling drawPolygon() with %d vertices", polygon.numVertices);
//		refresh();
//		drawPolygon( polygon, depthBuffer, fill, line);
//		drawPolygonShader( polygonProjected, polygonRestored, userData, depthBuffer, fragmentShader, line);
		drawPolygonShader( polygonProjected, polygonRestored, userData, line);
//		mvprintw(line++, 0, "Calling drawPolygon() Done");
//		refresh();
		
	}
}

void RenderPipeline::drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, int &line) {

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
	
	
//	char c = 'A';
	for ( int i = minIndex; ; ) {
//		mvprintw(line++, 0, "%d %c %f,%f %s", i, c, poly.vertices[i].x, poly.vertices[i].y, i == minIndex ? "Min" : (i == maxIndex ? "Max" : ""));
//		setWithDepthBuffer(onlyXY(poly.vertices[i]), 'o', poly.vertices[i].z, depthBuffer);
//		setFloatDotWithDepthBuffer(poly.vertices[i].x+0.5, poly.vertices[i].y+0.5, poly.vertices[i].z, depthBuffer);
		setFloatDotWithDepthBuffer(poly.vertices[i].x+0.5, poly.vertices[i].y+0.5, poly.vertices[i].z);
//		refresh();
		i = mod(i + 1, poly.numVertices);
		if (i == minIndex)
			break;
	}
	
	
	int indexRight = minIndex;
	int indexLeft = minIndex;
	
	Coordinates2D ar = onlyXY(poly.vertices[minIndex]);
	Coordinates2D br = onlyXY(poly.vertices[mod(minIndex+1, poly.numVertices)]);
	Coordinates3D arNormal = restored.normals[minIndex];
	Coordinates4D ar3dPoint = restored.vertices[minIndex];
	Coordinates3D brNormal = restored.normals[mod(minIndex+1, poly.numVertices)];
	Coordinates4D br3dPoint = restored.vertices[mod(minIndex+1, poly.numVertices)];
	double arDepth = poly.vertices[minIndex].z;
	double brDepth = poly.vertices[mod(minIndex+1, poly.numVertices)].z;
	
	Coordinates2D al = onlyXY(poly.vertices[minIndex]);
	Coordinates2D bl = onlyXY(poly.vertices[mod(minIndex-1, poly.numVertices)]);
	Coordinates3D alNormal = restored.normals[minIndex];
	Coordinates4D al3dPoint = restored.vertices[minIndex];
	Coordinates3D blNormal = restored.normals[mod(minIndex-1, poly.numVertices)];
	Coordinates4D bl3dPoint = restored.vertices[mod(minIndex-1, poly.numVertices)];
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
//	int lineY[3] = {0,0,0};// = a.y;
	double lineDepthStart;//[3] = {0,0,0};
	double lineDepthEnd;
	
	FragmentInfo fragStart, fragEnd;
	
	
//	Coordinates4D point3dr;
//	Coordinates4D point3dl;
//	Coordinates3D normalr;
//	Coordinates3D normall;
	
	bool skipRightLine = false;
	bool skipLeftLine = true;
	
	bool rightComplete = false;
	bool leftComplete = false;
	
	int lineCount = 0;
	int tr = 0;
	int tl = 0;
	
	double savearDepth, savebrDepth, savedMagnitudeR = 1;
	Coordinates2D savedBr;
	Coordinates4D savedBr3dPoint, savedAr3dPoint;
	
	double savealDepth, saveblDepth, savedMagnitude = 1;
	Coordinates2D savedBl;
	Coordinates4D savedBl3dPoint, savedAl3dPoint;
	
	
//	while (indexRight != maxIndex || indexLeft != maxIndex) {
	while (!leftComplete || !rightComplete) {
//		refresh();
//		usleep(500000);
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
//				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);

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
			depthsr[2] = 1.0/(1.0/brDepth + alphar*(1.0/arDepth - 1.0/brDepth));
//			lineDepthStart = depthsr[2];
//
//			normalr = interpolate(brNormal, arNormal, alphar);
//			point3dr = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, depthsr[2], alphar);
			

		
		
		if (ar.x == br.x && ar.y == br.y && indexRight != maxIndex) {
			// new line time
			
//			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);
			
			
			indexRight = mod(indexRight + 1, poly.numVertices);
			
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
			arNormal = restored.normals[indexRight];
			ar3dPoint = restored.vertices[indexRight];
			arDepth = poly.vertices[indexRight].z;
			
//			lineStartX[2] = ar.x;
			
			int nextIndex = mod(indexRight + 1, poly.numVertices);
			br = onlyXY(poly.vertices[nextIndex]);
			brNormal = restored.normals[nextIndex];
			br3dPoint = restored.vertices[nextIndex];
			brDepth = poly.vertices[nextIndex].z;
			
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
//				lineY[0] = lineY[1];
//				lineY[1] = lineY[2];
//				lineEndX[0] = lineEndX[1];
//				lineEndX[1] = lineEndX[2];
				
//				lineDepthStart[0] = lineDepthStart[1];
//				lineDepthStart[1] = lineDepthStart[2];
//				lineDepthEnd[0] = lineDepthEnd[1];
//				lineDepthEnd[1] = lineDepthEnd[2];
//				if (sx > 0 && sx2 > 0) {
//					lineStartX[2] = a.x + (a2.x < a.x ? 0 : 1);
//				} else {
				lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
				
//				savedMagnitudeR = lineMagnitudeSqr;
//				savedBr = br;
//				savearDepth = arDepth;
//				savebrDepth = brDepth;
//				savedAr3dPoint = ar3dPoint;
//				savedBr3dPoint = br3dPoint;
				
//				dxnr = lineStartX[1]-br.x;
//				dynr = (ar.y-1)-br.y;
//				alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr))/lineMagnitudeSqr;
//
//				lineDepthStart = 1.0/(1.0/brDepth + alphar*(1.0/arDepth - 1.0/brDepth));
////				normalr = interpolate(brNormal, arNormal, alphar);
////				point3dr = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, lineDepthStart, alphar);
//
//
//				fragStart.pixel.x = lineStartX[1];
//				fragStart.pixel.y = ar.y-1;
//				fragStart.normal = interpolate(brNormal, arNormal, alphar);
//				fragStart.location3D = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, lineDepthStart, alphar);
					
				
				dxnr = lineStartX[1]-savedBr.x;
				dynr = (ar.y-1)-savedBr.y;
				alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr))/savedMagnitudeR;
				
				lineDepthStart = 1.0/(1.0/savebrDepth + alphar*(1.0/savearDepth - 1.0/savebrDepth));
				
				
				fragStart.pixel.x = lineStartX[1];
				fragStart.pixel.y = ar.y-1;
				fragStart.normal = interpolate(brNormal, arNormal, alphar);
				fragStart.location3D = perspectiveInterpolate(savedBr3dPoint, savedAr3dPoint, savebrDepth, savearDepth, lineDepthStart, alphar);
				
				savedMagnitudeR = lineMagnitudeSqr;
				savedBr = br;
				savearDepth = arDepth;
				savebrDepth = brDepth;
				savedAr3dPoint = ar3dPoint;
				savedBr3dPoint = br3dPoint;

			} else {
				if (abs(al.x - ar.x) < abs(al.x - lineStartX[2])) {
//					mvaddstr(line++, 0, string);
					lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
//					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
					
					savedMagnitudeR = lineMagnitudeSqr;
					savedBr = br;
					savearDepth = arDepth;
					savebrDepth = brDepth;
					savedAr3dPoint = ar3dPoint;
					savedBr3dPoint = br3dPoint;
				}
				//				}
				errsr[2] -= dxr;
//				skipLeftLine = true;
			}
		}
		
		if (!skipLeftLine) {
			if (tl++ > 2) {
//				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1]);
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
//			depthsl[2] = blDepth - alphal*(blDepth - alDepth);
			depthsl[2] = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
//			lineDepthEnd[2] = depthsl[2];
			
			
//			normall[0] = normall[1];
//			normall[1] = normall[2];
//			normall[2] = interpolate(blNormal, alNormal, alphal);
			
//			point3dl[0] = point3dl[1];
//			point3dl[1] = point3dl[2];
////			point3dl[2] = interpolate(bl3dPoint, al3dPoint, alphal);
////
////			point3dl[2].z = 1.0/(1.0/bl3dPoint.z + alphal*(1.0/al3dPoint.z - 1.0/bl3dPoint.z));
//
//
//			point3dl[2] = perspectiveInterpolate(bl3dPoint, al3dPoint, blDepth, alDepth, depthsl[2], alphal);
//
//			if(al.x == bl.x && al.y == bl.y) {
//				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
//			}
//		}
		
		if (al.x == bl.x && al.y == bl.y && indexLeft != maxIndex) {
			// new line time
//			mvprintw(line++, 0, "New Left line!");
//
//			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1]);
			
			indexLeft = mod(indexLeft - 1, poly.numVertices);
			
			if (indexLeft == maxIndex) {
//				mvprintw(line++, 0, "Left line complete!");
//				refresh();
//				skipLeftLine = true;
//				skipRightLine = false;
				leftComplete = true;
				continue;
			}
//			mvprintw(line++, 0, "New Left line: %d - %d", indexLeft, nextIndex);
//			savealDepth = alDepth;
//			saveblDepth = blDepth;
			
//			int nextIndex = mod(indexLeft-1, poly.numVertices);
			al = onlyXY(poly.vertices[indexLeft]);
			alNormal = restored.normals[indexLeft];
			al3dPoint = restored.vertices[indexLeft];
			alDepth = poly.vertices[indexLeft].z;
			
//			if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
////						mvprintw(line++, 1, "Correcting lineEndX[1]");
//			lineEndX[1] = lineEndX[2] < al.x ? lineEndX[2] : al.x;
//			}
			
			int nextIndex = mod(indexLeft - 1, poly.numVertices);
			bl = onlyXY(poly.vertices[nextIndex]);
			blNormal = restored.normals[nextIndex];
			bl3dPoint = restored.vertices[nextIndex];
			blDepth = poly.vertices[nextIndex].z;
			
			
			
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
//				lineEndX[1] = lineEndX[2];
				
//				skipRightLine = false;
	
				
				if (lineCount++ > 0) {
					
					dxnl = lineEndX[1]-savedBl.x;
					dynl = (al.y-1)-savedBl.y;
					alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/savedMagnitude;
					lineDepthEnd = 1.0/(1.0/saveblDepth + alphal*(1.0/savealDepth - 1.0/saveblDepth));
					
					fragEnd.pixel.x = lineEndX[1];
					fragEnd.pixel.y = al.y-1;
					fragEnd.normal = interpolate(blNormal, alNormal, alphal);
					fragEnd.location3D = perspectiveInterpolate(savedBl3dPoint, savedAl3dPoint, saveblDepth, savealDepth, lineDepthEnd, alphal);
					
					
//						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
					//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//					drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], al.y-1, lineDepthStart, lineDepthEnd, point3dr,  point3dl, normalr, normall, modelView, userData, depthBuffer, fragmentShader);
					static bool shouldPlot = true;
					if (shouldPlot) {
//						shouldPlot = false;
						drawHorizonalLineWithShader(fragStart, fragEnd, lineDepthStart, lineDepthEnd, userData);
					} else {
						shouldPlot = true;
					}
					
					

				}
				
				lineEndX[1] = lineEndX[2];
				
				savedBl = bl;
				savedMagnitude = lineMagnitudeSql;
				savealDepth = alDepth;
				saveblDepth = blDepth;
				savedAl3dPoint = al3dPoint;
				savedBl3dPoint = bl3dPoint;
				
			} else {
				errsl[2] -= dxl;
//				skipRightLine = true;
				
				//if (!(ar.x == br.x && ar.y == br.y)) {
					if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
//						mvprintw(line++, 1, "Correcting lineEndX[1]");
						lineEndX[1] = al.x;
						
						savedMagnitude = lineMagnitudeSql;
						savedBl = bl;
						savealDepth = alDepth;
						saveblDepth = blDepth;
						savedAl3dPoint = al3dPoint;
						savedBl3dPoint = bl3dPoint;
//						dxnl = lineEndX[1]-bl.x;
//						savealDepth = alDepth;
//						saveblDepth = blDepth;
					}
				//}
			}
			
			
		}
		
		
	}
	if (tr > 2) {
		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);
	}
	if (tl > 2) {
		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1]);
	}
	// draw the final fill line
	
	if ( lineCount > 0 ) {
//		mvprintw(line++, 0, "final line: %d,%d", lineStartX[1],lineEndX[1]);
//		drawHorizonalLineWithDepthBuffer(lineStartX[1], lineEndX[1], lineY[1], fill, lineDepthStart[1], lineDepthEnd[1], depthBuffer);
		
//		mvprintw(line++, 1, "Final line from %d-%d @ row %d", lineStartX[1], lineEndX[1], lineY[1]);
		//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
	}
		
}


void RenderPipeline::setWithDepthBuffer( Coordinates2D pt, char c, double depth) {
	if (depthBuffer != NULL) {
		if ( pt.x < 0 || pt.y < 0 || depthBuffer->cols <= pt.x || depthBuffer->rows <= pt.y) {
			return;
		}
		if (1.0/depth > ((double*)depthBuffer->data)[pt.x + depthBuffer->cols*pt.y]) {
			ColorRGBA color;
			color.b = 255;
			color.g = 255;
			color.r = 255;
			color.a = c;
			setRenderBuffer(pt.x, pt.y, color);
			
			((double*)depthBuffer->data)[pt.x + depthBuffer->cols*pt.y] = 1.0/depth;
//			return;
//
//
//			if (depth > 0.4) {
//				set(pt, '`');
//			} else if (depth > 0.35) {
//				set(pt, '.');
//			} else if (depth > 0.3) {
//				set(pt, ',');
//			} else if (depth > 0.25) {
//				set(pt, '=');
//			} else if (depth > 0.2) {
//				set(pt, '+');
//			} else if (depth > 0.15) {
//				set(pt, '*');
//			} else if (depth > 0.10) {
//				set(pt, '#');
//			} else if (depth > 0.05) {
//				set(pt, '%');
//			} else {
//				set(pt, '@');
			//
			//			}
		}
	} else {
		ColorRGBA color;
		color.b = 255;
		color.g = 255;
		color.r = 255;
		color.a = c;
		setRenderBuffer(pt.x, pt.y, color);
//		set(pt, c);
	}
}

void RenderPipeline::setWithDepthBuffer( Coordinates4D pt, char c, double depth) {
	setWithDepthBuffer(onlyXY(pt), c, depth);
}

void RenderPipeline::drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData) {
	Coordinates2D pixel;
	int diff = end.pixel.x-start.pixel.x;
	if (abs(diff) < 2) {
		return;
	}
	if (start.pixel.y != end.pixel.y) {
//		mvprintw(20, 10, "start.pixel.y != end.pixel.y: %d,%d", start.pixel.y, end.pixel.y);
		return;
	}
	
	int increment = end.pixel.x > start.pixel.x ? 1 : -1;
	double factor;
	
	Coordinates4D point;
	Coordinates3D normal;
	double invDepth;
	
	pixel.y = start.pixel.y;
	for (int i = start.pixel.x+increment; i != end.pixel.x; i += increment) {
		pixel.x = i;
		
		factor = ((double)(i - start.pixel.x)/(double)diff);
		
//		point = vectorSubtract(point2, point1);
//		point.x *= factor;
//		point.y *= factor;
//		point.z *= factor;
//		point.w *= factor;
//		point = vectorAdd(point, point1);
//		point = interpolate(point1, point2, factor);
//		point.z = 1.0/(1.0/point1.z + factor*(1.0/point2.z - 1.0/point1.z));
		
		point = perspectiveInterpolate(start.location3D, end.location3D, 0, 0, 0, factor);
		
//		depth = 1.0/(1.0/depth1 + factor*(1.0/depth2 - 1.0/depth1));
		invDepth = (1.0/depth1 + factor*(1.0/depth2 - 1.0/depth1));
		
//		normal = vectorSubtract(normal2, normal1);
//		normal.x *= factor;
//		normal.y *= factor;
//		normal.z *= factor;
//		normal = vectorAdd(normal, normal1);
		
//		normal = interpolate(start.normal, end.normal, factor);
		normal = perspectiveInterpolate(start.normal, end.normal, start.location3D.z, end.location3D.z, point.z, factor);
		normal = normalizeVectorFast(normal);
		
		setWithShader(pixel, invDepth, point, normal, userData);
		
	}
	
}


void RenderPipeline::setWithShader( Coordinates2D& pixel, double invDepth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData) {
	if ( pixel.x < 0 || pixel.y < 0 || depthBuffer->cols <= pixel.x || depthBuffer->rows <= pixel.y) {
		return;
	}
	int depthIndex = pixel.x + depthBuffer->cols*pixel.y;
	if (invDepth > ((double*)depthBuffer->data)[depthIndex]) {
		((double*)depthBuffer->data)[depthIndex] = invDepth;
		
		
		ColorRGBA colorOutput;
		FragmentInfo fInfo;
		fInfo.pixel = pixel;
		fInfo.location3D = pt3D;
		fInfo.normal = normal;
		fInfo.data = userData;
		fInfo.colorOutput = &colorOutput;
		fragmentShader(fInfo);
		
		setRenderBuffer(pixel.x, pixel.y, colorOutput);
	}
}


void RenderPipeline::setFloatDotWithDepthBuffer( double x, double y, double depth) {
	if (depthBuffer != NULL) {
		if ( (int)x < 0 || (int)y < 0 || depthBuffer->cols <= (int)x || depthBuffer->rows <= (int)y) {
			return;
		}
		if (1.0/depth > ((double*)depthBuffer->data)[(int)x + depthBuffer->cols*(int)y]) {
			drawDotFloat(x, y);
			((double*)depthBuffer->data)[(int)x + depthBuffer->cols*(int)y] = 1.0/depth;
		}
	} else {
		drawDotFloat(x, y);
	}
}

void RenderPipeline::drawDotFloat(double x, double y) {
	char ch = '`';
	
	double remainder = fmod(y, 1);
	
	// ,.-*'` <- fonts differ in terminal
	
	// 28 pixel charaxter height
	// , -> 23
	// . -> 21  -> separation for above = (23-21)/2/28
	// - -> 16
	// * -> 11
	// ' -> 9
	// ` -> 8
	if (remainder > ((23.0+21.0)/(2.0*28.0))) {
		ch = ',';
	} else if (remainder > ((21.0+16.0)/(2.0*28.0))) {
		ch = '.';
	} else if (remainder > ((16.0+11.0)/(2.0*28.0))) {
		ch = '-';
	} else if(remainder >  ((11.0+9.0)/(2.0*28.0))) {
		ch = '*';
	} else if(remainder >  ((9.0+8.0)/(2.0*28.0))) {
		ch = '\'';
	}
	
	//mvaddch( y, x, ch);
	
		ColorRGBA color;
		color.b = 255;
		color.g = 255;
		color.r = 255;
		color.a = ch;
		setRenderBuffer(x, y, color);
}



void RenderPipeline::setRenderBuffer(int x, int y, ColorRGBA& color) {
	((ColorRGBA*)(fbo[0].data))[x + fbo[0].cols * y] = color;
}


void RenderPipeline::renderBufferToTerminal() {
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
#ifndef FB_SUPPORT
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
			
#else
			int offsetfb = (y * (fb_bytes_per_length) + x)*fb_bytes;
			uint16_t finalcolor = 0;
			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].r & 0xF8) << (11-3);
			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].g & 0xFC) << (5-2);
			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].b & 0xF8) >> (3);
//			fbdata[2 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].r;
//			fbdata[1 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].g;
//			fbdata[0 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].b;
			
			*(uint16_t*)&fbdata[0 + offsetfb] = finalcolor;
			
#endif
		}
	}
}

void RenderPipeline::depthBufferToTerminal() {
	Coordinates2D pixel;
	Coordinates3D color;
	int offset, index;
	
	for (int y = 0; y < depthBuffer->rows; y++) {
		offset = y * depthBuffer->cols;
		for (int x = 0; x < depthBuffer->cols; x++) {
			index = x + offset;
			
			pixel.x = x;
			pixel.y = y;
			
			uint8_t luminosity = (((double*)depthBuffer->data)[index]-1)*200000;
			
			color.x = (double)luminosity / 255.0;
			color.y = (double)luminosity / 255.0;
			color.z = (double)luminosity / 255.0;
#ifndef FB_SUPPORT
//			if (((ColorRGBA*)fbo[0].data)[index].a == 0) {
				setRGB(pixel, clipRGB(color));
//			} else {
//				Coordinates3D clippedRGB = clipRGB(color);
//				Coordinates3D hsl = rgbToHsv(clippedRGB);
//
//				int hueIndex = floor(hsl.x + 1);
//
//				if (hsl.y < 0.33) {
//					hueIndex = 7;
//				}
//
//				attron(COLOR_PAIR(hueIndex));
//				set(pixel, ((ColorRGBA*)fbo[0].data)[index].a);
//				attroff(COLOR_PAIR(hueIndex));
//			}
			
#else
			int offsetfb = (y * (fb_bytes_per_length) + x)*fb_bytes;
			uint16_t finalcolor = 0;
			finalcolor |= (luminosity & 0xF8) << (11-3);
			finalcolor |= (luminosity & 0xFC) << (5-2);
			finalcolor |= (luminosity & 0xF8) >> (3);
//			fbdata[2 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].r;
//			fbdata[1 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].g;
//			fbdata[0 + offsetfb] = ((ColorRGBA*)fbo[0].data)[index].b;
			
			*(uint16_t*)&fbdata[0 + offsetfb] = finalcolor;
			
#endif
		}
	}
}
