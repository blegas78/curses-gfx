#include "curses-gfx-3d.h"
#include "curses-gfx.h"

#include <cmath>
#include <unistd.h>
#include <stdint.h>

#include <ncurses.h>

void drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&)) {
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
	double depth;
	
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
		
		depth = 1.0/(1.0/depth1 + factor*(1.0/depth2 - 1.0/depth1));
		
//		normal = vectorSubtract(normal2, normal1);
//		normal.x *= factor;
//		normal.y *= factor;
//		normal.z *= factor;
//		normal = vectorAdd(normal, normal1);
		normal = interpolate(start.normal, end.normal, factor);
//		normal = normalizeVector(normal);
		normal = normalizeVectorFast(normal);
		
		setWithShader(pixel, depth, point, normal, userData, depthBuffer, fragmentShader);
		
	}
	
}

void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, DepthBuffer* depthBuffer, void (*fragmentShader)(const FragmentInfo&), int &line) {
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
//		set(onlyXY(poly.vertices[i]), c++);
//		setWithDepthBuffer(onlyXY(poly.vertices[i]), 'o', poly.vertices[i].z, depthBuffer);
		setFloatDotWithDepthBuffer(poly.vertices[i].x+0.5, poly.vertices[i].y+0.5, poly.vertices[i].z, depthBuffer);
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
			depthsr[2] = 1.0/(1.0/brDepth + alphar*(1.0/arDepth - 1.0/brDepth));
//			lineDepthStart = depthsr[2];
//
//			normalr = interpolate(brNormal, arNormal, alphar);
//			point3dr = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, depthsr[2], alphar);
			

		
		
		if (ar.x == br.x && ar.y == br.y && indexRight != maxIndex) {
			// new line time
			
			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
			
			
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
			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
			
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
					
//					dxnl = lineEndX[1]-bl.x;
//					dynl = (al.y-1)-bl.y;
//					alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/lineMagnitudeSql;
////					depthsl[0] = depthsl[1];
////					depthsl[1] = depthsl[2];
////		//			depthsl[2] = blDepth - alphal*(blDepth - alDepth);
////					depthsl[2] = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
////					lineDepthEnd[2] = depthsl[2];
//
//					lineDepthEnd = 1.0/(1.0/blDepth + alphal*(1.0/alDepth - 1.0/blDepth));
////
					
					
					dxnl = lineEndX[1]-savedBl.x;
					dynl = (al.y-1)-savedBl.y;
					alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl))/savedMagnitude;
					lineDepthEnd = 1.0/(1.0/saveblDepth + alphal*(1.0/savealDepth - 1.0/saveblDepth));
					
//					normall[0] = normall[1];
//					normall[1] = normall[2];
//					normall = interpolate(blNormal, alNormal, alphal);
//					point3dl = perspectiveInterpolate(bl3dPoint, al3dPoint, blDepth, alDepth, lineDepthEnd, alphal);
					
					fragEnd.pixel.x = lineEndX[1];
					fragEnd.pixel.y = al.y-1;
					fragEnd.normal = interpolate(blNormal, alNormal, alphal);
//					fragEnd.location3D = perspectiveInterpolate(bl3dPoint, al3dPoint, blDepth, alDepth, lineDepthEnd, alphal);
					
					fragEnd.location3D = perspectiveInterpolate(savedBl3dPoint, savedAl3dPoint, saveblDepth, savealDepth, lineDepthEnd, alphal);
					
					
//						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
					//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//					drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], al.y-1, lineDepthStart, lineDepthEnd, point3dr,  point3dl, normalr, normall, modelView, userData, depthBuffer, fragmentShader);
					static bool shouldPlot = true;
					if (shouldPlot) {
//						shouldPlot = false;
						drawHorizonalLineWithShader(fragStart, fragEnd, lineDepthStart, lineDepthEnd, userData, depthBuffer, fragmentShader);
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
		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
	}
	if (tl > 2) {
		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1], depthBuffer);
	}
	// draw the final fill line
	
	if ( lineCount > 0 ) {
//		mvprintw(line++, 0, "final line: %d,%d", lineStartX[1],lineEndX[1]);
//		drawHorizonalLineWithDepthBuffer(lineStartX[1], lineEndX[1], lineY[1], fill, lineDepthStart[1], lineDepthEnd[1], depthBuffer);
		
//		mvprintw(line++, 1, "Final line from %d-%d @ row %d", lineStartX[1], lineEndX[1], lineY[1]);
		//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
	}
		
}

