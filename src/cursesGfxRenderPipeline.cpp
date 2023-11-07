#include "curses-gfx-handler.h"

#include "curses-gfx-3d.h"

#include <ncurses.h>
RenderPipeline::RenderPipeline() {
	fbo = new FrameBuffer[1];
	depthBuffer = new FrameBuffer;
	
	asTexImage2d(fbo, FBT_RGBA, 1, 1);
	asTexImage2d(depthBuffer, FBT_DEPTH, 1, 1);
	depthClearColor = 0;
    backfaceCulling = true;
#ifdef FB_SUPPORT
	setupLinuxFb();
#endif
    
    mRasterizerThreadPool.Start(std::thread::hardware_concurrency());
}


void RenderPipeline::mvaddstring(int x, int y, const char* string) {

    
    if (fbo->type == FBT_RGBA) {
        for (int i=x; i < x+strlen(string); i++) {
//            ColorRGBA* location = fbo->at(i,y);
//            *location =
            fbo->set(i,y, ColorRGBA({255, 255, 255, string[i-x]}));
//            fbo->data[(i + y*fbo->cols)*4+0] = 255;
//            fbo->data[(i + y*fbo->cols)*4+1] = 255;
//            fbo->data[(i + y*fbo->cols)*4+2] = 255;
//            fbo->data[(i + y*fbo->cols)*4+3] = string[i-x];
//            ((double*)depthBuffer->data)[i + y*depthBuffer->cols] = 0;
            fbo->set(i, y, double(0));
        }
    }
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
	
//	fb_data_size = fb_width * fb_height * fb_bytes * fb_bytes_per_length;
	fb_data_size = fb_height * fb_bytes_per_length* fb_bytes;

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
//				fwrite( &((ColorRGBA*)fbo->data)[i], 1, 3, fp);
                fwrite( &fbo->at<ColorRGBA>(i, 0), 1, 3, fp);
			}
			break;
			
		case FBT_DEPTH:
			fprintf(fp, "P5\n%d %d\n255\n", fbo->cols, fbo->rows);
			for (int i = 0; i < fbo->cols*fbo->rows; i++) {
				uint8_t data;
//				if ( ((double*)fbo->data)[i] == 0 ) {
                if ( fbo->at<double>(i, 0) == 0 ) {
					data = 0;
				} else {
//					data = (((double*)fbo->data)[i]-1)*200000;
                    data = (fbo->at<double>(i, 0)-1)*200000;
				}
//				printf("%f == %d\n", ((double*)fbo->data)[i], data);
				fwrite( &data, 1, 1, fp);
			}
			
		default:
			break;
	}
	
	fclose(fp);
}

template <typename T> T perspectiveInterpolateBary(T& a, T& b, T& c, double& aInvDepth, double& bInvDepth, double& cInvDepth, double& correctDepth, double& alpha, double& beta, double& gamma) {
    return correctDepth*(a*(aInvDepth * alpha) + b*(bInvDepth * beta) + c*(cInvDepth * gamma));
}

template <typename T> T perspectiveIntBary(T& a, T& b, T& c, double& aInvDepthAlpha, double& bInvDepthBeta, double& cInvDepthGamma, double& correctDepth) {
    return correctDepth*(a*aInvDepthAlpha + b*bInvDepthBeta + c*cInvDepthGamma);
}
template <typename T, int count> void perspectiveIntBary2(void* result, const void* a, const void* b, const void* c, double& aInvDepthAlpha, double& bInvDepthBeta, double& cInvDepthGamma, double& correctDepth) {
//    T result;

    for(int i = count-1; i >= 0; i--)
        ((T*)result)[i] = correctDepth*( ((T*)a)[i] *aInvDepthAlpha + ((T*)b)[i]*bInvDepthBeta + ((T*)c)[i]*cInvDepthGamma);
//    return result;
}

RenderPipeline::~RenderPipeline() {
	// Fuck it, save the render and depth buffer to a file:
    mRasterizerThreadPool.Stop();
    
	saveFrameBufferToFile("render.ppm", &fbo[0]);
	saveFrameBufferToFile("depth.pgm", &depthBuffer[0]);
	
	delete [] fbo;
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
	clearColor.r = 0;   // Normal color Mode
	clearColor.g = 0;
	clearColor.b = 0;
	clearColor.a = 0;
//    clearColor.r = 255;   // Light Mode
//    clearColor.g = 255;
//    clearColor.b = 255;
//    clearColor.a = 0;
	fbo->clear(clearColor);
}

void RenderPipeline::setFragmentShader(void (*fragmentShader)(const FragmentInfo&)) {
//    RasterizerThreadPool::waitThreads(0);
    if(this->fragmentShader != fragmentShader)
        mRasterizerThreadPool.busyWait();
	this->fragmentShader = fragmentShader;
};


//void* RenderPipeline::renderThread(void* info) {
//	Renderable* This = (Renderable*)info;
//
//	//This->rasterizePolygonsShader();
//
//	This->pipeline->rasterizePolygonsShader(This->polygons,
//											This->count,
//											This->modelView,
//											This->projection,
//											This->viewPort,
//											This->userData,
//											This->line);
//
//
//	delete [] This->polygons;
//	delete This;
//	//return NULL;
//    pthread_exit(NULL);
//}
//
//void RenderPipeline::rasterizeThreaded(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line) {
//	Renderable* thisRenderable = new Renderable;
//	
//	thisRenderable->pipeline = this;
//	thisRenderable->polygons = new Polygon4D[count];
//	for (int i = 0; i < count; i++) {
//		thisRenderable->polygons[i].numVertices = polygons[i].numVertices;
//		
//		for (int v = 0; v < polygons[i].numVertices; v++) {
//			thisRenderable->polygons[i].normals[v] = polygons[i].normals[v];
//			thisRenderable->polygons[i].vertices[v] = polygons[i].vertices[v];
//			
//		}
//	}
////	memcpy(thisRenderable->polygons, polygons, sizeof(Polygon4D)*count);
//	
//	thisRenderable->count = count;
//	thisRenderable->modelView = modelView;
//	thisRenderable->projection = projection;
//	thisRenderable->viewPort = viewport;
//	thisRenderable->userData = userData;
//	thisRenderable->line = line;
//	
//	thisRenderable->fragmentShader = this->fragmentShader;	// hmmm
//	
//	pthread_t* thread = new pthread_t;
//	pthread_create(thread, NULL, renderThread, (void *)thisRenderable);
//	renderThreads.push(thread);
//	
//	waitForThreads();
//}
//
//void RenderPipeline::waitForThreads() {
//	while (!renderThreads.empty()) {
//		pthread_t* thread = renderThreads.front();
//		pthread_join(*thread, NULL);
//		renderThreads.pop();
//		delete thread;
//	}
//	
//}

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
			
//			rasterizeThreaded(&quadAsPolygon, 1, modelView, projection, viewport, userData, line);
			
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
		int valid = ::clipPolygon(polygonProjected, &polygonProjected, line);
		
		// now get back the original coordinates for separate handling
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
			polygonProjected.vertices[v] = matrixVectorMultiply(viewport, polygonProjected.vertices[v]);
		}
		
		// Draw:
////		mvprintw(line++, 0, "Calling drawPolygon() with %d vertices", polygon.numVertices);
////		refresh();
////		drawPolygon( polygon, depthBuffer, fill, line);
////		drawPolygonShader( polygonProjected, polygonRestored, userData, depthBuffer, fragmentShader, line);
//		drawPolygonShader( polygonProjected, polygonRestored, userData, line);
////		mvprintw(line++, 0, "Calling drawPolygon() Done");
////		refresh();
//
//
		drawPolygonShader( polygonProjected, polygonRestored, userData, line);
//		drawPolygonWithTriangles( polygonProjected, polygonRestored, userData);
		
	}
}

void RenderPipeline::drawPolygonWithTriangles( Polygon4D& poly, Polygon4D& restored, void* userData) {
//    this->userData = userData;
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
	fragments[0].normal = poly.normals[0];
	for (int i = 2; i < poly.numVertices; i++) {
		fragments[1].location3D = poly.vertices[i-1];
		fragments[2].location3D = poly.vertices[i];
		fragments[1].color = poly.colors[i-1];
		fragments[2].color = poly.colors[i];
		fragments[1].normal = poly.normals[i-1];
		fragments[2].normal = poly.normals[i];
		triangleFill(fragments);
	}
	
	
}

void RenderPipeline::triangleFill(FragmentInfo* fragments) {
 //https://web.archive.org/web/20050408192410/http://sw-shader.sourceforge.net/rasterizer.html
    // 28.4 fixed-point coordinates
    
    const int Y1 = round(16.0f * fragments[0].location3D.y);
        const int Y2 = round(16.0f * fragments[2].location3D.y);
        const int Y3 = round(16.0f * fragments[1].location3D.y);

    const int X1 = round(16.0f * fragments[0].location3D.x);
    const int X2 = round(16.0f * fragments[2].location3D.x);
    const int X3 = round(16.0f * fragments[1].location3D.x);
    
    
    double invDepth1 = 1.0/fragments[0].location3D.z;
    double invDepth2 = 1.0/fragments[2].location3D.z;
    double invDepth3 = 1.0/fragments[1].location3D.z;

        // Deltas
        const int DX12 = X1 - X2;
        const int DX23 = X2 - X3;
        const int DX31 = X3 - X1;

        const int DY12 = Y1 - Y2;
        const int DY23 = Y2 - Y3;
        const int DY31 = Y3 - Y1;

        // Fixed-point deltas
        const int FDX12 = DX12 << 4;
        const int FDX23 = DX23 << 4;
        const int FDX31 = DX31 << 4;

        const int FDY12 = DY12 << 4;
        const int FDY23 = DY23 << 4;
        const int FDY31 = DY31 << 4;

        // Bounding rectangle
        int minx = (std::min(X1, std::min(X2, X3)) + 0xF) >> 4;
        int maxx = (std::max(X1, std::max(X2, X3)) + 0xF) >> 4;
        int miny = (std::min(Y1, std::min(Y2, Y3)) + 0xF) >> 4;
        int maxy = (std::max(Y1, std::max(Y2, Y3)) + 0xF) >> 4;

        //barycentric parameters:
//    int A = X1*Y2 + X2*Y3 + X3*Y1 - X1*Y3 - X2*Y1 - X3*Y2;
    double A = 1.0/(double)(X1*(Y2-Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2));
    int aX2Y3mX3Y2 = X2*Y3 - X3*Y2;
    int bX3Y1mX1Y3 = X3*Y1 - X1*Y3;
//    int cX1Y2mX2Y1 = X1*Y2 - X2*Y1;
    double alpha, beta, gamma;
    
//        (char*&)colorBuffer += miny * stride;
//    ColorRGBA color = {255,255,0,'^'};
//    setFrameBufferRGBA(1, 1, fbo, color);
        // Half-edge constants
        int C1 = DY12 * X1 - DX12 * Y1;
        int C2 = DY23 * X2 - DX23 * Y2;
        int C3 = DY31 * X3 - DX31 * Y3;

        // Correct for fill convention
        if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
        if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
        if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

        int CY1 = C1 + DX12 * (miny << 4) - DY12 * (minx << 4);
        int CY2 = C2 + DX23 * (miny << 4) - DY23 * (minx << 4);
        int CY3 = C3 + DX31 * (miny << 4) - DY31 * (minx << 4);
    Coordinates2D pt;
    
        for(pt.y = miny; pt.y < maxy; pt.y++)
        {
            int CX1 = CY1;
            int CX2 = CY2;
            int CX3 = CY3;
       
            for(pt.x = minx; pt.x < maxx; pt.x++)
            {
                if(CX1 > 0 && CX2 > 0 && CX3 > 0)
                {
//                    colorBuffer[x] = 0x00FFFFFF;
                    int xp = pt.x << 4;
                    int yp = pt.y << 4;
                    alpha = (double)(xp*Y2 + aX2Y3mX3Y2 + X3*yp - xp*Y3 - X2*yp)*A;
                    beta = (double)(X1*yp + xp*Y3 + bX3Y1mX1Y3 - xp*Y1 - X3*yp)*A;
//                    gamma = (double)(cX1Y2mX2Y1 + X2*yp + xp*Y1 - X1*yp - xp*Y2)*A;
                    gamma = 1.0 - alpha - beta;
                    
                    ColorRGBA color;
//                    color.a = alpha * fragments[0].color.a + beta * fragments[2].color.a + gamma * fragments[1].color.a;
//                    color.r = alpha * fragments[0].color.r + beta * fragments[2].color.r + gamma * fragments[1].color.r;
//                    color.g = alpha * fragments[0].color.g + beta * fragments[2].color.g + gamma * fragments[1].color.g;
//                    color.b = alpha * fragments[0].color.b + beta * fragments[2].color.b + gamma * fragments[1].color.b;
                    
                    double correctInvDepth = invDepth1 * alpha + invDepth2 * beta + invDepth3 * gamma;
                    double correctDepth = 1.0/correctInvDepth;
//                    perspectiveInterpolateInv(<#T &a#>, <#T &b#>, <#T &c#>, <#double &aInvDepth#>, <#double &bInvDepth#>, <#double &cInvDepth#>, <#double &correctDepth#>, <#double &alpha#>, <#double &beta#>, <#double &gamma#>)
//                    (T& a, T& b, T& c, double& aInvDepth, double& bInvDepth, double& cInvDepth, double& correctDepth, double& alpha, double& beta, double& gamma)
//                    color.a = perspectiveInterpolateBary<uint8_t>(fragments[0].color.a, fragments[2].color.a, fragments[1].color.a, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//                    color.r = perspectiveInterpolateBary<uint8_t>(fragments[0].color.r, fragments[2].color.r, fragments[1].color.r, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//                    color.g = perspectiveInterpolateBary<uint8_t>(fragments[0].color.g, fragments[2].color.g, fragments[1].color.g, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//                    color.b = perspectiveInterpolateBary<uint8_t>(fragments[0].color.b, fragments[2].color.b, fragments[1].color.b, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
                    
                    alpha *= invDepth1;
                    beta *= invDepth2;
                    gamma *= invDepth3;
//                    color.a = perspectiveIntBary<uint8_t>(fragments[0].color.a, fragments[2].color.a, fragments[1].color.a, alpha, beta, gamma, correctDepth);
//                    color.r = perspectiveIntBary<uint8_t>(fragments[0].color.r, fragments[2].color.r, fragments[1].color.r, alpha, beta, gamma, correctDepth);
//                    color.g = perspectiveIntBary<uint8_t>(fragments[0].color.g, fragments[2].color.g, fragments[1].color.g, alpha, beta, gamma, correctDepth);
//                    color.b = perspectiveIntBary<uint8_t>(fragments[0].color.b, fragments[2].color.b, fragments[1].color.b, alpha, beta, gamma, correctDepth);
                    perspectiveIntBary2<uint8_t,4>(&color, &fragments[0].color, &fragments[2].color, &fragments[1].color, alpha, beta, gamma, correctDepth);
                    Coordinates3D normal;
                    perspectiveIntBary2<double,3>(&normal, &fragments[0].normal, &fragments[2].normal, &fragments[1].normal, alpha, beta, gamma, correctDepth);
                    Coordinates4D point;
                    perspectiveIntBary2<double,4>(&normal, &fragments[0].location3D, &fragments[2].location3D, &fragments[1].location3D, alpha, beta, gamma, correctDepth);
                    
                    
                    setFrameBufferRGBA(pt.x, pt.y, fbo, color);
//                    setFragmentShader(defaultFragment);
//                    setWithShader(pt, correctInvDepth, point, normal, this->userData);
                }

                CX1 -= FDY12;
                CX2 -= FDY23;
                CX3 -= FDY31;
            }

            CY1 += FDX12;
            CY2 += FDX23;
            CY3 += FDX31;

//            (char*&)colorBuffer += stride;
        }
}
//	ColorRGBA color1, color2;
//
//	FragmentInfo fragment;
//
//	std::vector<FragmentInfo> f(fragments, fragments+3);
//
////	for (int i = 0; i < 3; i++) {
////		fragment.location3D = vertices[i];
////		fragment.color = colors[i];
////		f.push_back(fragment);
////	}
//	std::sort(f.begin(), f.end(), compareCoordinatesFrag);
//
////	std::vector<Coordinates4D> v(vertices, vertices+3);
////	std::sort(v.begin(), v.end(), compareCoordinates);
//
//
////	drawHorizonalLineRGBA(v[0].x, v[0].x, v[0].y, color, fbo);
////	setFrameBufferRGBA(v[0].x, v[0].y, fbo, color);
////	setFrameBufferRGBA(v[1].x, v[1].y, fbo, color);
////	setFrameBufferRGBA(v[2].x, v[2].y, fbo, color);
//
//	double slope12 = (f[1].location3D.x - f[0].location3D.x)/(f[1].location3D.y-f[0].location3D.y);
//	if((int)f[1].location3D.y == (int)f[0].location3D.y) {
//		slope12 = 0;
//	}
//	double slope13 = (f[2].location3D.x - f[0].location3D.x)/(f[2].location3D.y-f[0].location3D.y);
//
//	double y = f[0].location3D.y;
//	double yp = 0;// = y - v[0].y; <- add v[0].y to get y
//	double x1, x2;
//
////	color.r = 0;
//
//	// y goes from v[0].y to v[2].y;
//	// yp goes from 0 to v[2].y - v[0].y
//	// cf2 goes from 0 to 1, or along yp/(v[2].y - v[0].y)
//	double colorFactor = 1.0/(f[2].location3D.y - f[0].location3D.y);
//
//	// cf goes from 0 to 1, or along yp/(v[1].y - v[0].y)
//	double colorFactor2 = 1.0/(f[1].location3D.y - f[0].location3D.y);
//
//
//	for( ; (yp+f[0].location3D.y) < f[1].location3D.y; yp += 1.0) {
//
////		x1 = slope12*(y-v[0].y) + v[0].x;
////		x2 = slope13*(y-v[0].y) + v[0].x;
//		x1 = slope12*(yp) + f[0].location3D.x;
//		x2 = slope13*(yp) + f[0].location3D.x;
////		drawHorizonalLineRGBA(x1, x2, y, color, fbo);
//
////		color.r = color[0].r + colorFactor*(color[2].r - color[0].r);
//		color2 = interpolate(f[0].color, f[2].color, colorFactor*yp);
//		color1 = interpolate(f[0].color, f[1].color, colorFactor2*yp);
//
//		drawHorizonalLineRGBA(x1, x2, yp+f[0].location3D.y, color1, color2);
//	}
//	yp =f[1].location3D.y - f[0].location3D.y;
//
//	double slope23 = (f[2].location3D.x - f[1].location3D.x)/(f[2].location3D.y-f[1].location3D.y);
//	if((int)f[2].location3D.y == (int)f[1].location3D.y) {
//		slope23 = 0;
//	}
//
//	// now yp goes from v[1].y-v[0].y to v[2].y - v[0].y
//	// cf2 goes from 0 to 1, or along (yp-(v[1].y-v[0].y))/((v[2].y - v[0].y) - (v[1].y-v[0].y))
////	colorFactor2 = 1.0/(f[2].location3D.y - (f[1].location3D.y-f[0].location3D.y));
//	colorFactor2 = 1.0/((f[2].location3D.y-f[0].location3D.y) - (f[1].location3D.y-f[0].location3D.y));
//
//	double end = f[2].location3D.y;
//	for ( ; (yp+f[0].location3D.y) < end; yp += 1.0) {
//
////		x1 = slope23*(yp+v[0].y-v[1].y) + v[1].x;
////		x2 = slope13*(yp+v[0].y-v[0].y) + v[0].x;
//		x1 = slope23*(yp+f[0].location3D.y-f[1].location3D.y) + f[1].location3D.x;
//		x2 = slope13*(yp) + f[0].location3D.x;
//
//
//		color2 = interpolate(f[0].color, f[2].color, colorFactor*yp);
//		color1 = interpolate(f[1].color, f[2].color, colorFactor2*(yp-(f[1].location3D.y-f[0].location3D.y)));
//		drawHorizonalLineRGBA(x1, x2, yp+f[0].location3D.y, color1, color2);
//	}
//
//}

void RenderPipeline::drawHorizonalLineRGBA(double x1, double x2, int y, ColorRGBA& color1, ColorRGBA& color2) {
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
		setFloatDotWithDepthBuffer(poly.vertices[i].x+0.5, poly.vertices[i].y+0.5, 1.0/poly.vertices[i].z);
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
//	double arDepth = poly.vertices[minIndex].z;
//	double brDepth = poly.vertices[mod(minIndex+1, poly.numVertices)].z;
	double arInvDepth = 1.0/poly.vertices[minIndex].z;
	double brInvDepth = 1.0/poly.vertices[mod(minIndex+1, poly.numVertices)].z;
	
	Coordinates2D al = onlyXY(poly.vertices[minIndex]);
	Coordinates2D bl = onlyXY(poly.vertices[mod(minIndex-1, poly.numVertices)]);
	Coordinates3D alNormal = restored.normals[minIndex];
	Coordinates4D al3dPoint = restored.vertices[minIndex];
	Coordinates3D blNormal = restored.normals[mod(minIndex-1, poly.numVertices)];
	Coordinates4D bl3dPoint = restored.vertices[mod(minIndex-1, poly.numVertices)];
//	double alDepth = poly.vertices[minIndex].z;
//	double blDepth = poly.vertices[mod(minIndex-1, poly.numVertices)].z;
	double alInvDepth = 1.0/poly.vertices[minIndex].z;
	double blInvDepth = 1.0/poly.vertices[mod(minIndex-1, poly.numVertices)].z;
	
	// The right line:
	int dxr = abs(br.x - ar.x);
	int sxr = ar.x < br.x ? 1 : -1;
	int dyr = -abs(br.y - ar.y);
	int syr = ar.y < br.y ? 1 : -1;
	int errr = dxr + dyr;
	int e2r;
	double lineInvMagnitudeSqr = 1.0/sqrt(dxr*dxr + dyr*dyr);
	
	
	int errsr[3];
	double errsNormalizedr[3];
	Coordinates2D ptsr[3];
	
	double alphar;
	double dxnr, dynr;
//	double depthsr[3];
	double invDepthsr[3];
	
	// The left line:
	int dxl = abs(bl.x - al.x);
	int sxl = al.x < bl.x ? 1 : -1;
	int dyl = -abs(bl.y - al.y);
	int syl = al.y < bl.y ? 1 : -1;
	int errl = dxl + dyl;
	int e2l;
	double lineInvMagnitudeSql = 1.0/sqrt(dxl*dxl + dyl*dyl);
	
	int errsl[3];
	double errsNormalizedl[3];
	Coordinates2D ptsl[3];
	
	double alphal;
	double dxnl, dynl;
//	double depthsl[3];
	double invDepthsl[3];
	
	
	
	int lineStartX[3] = {0,0,0};// = a.x - sx;
	int lineEndX[3] = {0,0,0};
//	int lineY[3] = {0,0,0};// = a.y;
	//double lineDepthStart;//[3] = {0,0,0};
//	double lineDepthEnd;
	double lineInvDepthStart;//[3] = {0,0,0};
	double lineInvDepthEnd;
	
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
	
	double savearInvDepth, savebrInvDepth, savedInvMagnitudeR = 1;
	Coordinates2D savedBr;
	Coordinates4D savedBr3dPoint, savedAr3dPoint;
	
	double savealInvDepth, saveblInvDepth, savedInvMagnitude = 1;
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
//				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);
				setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), invDepthsr[1]);
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
//			alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr))/lineMagnitudeSqr;
			alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr)) * lineInvMagnitudeSqr;
			
			invDepthsr[0] = invDepthsr[1];
			invDepthsr[1] = invDepthsr[2];
			invDepthsr[2] = brInvDepth + alphar*(arInvDepth - brInvDepth);
//			lineDepthStart = depthsr[2];
//
//			normalr = interpolate(brNormal, arNormal, alphar);
//			point3dr = perspectiveInterpolate(br3dPoint, ar3dPoint, brDepth, arDepth, depthsr[2], alphar);
			

		
		
		if (ar.x == br.x && ar.y == br.y && indexRight != maxIndex) {
			// new line time
			
//			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1], depthBuffer);
//			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);
			setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), invDepthsr[1]);
			
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
//			arDepth = poly.vertices[indexRight].z;
			arInvDepth = 1.0/poly.vertices[indexRight].z;
//			arInvDepth = 1.0/restored.vertices[indexRight].z;
			
//			lineStartX[2] = ar.x;
			
			int nextIndex = mod(indexRight + 1, poly.numVertices);
			br = onlyXY(poly.vertices[nextIndex]);
			brNormal = restored.normals[nextIndex];
			br3dPoint = restored.vertices[nextIndex];
//			brDepth = poly.vertices[nextIndex].z;
			brInvDepth = 1.0/poly.vertices[nextIndex].z;
//			brInvDepth = 1.0/restored.vertices[nextIndex].z;
			
			dxr = abs(br.x - ar.x);
			sxr = ar.x < br.x ? 1 : -1;
			dyr = -abs(br.y - ar.y);
			syr = ar.y < br.y ? 1 : -1;
			errr = dxr + dyr;
//			e2r;
			lineInvMagnitudeSqr = 1.0/sqrt(dxr*dxr + dyr*dyr);
			
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
				lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
					
				
				dxnr = lineStartX[1]-savedBr.x;
				dynr = (ar.y-1)-savedBr.y;
				alphar = ((double)sqrt(dxnr*dxnr + dynr*dynr)) * savedInvMagnitudeR;
				
				lineInvDepthStart = savebrInvDepth + alphar*(savearInvDepth - savebrInvDepth);
				
				fragStart.pixel.x = lineStartX[1];
				fragStart.pixel.y = ar.y-1;
				fragStart.normal = interpolate(brNormal, arNormal, alphar);
				
//				fragStart.location3D = perspectiveInterpolateInv(savedBr3dPoint, savedAr3dPoint, alphar);
				fragStart.location3D = perspectiveInterpolateInv(savedBr3dPoint, savedAr3dPoint, 1.0/savedBr3dPoint.z, 1.0/savedAr3dPoint.z, 0, alphar);
				
				savedInvMagnitudeR = lineInvMagnitudeSqr;
				savedBr = br;
//				savearDepth = arDepth;
//				savebrDepth = brDepth;
				savearInvDepth = arInvDepth;
				savebrInvDepth = brInvDepth;
				savedAr3dPoint = ar3dPoint;
				savedBr3dPoint = br3dPoint;

			} else {
				if (abs(al.x - ar.x) < abs(al.x - lineStartX[2])) {
//					mvaddstr(line++, 0, string);
					lineStartX[2] = ar.x;// + (a2.x < a.x ? -1 : 1);
//					sprintf(string, "lineStartX[1] corrected = %d", lineStartX[2]);
					
					savedInvMagnitudeR = lineInvMagnitudeSqr;
					savedBr = br;
//					savearDepth = arDepth;
//					savebrDepth = brDepth;
					savearInvDepth = arInvDepth;
					savebrInvDepth = brInvDepth;
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
				
				setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), invDepthsl[1]);
				
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
			alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl)) * lineInvMagnitudeSql;
			
			invDepthsl[0] = invDepthsl[1];
			invDepthsl[1] = invDepthsl[2];
			invDepthsl[2] = (blInvDepth + alphal*(alInvDepth - blInvDepth));
			
		
		if (al.x == bl.x && al.y == bl.y && indexLeft != maxIndex) {
			// new line time
//			mvprintw(line++, 0, "New Left line!");
			
			setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), invDepthsl[1]);
			
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
//			alDepth = poly.vertices[indexLeft].z;
			alInvDepth = 1.0/poly.vertices[indexLeft].z;
			
//			if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
////						mvprintw(line++, 1, "Correcting lineEndX[1]");
//			lineEndX[1] = lineEndX[2] < al.x ? lineEndX[2] : al.x;
//			}
			
			int nextIndex = mod(indexLeft - 1, poly.numVertices);
			bl = onlyXY(poly.vertices[nextIndex]);
			blNormal = restored.normals[nextIndex];
			bl3dPoint = restored.vertices[nextIndex];
//			blDepth = poly.vertices[nextIndex].z;
			blInvDepth = 1.0/poly.vertices[nextIndex].z;
			
			
			
			dxl = abs(bl.x - al.x);
			sxl = al.x < bl.x ? 1 : -1;
			dyl = -abs(bl.y - al.y);
			syl = al.y < bl.y ? 1 : -1;
			errl = dxl + dyl;
//			e2l;
			lineInvMagnitudeSql = 1.0/sqrt(dxl*dxl + dyl*dyl);
			
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
					alphal = ((double)sqrt(dxnl*dxnl + dynl*dynl)) * savedInvMagnitude;
					
					lineInvDepthEnd = saveblInvDepth + alphal*(savealInvDepth - saveblInvDepth);
					
					fragEnd.pixel.x = lineEndX[1];
					fragEnd.pixel.y = al.y-1;
					fragEnd.normal = interpolate(blNormal, alNormal, alphal);
					
//					fragEnd.location3D = perspectiveInterpolateInv(savedBl3dPoint, savedAl3dPoint, alphal);
//					fragEnd.location3D = perspectiveInterpolateInv(savedBl3dPoint, savedAl3dPoint, saveblInvDepth, savealInvDepth, 0, alphal);
					fragEnd.location3D = perspectiveInterpolateInv(savedBl3dPoint, savedAl3dPoint, 1.0/savedBl3dPoint.z, 1.0/savedAl3dPoint.z, 0, alphal);
					
//						drawHorizonalLineWithDepthBuffer(lineStartX[0], lineEndX[0], lineY[0], fill, lineDepthStart[0], lineDepthEnd[0], depthBuffer);
					//drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], lineY[1], lineDepthStart[1], lineDepthEnd[1], point3dr[1],  point3dl[1], normalr[1], normall[1], modelView, userData, depthBuffer, fragmentShader);
//					drawHorizonalLineWithShader(lineStartX[1], lineEndX[1], al.y-1, lineDepthStart, lineDepthEnd, point3dr,  point3dl, normalr, normall, modelView, userData, depthBuffer, fragmentShader);
//					static bool shouldPlot = true;
//					if (shouldPlot) {
////						shouldPlot = false;
////						drawHorizonalLineWithShader(fragStart, fragEnd, lineDepthStart, lineDepthEnd, userData);
						this->drawHorizonalLineWithShader(fragStart, fragEnd, lineInvDepthStart, lineInvDepthEnd, userData);
//					} else {
//						shouldPlot = true;
//					}
					
					

				}
				
				lineEndX[1] = lineEndX[2];
				
				savedBl = bl;
				savedInvMagnitude = lineInvMagnitudeSql;
//				savealDepth = alDepth;
//				saveblDepth = blDepth;
				savealInvDepth = alInvDepth;
				saveblInvDepth = blInvDepth;
				savedAl3dPoint = al3dPoint;
				savedBl3dPoint = bl3dPoint;
				
			} else {
				errsl[2] -= dxl;
//				skipRightLine = true;
				
				//if (!(ar.x == br.x && ar.y == br.y)) {
					if (abs(al.x - lineStartX[1]) < abs(lineEndX[1] - lineStartX[1])) {
//						mvprintw(line++, 1, "Correcting lineEndX[1]");
						lineEndX[1] = al.x;
						
						savedInvMagnitude = lineInvMagnitudeSql;
						savedBl = bl;
//						savealDepth = alDepth;
//						saveblDepth = blDepth;
						savealInvDepth = alInvDepth;
						saveblInvDepth = blInvDepth;
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
//		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), depthsr[1]);
		setWithDepthBuffer( ptsr[1], getp( ptsr, errsNormalizedr[1] + 0.5), invDepthsr[1]);
	}
	if (tl > 2) {
//		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), depthsl[1]);
		setWithDepthBuffer( ptsl[1], getp( ptsl, errsNormalizedl[1] + 0.5), invDepthsl[1]);
	}
	// draw the final fill line
		
}


void RenderPipeline::setWithDepthBuffer( Coordinates2D pt, char c, double invDepth) {
	if (depthBuffer != NULL) {
		if ( pt.x < 0 || pt.y < 0 || depthBuffer->cols <= pt.x || depthBuffer->rows <= pt.y) {
			return;
		}
		int depthIndex = pt.x + depthBuffer->cols*pt.y;
//		if (invDepth > ((double*)depthBuffer->data)[depthIndex]) {
        if (invDepth > depthBuffer->at<double>(depthIndex)) {
			ColorRGBA color;
			color.b = 255;
			color.g = 255;
			color.r = 255;
			color.a = c;
//			setRenderBuffer(pt.x, pt.y, color);
			setRenderBuffer(depthIndex, color);
			
//			((double*)depthBuffer->data)[depthIndex] = invDepth;
            depthBuffer->at<double>(depthIndex) = invDepth;
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

void RenderPipeline::setWithDepthBuffer( Coordinates4D pt, char c, double invDepth) {
	setWithDepthBuffer(onlyXY(pt), c, invDepth);
}

//double perspectiveInterpolateInv(double& a, double& b, double &c, double& aInvDepth, double& bInvDepth, double& cInvDepth, double& correctDepth, double& alpha, double& beta, double& gamma) {
//    return correctDepth*(a*aInvDepth * alpha + a*bInvDepth * beta + c*cInvDepth * gamma);
//}



//double perspectiveInterpolateInv(double& a, double& b, double &c, double& aInvDepth, double& bInvDepth, double& cInvDepth, double& correctDepth, double& alpha, double& beta, double& gamma) {
//    return correctDepth*(a*aInvDepth * alpha + a*bInvDepth * beta + c*cInvDepth * gamma);
//}

Coordinates4D perspectiveInterpolateBary(Coordinates4D& a, Coordinates4D& b, Coordinates4D& c, double aInvDepth, double bInvDepth, double cInvDepth, double correctDepth, double& alpha, double& beta, double& gamma) {
    Coordinates4D result;

//     aInvDepth = 1.0/a.z;
//     bInvDepth = 1.0/b.z;
    result.z = 1.0/(aInvDepth *alpha + bInvDepth*beta + cInvDepth*gamma);    // This is Zt, where Z1 is with "a" and Z2 is with "b"
    result.x = perspectiveInterpolateBary<double>(a.x, b.x, c.x, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    result.y = perspectiveInterpolateBary<double>(a.y, b.y, c.y, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    result.w = perspectiveInterpolateBary<double>(a.w, b.w, c.w, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    
    
    return result;
}

Coordinates3D perspectiveInterpolateBary(Coordinates3D& a, Coordinates3D& b, Coordinates3D& c, double aInvDepth, double bInvDepth, double cInvDepth, double correctDepth, double& alpha, double& beta, double& gamma) {
    Coordinates3D result;

    result.x = perspectiveInterpolateBary<double>(a.x, b.x, c.x, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    result.y = perspectiveInterpolateBary<double>(a.y, b.y, c.y, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    result.y = perspectiveInterpolateBary<double>(a.z, b.z, c.z, aInvDepth, bInvDepth, cInvDepth, result.z, alpha, beta, gamma);
    
    return result;
}

void RenderPipeline::drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double invDepth1, double invDepth2, void* userData) {
	Coordinates2D pixel;
	int diff = end.pixel.x-start.pixel.x;
	if (abs(diff) < 2) {
		return;
	}
	double invDiff = 1.0/(double)diff;
	
	if (start.pixel.y != end.pixel.y) {
//		mvprintw(20, 10, "start.pixel.y != end.pixel.y: %d,%d", start.pixel.y, end.pixel.y);
		return;
	}
	
	int increment = end.pixel.x > start.pixel.x ? 1 : -1;
	double factor;
	
	Coordinates4D point;
	Coordinates3D normal;
	double invDepth;
	
	double startLocation3DInvDepth = 1.0/start.location3D.z;
	double endLocation3DInvDepth = 1.0/end.location3D.z;
	
	pixel.y = start.pixel.y;
	for (int i = start.pixel.x+increment; i != end.pixel.x; i += increment) {
		pixel.x = i;
		
		factor = ((double)(i - start.pixel.x))*invDiff;
		
//		point = vectorSubtract(point2, point1);
//		point.x *= factor;
//		point.y *= factor;
//		point.z *= factor;
//		point.w *= factor;
//		point = vectorAdd(point, point1);
//		point = interpolate(point1, point2, factor);
//		point.z = 1.0/(1.0/point1.z + factor*(1.0/point2.z - 1.0/point1.z));
		
//		point = perspectiveInterpolate(start.location3D, end.location3D, factor);
		point = perspectiveInterpolateInv(start.location3D, end.location3D, startLocation3DInvDepth, endLocation3DInvDepth, 0, factor);
		
		invDepth = (invDepth1 + factor*(invDepth2 - invDepth1));
		
//		normal = vectorSubtract(normal2, normal1);
//		normal.x *= factor;
//		normal.y *= factor;
//		normal.z *= factor;
//		normal = vectorAdd(normal, normal1);
		
//		normal = interpolate(start.normal, end.normal, factor);
		normal = perspectiveInterpolateInv(start.normal, end.normal, startLocation3DInvDepth, endLocation3DInvDepth, point.z, factor);
		normal = normalizeVectorFast(normal);
		
		setWithShader(pixel, invDepth, point, normal, userData);
		
	}
	
}


void RenderPipeline::setWithShader( Coordinates2D& pixel, double invDepth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData) {
	if ( pixel.x < 0 || pixel.y < 0 || depthBuffer->cols <= pixel.x || depthBuffer->rows <= pixel.y) {
		return;
	}
	int depthIndex = pixel.x + depthBuffer->cols*pixel.y;
//	if (invDepth > ((double*)depthBuffer->data)[depthIndex]) {
//		((double*)depthBuffer->data)[depthIndex] = invDepth;
        if (invDepth > depthBuffer->at<double>(depthIndex)) {
            depthBuffer->at<double>(depthIndex) = invDepth;
		
//        ColorRGBA colorOutput = ((ColorRGBA*)(fbo[0].data))[depthIndex];
		FragmentInfo fInfo;
		fInfo.pixel = pixel;
		fInfo.location3D = pt3D;
		fInfo.normal = normal;
		fInfo.data = userData;
//		fInfo.colorOutput = &colorOutput;
//        fInfo.colorOutput = &((ColorRGBA*)(fbo[0].data))[depthIndex];
        fInfo.colorOutput = &fbo[0].at<ColorRGBA>(depthIndex);
		fragmentShader(fInfo);
		
//		setRenderBuffer(depthIndex, colorOutput);
	}
}

void RenderPipeline::setWithShader2( Coordinates2D& pixel, double invDepth, void* interpolatedData, void* userData) {
//    if ( pixel.x < 0 || pixel.y < 0 || depthBuffer->cols <= pixel.x || depthBuffer->rows <= pixel.y) {
//        return;
//    }
    int depthIndex = pixel.x + depthBuffer->cols*pixel.y;
//    depthBuffer->mutex.lock();
    
//    if (invDepth > ((double*)depthBuffer->data)[depthIndex]) {
//        ((double*)depthBuffer->data)[depthIndex] = invDepth;
        if (invDepth > depthBuffer->at<double>(depthIndex)) {
            depthBuffer->at<double>(depthIndex) = invDepth;
//        depthBuffer->mutex.unlock();
        
//        ColorRGBA colorOutput = ((ColorRGBA*)(fbo[0].data))[depthIndex];
        FragmentInfo fInfo;
        fInfo.pixel = pixel;
        fInfo.data = userData;
        fInfo.interpolated = interpolatedData;
//        fInfo.colorOutput = &colorOutput;
//        fInfo.colorOutput = &((ColorRGBA*)(fbo[0].data))[depthIndex];
        fInfo.colorOutput = &fbo[0].at<ColorRGBA>(depthIndex);
        fragmentShader(fInfo);
        
//        setRenderBuffer(depthIndex, colorOutput);
//        return;
    }
//    depthBuffer->mutex.unlock();
}


void RenderPipeline::setFloatDotWithDepthBuffer( double x, double y, double invDepth) {
	if (depthBuffer != NULL) {
		if ( (int)x < 0 || (int)y < 0 || depthBuffer->cols <= (int)x || depthBuffer->rows <= (int)y) {
			return;
		}
        int depthIndex = (int)x + depthBuffer->cols*(int)y;
        if (invDepth > depthBuffer->at<double>(depthIndex)) {
            depthBuffer->at<double>(depthIndex) = invDepth;
//		if (invDepth > ((double*)depthBuffer->data)[(int)x + depthBuffer->cols*(int)y]) {
//			drawDotFloat(x, y);
//			((double*)depthBuffer->data)[(int)x + depthBuffer->cols*(int)y] = invDepth;
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



void RenderPipeline::setRenderBuffer(const int& x, const int& y, ColorRGBA& color) {
//	((ColorRGBA*)(fbo[0].data))[x + fbo[0].cols * y] = color;
    fbo[0].at<ColorRGBA>(x, y) = color;
}
void RenderPipeline::setRenderBuffer(const int& index, ColorRGBA& color) {
//	((ColorRGBA*)(fbo[0].data))[index] = color;
    fbo[0].at<ColorRGBA>(index) = color;
}


void RenderPipeline::renderBufferToTerminal(int fboIndex) {
//	waitForThreads();
//    RasterizerThreadPool::waitThreads(0);
    mRasterizerThreadPool.busyWait();
	
	Coordinates2D pixel;
	Coordinates3D color;
	int offset, index;
	
	for (int y = 0; y < fbo[fboIndex].rows; y++) {
		offset = y * fbo[fboIndex].cols;
		for (int x = 0; x < fbo[fboIndex].cols; x++) {
			index = x + offset;
			
			pixel.x = x;
			pixel.y = y;
            
			color.x = (double)fbo[fboIndex].at<ColorRGBA>(index).r / 255.0;
			color.y = (double)fbo[fboIndex].at<ColorRGBA>(index).g / 255.0;
			color.z = (double)fbo[fboIndex].at<ColorRGBA>(index).b / 255.0;
#ifndef FB_SUPPORT
			if (fbo[fboIndex].at<ColorRGBA>(index).a == 0) {
				CursesGfxTerminal::setRGB(pixel, color);
			} else {
                double dummyLevel;
                CursesGfxTerminal::enableColor(color, dummyLevel);
//				set(pixel, ((ColorRGBA*)fbo[0].data)[index].a);
                set(pixel, fbo[fboIndex].at<ColorRGBA>(index).a);
                CursesGfxTerminal::disableColor();
			}
			
#else
			int offsetfb = (y * (fb_bytes_per_length) + x)*fb_bytes;
//			uint16_t finalcolor = 0;
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].r & 0xF8) << (11-3);
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].g & 0xFC) << (5-2);
//			finalcolor |= (((ColorRGBA*)fbo[0].data)[index].b & 0xF8) >> (3);
//			*(uint16_t*)&fbdata[0 + offsetfb] = finalcolor;
			
			fbdata[2 + offsetfb] = ((ColorRGBA*)fbo[fboIndex].data)[index].r;
   			fbdata[1 + offsetfb] = ((ColorRGBA*)fbo[fboIndex].data)[index].g;
   			fbdata[0 + offsetfb] = ((ColorRGBA*)fbo[fboIndex].data)[index].b;
			
#endif
		}
	}
}

void RenderPipeline::depthBufferToTerminal() {
//	waitForThreads();
//    RasterizerThreadPool::waitThreads(0);
    mRasterizerThreadPool.busyWait();
	
	Coordinates2D pixel;
	Coordinates3D color;
	int offset, index;
    
    double minDepth = std::numeric_limits<double>::max(), maxDepth = std::numeric_limits<double>::min();
    for (int y = 0; y < depthBuffer->rows; y++) {
        offset = y * depthBuffer->cols;
        for (int x = 0; x < depthBuffer->cols; x++) {
            index = x + offset;
//            if(((double*)depthBuffer->data)[index] < minDepth) {
//              minDepth = ((double*)depthBuffer->data)[index];
            if(depthBuffer->at<double>(index) < minDepth) {
                minDepth = depthBuffer->at<double>(index);
            }
            if(depthBuffer->at<double>(index) > maxDepth) {
//                maxDepth = ((double*)depthBuffer->data)[index];
                maxDepth = depthBuffer->at<double>(index);
            }
        }
    }
    //minDepth = 0;
    double scale = 255.0/(1-(1.0/100.0));
    scale = 2500;
    minDepth = 1;
	for (int y = 0; y < depthBuffer->rows; y++) {
		offset = y * depthBuffer->cols;
		for (int x = 0; x < depthBuffer->cols; x++) {
			index = x + offset;
			
			pixel.x = x;
			pixel.y = y;
			
//			double zFar = 100;
//			double zNear = 0.001;
//
//			double z = 1.0/((((double*)depthBuffer->data)[index] - (zFar+zNear)/(zFar-zNear)) * );
			
//			uint8_t luminosity = (((double*)depthBuffer->data)[index]-1)*150000;
//			int L = ((((double*)depthBuffer->data)[index])-1)*25500000.0;
//            int L = ((((double*)depthBuffer->data)[index])-minDepth)*scale;
            int L = (depthBuffer->at<double>(index)-minDepth)*scale;
//            int L = ((((double*)depthBuffer->data)[index]))*scale;
			//printf("d = %f\n",  1.0/((double*)depthBuffer->data)[index]);
			L = L > 255 ? 255 : L;
			L = L < 0 ? 0 : L;
			uint8_t luminosity = L;//(((double*)depthBuffer->data)[index]+1)*255.0/2.0;
			
			color.x = (double)luminosity / 255.0;
			color.y = (double)luminosity / 255.0;
			color.z = (double)luminosity / 255.0;
#ifndef FB_SUPPORT
//			if (((ColorRGBA*)fbo[0].data)[index].a == 0) {
            CursesGfxTerminal::setRGB(pixel, clipRGB(color));
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
//			uint16_t finalcolor = 0;
//			finalcolor |= (luminosity & 0xF8) << (11-3);
//			finalcolor |= (luminosity & 0xFC) << (5-2);
//			finalcolor |= (luminosity & 0xF8) >> (3);
//			*(uint16_t*)&fbdata[0 + offsetfb] = finalcolor;
			
			
			fbdata[2 + offsetfb] = luminosity;
			fbdata[1 + offsetfb] = luminosity;
			fbdata[0 + offsetfb] = luminosity;
#endif
		}
	}
}
