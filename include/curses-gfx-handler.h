#ifndef CURSES_GFX_HANDLER_H
#define CURSES_GFX_HANDLER_H

#include <pthread.h>
#include <queue>
#include <tuple>
#ifdef FB_SUPPORT
#include <linux/fb.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <sys/mman.h>
#endif

#include "curses-gfx.h"
#include "curses-gfx-3d.h"

class RenderPipeline;

typedef struct _Renderable {
	RenderPipeline* pipeline;
	Polygon4D* polygons;
	int count;
	Mat4D modelView;
	Mat4D projection;
	Mat4D viewPort;
	void* userData;
//	DepthBuffer* depthBuffer;
	void (*fragmentShader)(const FragmentInfo&);
	int line;
} Renderable;






class CursesGfxHandler {
private:
	pthread_t* threads;
	pthread_mutex_t mutex;
	pthread_mutex_t condMutex;
	pthread_cond_t cond;
	
	int totalThreads;
	bool busyThreads;
	
public:
	CursesGfxHandler();
	~CursesGfxHandler();
	
	void setTotalThreads(int numthreads);
	
};


class RenderPipeline {
public:
	FrameBuffer* fbo;
	FrameBuffer* depthBuffer;
	
	double depthClearColor;
	
	std::queue<pthread_t*> renderThreads;
	
//	DepthBuffer tempD;
    void* userData;
    
	void (*fragmentShader)(const FragmentInfo&);
	
	void drawPolygonWithTriangles( Polygon4D& poly, Polygon4D& restored, void* userData);
	void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, int &line);
	
	void triangleFill(FragmentInfo* fragments);
//    template <class T> void triangleFill(T* vertexInfo);
    template <class T> void triangleFill(T* fragment1, T* fragment2, T* fragment3);
    template <class T> void trianglesFill(T* vertexInfo, int edgeIndices[][3], int count);
	void drawHorizonalLineRGBA(double x1, double x2, int y, ColorRGBA& color1, ColorRGBA& color2);
	
	void setWithDepthBuffer( Coordinates2D pt, char c, double invDepth);
	void setWithDepthBuffer( Coordinates4D pt, char c, double invDepth);
	
	void drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData);
	
	void setFloatDotWithDepthBuffer( double x, double y, double depth);
	void drawDotFloat(double x, double y);
	
	void setWithShader( Coordinates2D& pixel, double invDepth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData);
    void setWithShader2( Coordinates2D& pixel, double invDepth, void* userData, void* interpolatedData);
    
	static void* renderThread(void* info);
	void waitForThreads();
	
#ifdef FB_SUPPORT
	int fbfd;
	struct fb_var_screeninfo vinfo;
	struct fb_fix_screeninfo finfo;
	int fb_width;
	int fb_height;
	int fb_bpp;
	int fb_bytes;
	int fb_bytes_per_length;
	
	int fb_data_size;

	char *fbdata;
	
	void setupLinuxFb();
#endif
public:
	RenderPipeline();
	~RenderPipeline();
	
	void resize(int width, int height);
	void reset(); // resets depth buffer
	
	void setFragmentShader(void (*fragmentShader)(const FragmentInfo&));
	
	void rasterizeQuadsShader(Coordinates4D* vertices, int quadIndices[][4], int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	void rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	void rasterizeThreaded(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	
	void setRenderBuffer(const int& x, const int& y, ColorRGBA& color);
	void setRenderBuffer(const int& index, ColorRGBA& color);
	
	void renderBufferToTerminal();
	void depthBufferToTerminal();
};


template<class T> void baryInterpolate(T& output, const T& input, const T& input2, const T& input3, const double& alpha, const double& beta, const double& gamma, std::integral_constant<int, 0>)
{
}
    
template<class T, int index = std::tuple_size_v<T>>
void baryInterpolate(T& output, const T& input, const T& input2, const T& input3, const double& alpha, const double& beta, const double& gamma, std::integral_constant<int, index> = std::integral_constant<int, std::tuple_size_v<T>>())
{
    *std::get<index-1>(output) = *std::get<index-1>(input)*alpha + *std::get<index-1>(input2)*beta + *std::get<index-1>(input3)*gamma;
    baryInterpolate(output, input, input2, input3, alpha, beta, gamma, std::integral_constant<int, index - 1>());
}


template <class T> void RenderPipeline::trianglesFill(T* vertexInfo, int edgeIndices[][3], int count) {
    
    for (int i = 0; i < count; i++) {
        this->triangleFill(&vertexInfo[edgeIndices[i][0]], &vertexInfo[edgeIndices[i][1]], &vertexInfo[edgeIndices[i][2]]);
    }
}

//template <class T> void RenderPipeline::triangleFill(T* fragments) {
template <class T> void RenderPipeline::triangleFill(T* fragment1, T* fragment2, T* fragment3) {
    //https://web.archive.org/web/20050408192410/http://sw-shader.sourceforge.net/rasterizer.html
       // 28.4 fixed-point coordinates
       
//       const int Y1 = round(16.0f * fragments[0].vertex.y);
//           const int Y2 = round(16.0f * fragments[1].vertex.y);
//           const int Y3 = round(16.0f * fragments[2].vertex.y);
//
//       const int X1 = round(16.0f * fragments[0].vertex.x);
//       const int X2 = round(16.0f * fragments[1].vertex.x);
//       const int X3 = round(16.0f * fragments[2].vertex.x);
//
//
//       double invDepth1 = 1.0/fragments[0].vertex.z;
//       double invDepth2 = 1.0/fragments[1].vertex.z;
//       double invDepth3 = 1.0/fragments[2].vertex.z;
    
    
    const int Y1 = round(16.0f * fragment1->vertex.y);
        const int Y2 = round(16.0f * fragment2->vertex.y);
        const int Y3 = round(16.0f * fragment3->vertex.y);

    const int X1 = round(16.0f * fragment1->vertex.x);
    const int X2 = round(16.0f * fragment2->vertex.x);
    const int X3 = round(16.0f * fragment3->vertex.x);
    
    
    double invDepth1 = 1.0/fragment1->vertex.z;
    double invDepth2 = 1.0/fragment2->vertex.z;
    double invDepth3 = 1.0/fragment3->vertex.z;
    double invW1 = 1.0/fragment1->vertex.w;
    double invW2 = 1.0/fragment2->vertex.w;
    double invW3 = 1.0/fragment3->vertex.w;

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
    T fragment;
    auto f1 = regist(*fragment1);
    auto f2 = regist(*fragment2);
    auto f3 = regist(*fragment3);
    auto output = regist(fragment);
       
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
                       
//                       ColorRGBA color;
   //                    color.a = alpha * fragments[0].color.a + beta * fragments[2].color.a + gamma * fragments[1].color.a;
   //                    color.r = alpha * fragments[0].color.r + beta * fragments[2].color.r + gamma * fragments[1].color.r;
   //                    color.g = alpha * fragments[0].color.g + beta * fragments[2].color.g + gamma * fragments[1].color.g;
   //                    color.b = alpha * fragments[0].color.b + beta * fragments[2].color.b + gamma * fragments[1].color.b;
                       
                       double correctInvDepth = invDepth1 * alpha + invDepth2 * beta + invDepth3 * gamma;
//                       double correctDepth = 1.0/correctInvDepth;
////                       double correctDepth = 1.0/invDepth1 * alpha + 1.0/invDepth2 * beta + 1.0/invDepth3 * gamma;
//   //                    perspectiveInterpolateInv(<#T &a#>, <#T &b#>, <#T &c#>, <#double &aInvDepth#>, <#double &bInvDepth#>, <#double &cInvDepth#>, <#double &correctDepth#>, <#double &alpha#>, <#double &beta#>, <#double &gamma#>)
//   //                    (T& a, T& b, T& c, double& aInvDepth, double& bInvDepth, double& cInvDepth, double& correctDepth, double& alpha, double& beta, double& gamma)
//   //                    color.a = perspectiveInterpolateBary<uint8_t>(fragments[0].color.a, fragments[2].color.a, fragments[1].color.a, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//   //                    color.r = perspectiveInterpolateBary<uint8_t>(fragments[0].color.r, fragments[2].color.r, fragments[1].color.r, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//   //                    color.g = perspectiveInterpolateBary<uint8_t>(fragments[0].color.g, fragments[2].color.g, fragments[1].color.g, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//   //                    color.b = perspectiveInterpolateBary<uint8_t>(fragments[0].color.b, fragments[2].color.b, fragments[1].color.b, invDepth1, invDepth2, invDepth3, correctDepth, alpha, beta, gamma);
//
//                       alpha *= invDepth1*correctDepth;
//                       beta *= invDepth2*correctDepth;
//                       gamma *= invDepth3*correctDepth;
                       
                       double pBaryDivisor = 1.0/(alpha*invW1 + beta*invW2 + gamma*invW3);
                       alpha *= invW1 * pBaryDivisor;
                       beta  *= invW2 * pBaryDivisor;
                       gamma *= invW3 * pBaryDivisor;
                       
                       
   //                    color.a = perspectiveIntBary<uint8_t>(fragments[0].color.a, fragments[2].color.a, fragments[1].color.a, alpha, beta, gamma, correctDepth);
   //                    color.r = perspectiveIntBary<uint8_t>(fragments[0].color.r, fragments[2].color.r, fragments[1].color.r, alpha, beta, gamma, correctDepth);
   //                    color.g = perspectiveIntBary<uint8_t>(fragments[0].color.g, fragments[2].color.g, fragments[1].color.g, alpha, beta, gamma, correctDepth);
   //                    color.b = perspectiveIntBary<uint8_t>(fragments[0].color.b, fragments[2].color.b, fragments[1].color.b, alpha, beta, gamma, correctDepth);
//                       perspectiveIntBary2<uint8_t,4>(&color, &fragments[0].color, &fragments[2].color, &fragments[1].color, alpha, beta, gamma, correctDepth);
//                       Coordinates3D normal;
//                       perspectiveIntBary2<double,3>(&normal, &fragments[0].normal, &fragments[2].normal, &fragments[1].normal, alpha, beta, gamma, correctDepth);
//                       Coordinates4D point;
//                       perspectiveIntBary2<double,4>(&normal, &fragments[0].location3D, &fragments[2].location3D, &fragments[1].location3D, alpha, beta, gamma, correctDepth);
                       
                       baryInterpolate(output, f1, f2, f3, alpha, beta, gamma);
//                       color.r = 255;
//                       color.g = 0;
//                       color.b = 255;
//                       color.a = 'V';
                       
//                       color = fragment.color;
//                       setFrameBufferRGBA(pt.x, pt.y, fbo, color);
   //                    setFragmentShader(defaultFragment);
   //                    setWithShader(pt, correctInvDepth, point, normal, this->userData);
                       
                       
                       setWithShader2( pt, correctInvDepth, NULL, (void*) &fragment);
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


#endif
