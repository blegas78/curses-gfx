#ifndef CURSES_GFX_HANDLER_H
#define CURSES_GFX_HANDLER_H

#include <pthread.h>
#include <queue>
#include <future>
#include <list>
#include <tuple>
#include <chrono>
#include <thread>

#include <unistd.h>

//#ifdef FB_SUPPORT
//#include <linux/fb.h>
//#include <sys/ioctl.h>
//#include <fcntl.h>
//#include <sys/mman.h>
//#endif

#include <chrono>

#include "curses-gfx.h"
#include "curses-gfx-3d.h"

class RenderPipeline;
template <class T> class BlockRenderer;

//typedef struct _Renderable {
//	RenderPipeline* pipeline;
//	Polygon4D* polygons;
//	int count;
//	Mat4D modelView;
//	Mat4D projection;
//	Mat4D viewPort;
//	void* userData;
////	DepthBuffer* depthBuffer;
//	void (*fragmentShader)(const FragmentInfo&);
//	int line;
//} Renderable;


typedef struct _RenderStats {
    double timeVertexShading;
    double timeClipping;
    double timeDrawing;
} RenderStats;

class RasterizerThreadPool {
private:
    void ThreadLoop();

    int activeThreadCount;
    bool should_terminate;           // Tells threads to stop looking for jobs
    std::mutex queue_mutex;                  // Prevents data races to the job queue
    std::condition_variable mutex_condition; // Allows threads to wait on new jobs or termination
    std::mutex busy_mutex;                  // Prevents data races to the job queue
    std::condition_variable mutex_cv_busy; // Allows threads to wait on new jobs or termination
    std::vector<std::thread> threads;
//    std::queue<std::function<void()>> jobs;
    std::queue<std::pair<std::function<void(void*)>,void*>> jobs;
    
public:

    // Following this: https://stackoverflow.com/questions/15752659/thread-pooling-in-c11
    
    void Start(uint32_t numThreads);
//    void QueueJob(const std::function<void()>& job);
    void QueueJob(const std::function<void(void*)>& job, void* data);
    void Stop();
    bool busy();
    void busyWait();
};


class RenderPipeline {
public:
	FrameBuffer* fbo;
	FrameBuffer* depthBuffer;
	
    bool backfaceCulling;
	double depthClearColor;
	
    RasterizerThreadPool mRasterizerThreadPool;
	std::queue<pthread_t*> renderThreads;
	
//	DepthBuffer tempD;
//    void* userData;
    
    Mat4D viewport;
    
//    void (*vertexShader)(...);
	void (*fragmentShader)(const FragmentInfo&);
	
	void drawPolygonWithTriangles( Polygon4D& poly, Polygon4D& restored, void* userData);
	void drawPolygonShader( Polygon4D& poly, Polygon4D& restored, void* userData, int &line);
	
	void triangleFill(FragmentInfo* fragments);
//    template <class T> void triangleFill(T* vertexInfo);
    template <class T> void triangleFill(T* fragment1, T* fragment2, T* fragment3, void* userData);
//    template <class T> void trianglesFill(T* vertexInfo, int edgeIndices[][3], int count);
	void drawHorizonalLineRGBA(double x1, double x2, int y, ColorRGBA& color1, ColorRGBA& color2);
	
	void setWithDepthBuffer( Coordinates2D pt, char c, double invDepth);
	void setWithDepthBuffer( Coordinates4D pt, char c, double invDepth);
	
	void drawHorizonalLineWithShader( FragmentInfo& start, FragmentInfo& end, double depth1, double depth2, void* userData);
	
	void setFloatDotWithDepthBuffer( double x, double y, double depth);
	void drawDotFloat(double x, double y);
	
	void setWithShader( Coordinates2D& pixel, double invDepth, Coordinates4D& pt3D, Coordinates3D& normal, void* userData);
    void setWithShader2( Coordinates2D& pixel, double invDepth, void* interpolatedData, void* userData);
    
//	static void *renderThread(void* info);
//	void waitForThreads();
	
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
    template <class T, class U> RenderStats rasterizeShader(T* vertexInfo, U* uniformInfo, int triangleLayout[][3], int numTriangles, void* userData, void (*vertexShader)(U*, T&, const T&));
    template <class T, class U> RenderStats rasterizeShader(T* vertexInfo, U* uniformInfo, int numTriangles, void* userData, void (*vertexShader)(U*, T&, const T&));
	void rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
//	void rasterizeThreaded(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	
	void setRenderBuffer(const int& x, const int& y, ColorRGBA& color);
	void setRenderBuffer(const int& index, ColorRGBA& color);
	
	void renderBufferToTerminal();
	void depthBufferToTerminal();
    
    template <class T> void clipPolygon(T* input, int inputCount, T* output, int& outputCountResult);
    template <class T> void clipW(T* input, int inputCount, T* output, int& outputCountResult);
    template <class T> void clipPlane(T* input, int inputCount, T* output, int& outputCountResult, int axis, int plane);
    
    
    void mvaddstring(int x, int y, const char* string);
};




template <class T, class U> RenderStats RenderPipeline::rasterizeShader(T* vertexInfo, U* uniformInfo, int triangleLayout[][3], int numTriangles, void* userData, void (*vertexShader)(U*, T&, const T&)) {
    
    auto now = std::chrono::high_resolution_clock::now();
    auto before = now;
    std::chrono::duration<double, std::milli> float_ms;
    
    T scratch[3];
//    int scratchLayout[][3] = {
//        {0, 1, 2}
//    };
    T scratchClipped[20];
    int clippedVertexCount;
    
    RenderStats mRenderStats = {0,0,0};
    
//    this->userData = userData;
    BlockRenderer<T> br(this);
    br.rtp = &mRasterizerThreadPool;
    for (int t = 0; t < numTriangles; t++) {
        
        for(int i = 0; i < 3; i++) {
            vertexShader(uniformInfo, scratch[i], vertexInfo[triangleLayout[t][i]]);
        }
        
        
        vertexShader(uniformInfo, scratch[0],  vertexInfo[triangleLayout[t][0]]);
        if(backfaceCulling) {
            vertexShader(uniformInfo, scratch[1],  vertexInfo[triangleLayout[t][1]]);
            vertexShader(uniformInfo, scratch[2],  vertexInfo[triangleLayout[t][2]]);
        } else {
            vertexShader(uniformInfo, scratch[2],  vertexInfo[triangleLayout[t][1]]);   // indices are flipped
            vertexShader(uniformInfo, scratch[1],  vertexInfo[triangleLayout[t][2]]);
        }
        
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeVertexShading += float_ms.count()/1000.0;
        
        // Then.. clip geometry
        clipPolygon(scratch, 3, scratchClipped, clippedVertexCount);
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeClipping += float_ms.count()/1000.0;
            
        for(int i = 0; i < clippedVertexCount; i++) {
            // Then.. apply a viewport and provide to triangle rasterizer
            scratchClipped[i].vertex = matrixVectorMultiply(viewport, scratchClipped[i].vertex);
        }
        
        
//            triangleFill(&scratch[0], &scratch[1], &scratch[2]);
        for(int i = 2; i < clippedVertexCount; i++) {
                br.triangleFill(&scratchClipped[0], &scratchClipped[i-1], &scratchClipped[i], userData);
        }
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeDrawing += float_ms.count()/1000.0;
//        fbo->data[20*4 + 3] = '0'+clippedVertexCount;

    }
//    mRasterizerThreadPool.busyWait();
    return mRenderStats;
}

template <class T, class U> RenderStats RenderPipeline::rasterizeShader(T* vertexInfo, U* uniformInfo, int numTriangles, void* userData, void (*vertexShader)(U*, T&, const T&)) {
    
    auto now = std::chrono::high_resolution_clock::now();
    auto before = now;
    std::chrono::duration<double, std::milli> float_ms;
    
    T scratch[3];
//    int scratchLayout[][3] = {
//        {0, 1, 2}
//    };
    T scratchClipped[20];
    int clippedVertexCount;
    
    RenderStats mRenderStats = {0,0,0};
    
//    this->userData = userData;
    BlockRenderer<T> br(this);
    br.rtp = &mRasterizerThreadPool;
    for (int t = 0; t < numTriangles; t++) {
        
        vertexShader(uniformInfo, scratch[0],  *vertexInfo++);
        if(backfaceCulling) {
            vertexShader(uniformInfo, scratch[1],  *vertexInfo++);
            vertexShader(uniformInfo, scratch[2],  *vertexInfo++);
        } else {
            vertexShader(uniformInfo, scratch[2],  *vertexInfo++);
            vertexShader(uniformInfo, scratch[1],  *vertexInfo++);
        }
        
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeVertexShading += float_ms.count()/1000.0;
        
        // Then.. clip geometry?
        clipPolygon(scratch, 3, scratchClipped, clippedVertexCount);
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeClipping += float_ms.count()/1000.0;
            
        for(int i = 0; i < clippedVertexCount; i++) {
            // Then.. apply a viewport and provide to triangle rasterizer
//            scratch[i].vertex = matrixVectorMultiply(viewport, scratch[i].vertex);
            scratchClipped[i].vertex = matrixVectorMultiply(viewport, scratchClipped[i].vertex);
        }
//        now = std::chrono::high_resolution_clock::now();
//        float_ms = (now - before);
//        before = now;
//        mRenderStats.timeVertexShading += float_ms.count()/1000.0;
        
        
//            triangleFill(&scratch[0], &scratch[1], &scratch[2]);
        for(int i = 2; i < clippedVertexCount; i++) {
                br.triangleFill(&scratchClipped[0], &scratchClipped[i-1], &scratchClipped[i], userData);
        }
        now = std::chrono::high_resolution_clock::now();
        float_ms = (now - before);
        before = now;
        mRenderStats.timeDrawing += float_ms.count()/1000.0;
//        fbo->data[20*4 + 3] = '0'+clippedVertexCount;

    }
//    RasterizerThreadPool::busyWait();   // finish rendering
//    mRasterizerThreadPool.busyWait();
    return mRenderStats;
}

template<class T> void baryInterpolate2d(T& output, const T& input, const T& input2, const double& alpha, const double& beta, std::integral_constant<int, 0>)
{
}
    
template<class T, int index = std::tuple_size_v<T>>
void baryInterpolate2d(T& output, const T& input, const T& input2, const double& alpha, const double& beta, std::integral_constant<int, index> = std::integral_constant<int, std::tuple_size_v<T>>())
{
    *std::get<index-1>(output) = *std::get<index-1>(input)*alpha + *std::get<index-1>(input2)*beta;
    baryInterpolate2d(output, input, input2, alpha, beta, std::integral_constant<int, index - 1>());
}

template<class T> void baryInterpolate(T& output, const T& input, const T& input2, const T& input3, const double& alpha, const double& beta, const double& gamma, std::integral_constant<int, 0>)
{
}
    
template<class T, int index = std::tuple_size_v<T>>
void baryInterpolate(T& output, const T& input, const T& input2, const T& input3, const double& alpha, const double& beta, const double& gamma, std::integral_constant<int, index> = std::integral_constant<int, std::tuple_size_v<T>>())
{
    *std::get<index-1>(output) = *std::get<index-1>(input)*alpha + *std::get<index-1>(input2)*beta + *std::get<index-1>(input3)*gamma;
    baryInterpolate(output, input, input2, input3, alpha, beta, gamma, std::integral_constant<int, index - 1>());
}


//template <class T> void RenderPipeline::trianglesFill(T* vertexInfo, int edgeIndices[][3], int count) {
//
//    for (int i = 0; i < count; i++) {
//        this->triangleFill(&vertexInfo[edgeIndices[i][0]], &vertexInfo[edgeIndices[i][1]], &vertexInfo[edgeIndices[i][2]]);
//    }
//}




//template <class T> std::list<std::future<void>> BlockRenderer<T>::threadStatus;
//template <class T> int BlockRenderer<T>::numberRenderThreads;

//std::list<std::future<void>> TriangleRasterizer::threadStatus;
//int TriangleRasterizer::numberRenderThreads = 1;

template <class T> class BlockRenderer {
private:
//    static int numberRenderThreads;
    struct RenderInfo {
//        BlockRenderer<T>* This;
        Coordinates2D pt;
        T fragment1, fragment2, fragment3;
        int CY1;
        int CY2;
        int CY3;
        // Uniform for full triangle:
        int q;
        decltype(regist(std::declval<T&>())) f1;
        decltype(regist(std::declval<T&>())) f2;
        decltype(regist(std::declval<T&>())) f3;
        int X1, X2, X3;
        int Y1, Y2, Y3;
        int aX2Y3mX3Y2, bX3Y1mX1Y3;
        double A;
        double invDepth1, invDepth2, invDepth3;
        double invW1, invW2, invW3;
        int FDX12;
        int FDX23;
        int FDX31;
        
        int FDY12;
        int FDY23;
        int FDY31;
        RenderPipeline* p;
        void* userData;
    };
    
    
//    static void renderThread(RenderInfo rip) {
      static void renderThread(void* data) {
        RenderInfo* rip = (RenderInfo*)data;  // TODO make this not another copy
//        RenderInfo* ri = ;
        Coordinates2D pt = rip->pt;
        //BlockRenderer<T>* This = (BlockRenderer<T>*) rip->This;
        int xp;
        int yp;
        double alpha, beta, gamma;
        double correctInvDepth;
        double pBaryDivisor;
        Coordinates2D ipt;
        
        T fragment;
        rip->f1 = regist(rip->fragment1);
        rip->f2 = regist(rip->fragment2);
        rip->f3 = regist(rip->fragment3);
        auto output = regist(fragment);
        for(ipt.y = pt.y; ipt.y < pt.y + rip->q; ipt.y++)
        {
            yp = ipt.y << 4;
            for(ipt.x = pt.x; ipt.x < pt.x + rip->q; ipt.x++)
            {
                xp = ipt.x << 4;
                alpha = (double)(xp*rip->Y2 + rip->aX2Y3mX3Y2 + rip->X3*yp - xp*rip->Y3 - rip->X2*yp)*rip->A;
                beta = (double)(rip->X1*yp + xp*rip->Y3 + rip->bX3Y1mX1Y3 - xp*rip->Y1 - rip->X3*yp)*rip->A;
                gamma = 1.0 - alpha - beta;
                
                correctInvDepth = rip->invDepth1 * alpha + rip->invDepth2 * beta + rip->invDepth3 * gamma;
                pBaryDivisor = 1.0/(alpha*rip->invW1 + beta*rip->invW2 + gamma*rip->invW3);
                alpha *= rip->invW1 * pBaryDivisor;
                beta  *= rip->invW2 * pBaryDivisor;
                gamma *= rip->invW3 * pBaryDivisor;
                
                baryInterpolate(output, rip->f1, rip->f2, rip->f3, alpha, beta, gamma);
                rip->p->setWithShader2( ipt, correctInvDepth, (void*)&fragment, rip->userData);
            }
        }
//        pthread_exit(NULL);
          
          delete rip;
    }
    
//    static void renderThread2(RenderInfo rip) {
        static void renderThread2(void* data) {
        RenderInfo* rip = (RenderInfo*)data;  // TODO make this not another copy
//        RenderInfo* ri = ;
        Coordinates2D pt = rip->pt;
//        BlockRenderer<T>* This = (BlockRenderer<T>*) rip->This;
        int xp;
        int yp;
        double alpha, beta, gamma;
        double correctInvDepth;
        double pBaryDivisor;
        Coordinates2D ipt;
        
        T fragment;
        rip->f1 = regist(rip->fragment1);
        rip->f2 = regist(rip->fragment2);
        rip->f3 = regist(rip->fragment3);
        auto output = regist(fragment);
        //                continue;
        int CY1 = rip->CY1; //C1 + DX12 * y0 - DY12 * x0;
        int CY2 = rip->CY2; //C2 + DX23 * y0 - DY23 * x0;
        int CY3 = rip->CY3; //C3 + DX31 * y0 - DY31 * x0;
        
        for(ipt.y = pt.y; ipt.y < pt.y + rip->q; ipt.y++)
        {
            int CX1 = CY1;
            int CX2 = CY2;
            int CX3 = CY3;
            
            for(ipt.x = pt.x; ipt.x < pt.x + rip->q; ipt.x++)
            {
                if(CX1 > 0 && CX2 > 0 && CX3 > 0)
                {
                    //                            buffer[ix] = 0x0000007F;   // Blue
                    xp = ipt.x << 4;
                    yp = ipt.y << 4;
                    alpha = (double)(xp*rip->Y2 + rip->aX2Y3mX3Y2 + rip->X3*yp - xp*rip->Y3 - rip->X2*yp)*rip->A;
                    beta = (double)(rip->X1*yp + xp*rip->Y3 + rip->bX3Y1mX1Y3 - xp*rip->Y1 - rip->X3*yp)*rip->A;
                    gamma = 1.0 - alpha - beta;
                    
                    correctInvDepth = rip->invDepth1 * alpha + rip->invDepth2 * beta + rip->invDepth3 * gamma;
                    pBaryDivisor = 1.0/(alpha*rip->invW1 + beta*rip->invW2 + gamma*rip->invW3);
                    alpha *= rip->invW1 * pBaryDivisor;
                    beta  *= rip->invW2 * pBaryDivisor;
                    gamma *= rip->invW3 * pBaryDivisor;
                    
                    baryInterpolate(output, rip->f1, rip->f2, rip->f3, alpha, beta, gamma);
                    rip->p->setWithShader2( ipt, correctInvDepth, (void*)&fragment, rip->userData);
                }
                
                CX1 -= rip->FDY12;
                CX2 -= rip->FDY23;
                CX3 -= rip->FDY31;
            }
            
            CY1 += rip->FDX12;
            CY2 += rip->FDX23;
            CY3 += rip->FDX31;
            
            //                    (char*&)buffer += stride;
        }
//        pthread_exit(NULL);
            delete rip;
    }
  
public:

    RenderPipeline* p;
    RasterizerThreadPool* rtp;

    BlockRenderer(RenderPipeline* p)
    :p(p) {
    }
    
    
    void render(RenderInfo& ri) {
        RenderInfo* copy = new RenderInfo(ri);
        rtp->QueueJob(renderThread, copy);
        
    }
    void render2(RenderInfo& ri) {
        RenderInfo* copy = new RenderInfo(ri);
        rtp->QueueJob(renderThread2, copy);
    }
    
    
    void triangleFill(T* fragment1, T* fragment2, T* fragment3, void* userData) {
        //https://web.archive.org/web/20050408192410/http://sw-shader.sourceforge.net/rasterizer.html
        // 28.4 fixed-point coordinates
        
        RenderInfo ri;
        
        ri.invW1 = 1.0/fragment1->vertex.w;
        ri.invW2 = 1.0/fragment2->vertex.w;
        ri.invW3 = 1.0/fragment3->vertex.w;
        
        
        ri.Y1 = round(16.0f * fragment1->vertex.y*ri.invW1);
        ri.Y2 = round(16.0f * fragment2->vertex.y*ri.invW2);
        ri.Y3 = round(16.0f * fragment3->vertex.y*ri.invW3);
        
        ri.X1 = round(16.0f * fragment1->vertex.x*ri.invW1);
        ri.X2 = round(16.0f * fragment2->vertex.x*ri.invW2);
        ri.X3 = round(16.0f * fragment3->vertex.x*ri.invW3);
        
        
        
        
        // Deltas
        const int DX12 = ri.X1 - ri.X2;
        const int DX23 = ri.X2 - ri.X3;
        const int DX31 = ri.X3 - ri.X1;
        
        const int DY12 = ri.Y1 - ri.Y2;
        const int DY23 = ri.Y2 - ri.Y3;
        const int DY31 = ri.Y3 - ri.Y1;
        
        // Cross product check for CCW, otherwise return:
        if(DX12*DY23 - DX23*DY12 > 0)
            return;
        
        // now copy fragments:
        ri.fragment1 = *fragment1;
        ri.fragment2 = *fragment2;
        ri.fragment3 = *fragment3;
        
        ri.p = p;
        ri.userData = userData;
//        ri.This = this;
        
        ri.invDepth1 = ri.fragment1.vertex.w/ri.fragment1.vertex.z;
        ri.invDepth2 = ri.fragment2.vertex.w/ri.fragment2.vertex.z;
        ri.invDepth3 = ri.fragment3.vertex.w/ri.fragment3.vertex.z;
        // Fixed-point deltas
        ri.FDX12 = DX12 << 4;
        ri.FDX23 = DX23 << 4;
        ri.FDX31 = DX31 << 4;
        
        ri.FDY12 = DY12 << 4;
        ri.FDY23 = DY23 << 4;
        ri.FDY31 = DY31 << 4;
        
        // Bounding rectangle
        int minx = (std::min(ri.X1, std::min(ri.X2, ri.X3)) + 0xF) >> 4;
        int maxx = (std::max(ri.X1, std::max(ri.X2, ri.X3)) + 0xF) >> 4;
        int miny = (std::min(ri.Y1, std::min(ri.Y2, ri.Y3)) + 0xF) >> 4;
        int maxy = (std::max(ri.Y1, std::max(ri.Y2, ri.Y3)) + 0xF) >> 4;
        
        // Block size, standard 8x8 (must be power of two)
        ri.q = 8;
        
        // Start in corner of 8x8 block
        minx &= ~(ri.q - 1);
        miny &= ~(ri.q - 1);
        
        //barycentric parameters:
        //    int A = X1*Y2 + X2*Y3 + X3*Y1 - X1*Y3 - X2*Y1 - X3*Y2;
        //       double A = 1.0/(double)(X1*(Y2-Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2));
        ri.A = 1.0/(double)(ri.X1*(DY23) + ri.X2*(DY31) + ri.X3*(DY12));
        ri.aX2Y3mX3Y2 = ri.X2*ri.Y3 - ri.X3*ri.Y2;
        ri.bX3Y1mX1Y3 = ri.X3*ri.Y1 - ri.X1*ri.Y3;
        //    int cX1Y2mX2Y1 = X1*Y2 - X2*Y1;
//        double alpha, beta, gamma;
        
        //        (char*&)colorBuffer += miny * stride;
        //    ColorRGBA color = {255,255,0,'^'};
        //    setFrameBufferRGBA(1, 1, fbo, color);
        // Half-edge constants
        int C1 = DY12 * ri.X1 - DX12 * ri.Y1;
        int C2 = DY23 * ri.X2 - DX23 * ri.Y2;
        int C3 = DY31 * ri.X3 - DX31 * ri.Y3;
        
        // Correct for fill convention
        if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
        if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
        if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;
        
//        pt;
        Coordinates2D pt;
//        T fragment;
//        ri.f1 = regist(ri.fragment1);
//        ri.f2 = regist(ri.fragment2);
//        ri.f3 = regist(ri.fragment3);
//        auto output = regist(fragment);
        
//        int numThreads = 1;
        for(pt.y = miny; pt.y < maxy; pt.y += ri.q)
        {
            for(pt.x = minx; pt.x < maxx; pt.x += ri.q)
            {
                // Corners of block
                int x0 = pt.x << 4;
                int x1 = (pt.x + ri.q - 1) << 4;
                int y0 = pt.y << 4;
                int y1 = (pt.y + ri.q - 1) << 4;
                
                // Evaluate half-space functions
                bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
                bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
                bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
                bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
                int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
                
                bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
                bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
                bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
                bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
                int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
                
                bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
                bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
                bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
                bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
                int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
                
                // Skip block when outside an edge
                if(a == 0x0 || b == 0x0 || c == 0x0) continue;
                
                // Accept whole block when totally covered
                if(a == 0xF && b == 0xF && c == 0xF)
                {
//                    ri.This = this;
                    ri.pt = pt;
                    
                    
                    render(ri);
                }
                else   // Partially covered block
                {
                    //                continue;
                    ri.CY1 = C1 + DX12 * y0 - DY12 * x0;
                    ri.CY2 = C2 + DX23 * y0 - DY23 * x0;
                    ri.CY3 = C3 + DX31 * y0 - DY31 * x0;
                    
                    ri.pt = pt;
//                    ri.This = this;
//                    ri.CY1 = CY1;
//                    ri.CY2 = CY2;
//                    ri.CY3 = CY3;
                    render2(ri);
                }
                
                
            }
        }
//        waitThreads(0);
    }
};



//template <class T> void RenderPipeline::triangleFill(T* fragments) {
template <class T> void RenderPipeline::triangleFill(T* fragment1, T* fragment2, T* fragment3, void* userData) {
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
    
    double invW1 = 1.0/fragment1->vertex.w;
    double invW2 = 1.0/fragment2->vertex.w;
    double invW3 = 1.0/fragment3->vertex.w;
    
    const int Y1 = round(16.0f * fragment1->vertex.y*invW1);
    const int Y2 = round(16.0f * fragment2->vertex.y*invW2);
    const int Y3 = round(16.0f * fragment3->vertex.y*invW3);
    
    const int X1 = round(16.0f * fragment1->vertex.x*invW1);
    const int X2 = round(16.0f * fragment2->vertex.x*invW2);
    const int X3 = round(16.0f * fragment3->vertex.x*invW3);
    
    
    
    
    // Deltas
    const int DX12 = X1 - X2;
    const int DX23 = X2 - X3;
    const int DX31 = X3 - X1;
    
    const int DY12 = Y1 - Y2;
    const int DY23 = Y2 - Y3;
    const int DY31 = Y3 - Y1;
    
    // Cross product check for CCW, otherwise return:
    if(DX12*DY23 - DX23*DY12 > 0)
        return;
    
    
    double invDepth1 = fragment1->vertex.w/fragment1->vertex.z;
    double invDepth2 = fragment2->vertex.w/fragment2->vertex.z;
    double invDepth3 = fragment3->vertex.w/fragment3->vertex.z;
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
    
    // Block size, standard 8x8 (must be power of two)
    const int q = 8;
    
    // Start in corner of 8x8 block
    minx &= ~(q - 1);
    miny &= ~(q - 1);
    
    //barycentric parameters:
    //    int A = X1*Y2 + X2*Y3 + X3*Y1 - X1*Y3 - X2*Y1 - X3*Y2;
    //       double A = 1.0/(double)(X1*(Y2-Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2));
    double A = 1.0/(double)(X1*(DY23) + X2*(DY31) + X3*(DY12));
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
    
    //           int CY1 = C1 + DX12 * (miny << 4) - DY12 * (minx << 4);
    //           int CY2 = C2 + DX23 * (miny << 4) - DY23 * (minx << 4);
    //           int CY3 = C3 + DX31 * (miny << 4) - DY31 * (minx << 4);
    Coordinates2D pt;
    Coordinates2D ipt;
    T fragment;
    auto f1 = regist(*fragment1);
    auto f2 = regist(*fragment2);
    auto f3 = regist(*fragment3);
    auto output = regist(fragment);
    
    for(pt.y = miny; pt.y < maxy; pt.y += q)
    {
        //               int CX1 = CY1;
        //               int CX2 = CY2;
        //               int CX3 = CY3;
        
        for(pt.x = minx; pt.x < maxx; pt.x += q)
        {
            // Corners of block
            int x0 = pt.x << 4;
            int x1 = (pt.x + q - 1) << 4;
            int y0 = pt.y << 4;
            int y1 = (pt.y + q - 1) << 4;
            
            // Evaluate half-space functions
            bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
            bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
            bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
            bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
            int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
            
            bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
            bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
            bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
            bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
            int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
            
            bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
            bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
            bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
            bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
            int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
            
            // Skip block when outside an edge
            if(a == 0x0 || b == 0x0 || c == 0x0) continue;
            
            // Accept whole block when totally covered
            if(a == 0xF && b == 0xF && c == 0xF)
            {
                for(ipt.y = pt.y; ipt.y < pt.y + q; ipt.y++)
                {
                    for(ipt.x = pt.x; ipt.x < pt.x + q; ipt.x++)
                    {
                        //                               buffer[ix] = 0x00007F00;   // Green
                        int xp = ipt.x << 4;
                        int yp = ipt.y << 4;
                        alpha = (double)(xp*Y2 + aX2Y3mX3Y2 + X3*yp - xp*Y3 - X2*yp)*A;
                        beta = (double)(X1*yp + xp*Y3 + bX3Y1mX1Y3 - xp*Y1 - X3*yp)*A;
                        gamma = 1.0 - alpha - beta;
                        
                        double correctInvDepth = invDepth1 * alpha + invDepth2 * beta + invDepth3 * gamma;
                        double pBaryDivisor = 1.0/(alpha*invW1 + beta*invW2 + gamma*invW3);
                        alpha *= invW1 * pBaryDivisor;
                        beta  *= invW2 * pBaryDivisor;
                        gamma *= invW3 * pBaryDivisor;
                        
                        baryInterpolate(output, f1, f2, f3, alpha, beta, gamma);
                        setWithShader2( ipt, correctInvDepth, (void*) &fragment, userData);
                    }
                    
                    //                           (char*&)buffer += stride;
                }
            }
            else   // Partially covered block
            {
                //                continue;
                int CY1 = C1 + DX12 * y0 - DY12 * x0;
                int CY2 = C2 + DX23 * y0 - DY23 * x0;
                int CY3 = C3 + DX31 * y0 - DY31 * x0;
                
                for(ipt.y = pt.y; ipt.y < pt.y + q; ipt.y++)
                {
                    int CX1 = CY1;
                    int CX2 = CY2;
                    int CX3 = CY3;
                    
                    for(ipt.x = pt.x; ipt.x < pt.x + q; ipt.x++)
                    {
                        if(CX1 > 0 && CX2 > 0 && CX3 > 0)
                        {
                            //                            buffer[ix] = 0x0000007F;   // Blue
                            int xp = ipt.x << 4;
                            int yp = ipt.y << 4;
                            alpha = (double)(xp*Y2 + aX2Y3mX3Y2 + X3*yp - xp*Y3 - X2*yp)*A;
                            beta = (double)(X1*yp + xp*Y3 + bX3Y1mX1Y3 - xp*Y1 - X3*yp)*A;
                            gamma = 1.0 - alpha - beta;
                            
                            double correctInvDepth = invDepth1 * alpha + invDepth2 * beta + invDepth3 * gamma;
                            double pBaryDivisor = 1.0/(alpha*invW1 + beta*invW2 + gamma*invW3);
                            alpha *= invW1 * pBaryDivisor;
                            beta  *= invW2 * pBaryDivisor;
                            gamma *= invW3 * pBaryDivisor;
                            
                            baryInterpolate(output, f1, f2, f3, alpha, beta, gamma);
                            setWithShader2( ipt, correctInvDepth, (void*) &fragment, userData);
                        }
                        
                        CX1 -= FDY12;
                        CX2 -= FDY23;
                        CX3 -= FDY31;
                    }
                    
                    CY1 += FDX12;
                    CY2 += FDX23;
                    CY3 += FDX31;
                    
                    //                    (char*&)buffer += stride;
                }
            }
            
            
        }
    }
}



#define CLIP_PLANE_W_2 (0.0001)
template <class T> void RenderPipeline::clipW(T* input, int inputCount, T* output, int& outputCountResult) {
//    output->numVertices = 0;
    outputCountResult = 0;
    
    double factor;
    
//    char previousDot;
//    char currentDot;
    bool previousDot;
    bool currentDot;
    
//    Coordinates4D *current = &input.vertices[0];
    T* current = &input[0];
//    Coordinates4D *prior = &input.vertices[input.numVertices-1];
    T* prior = &input[inputCount-1];
//    previousDot = prior->w < CLIP_PLANE_W_2 ? -1 : 1;
    previousDot = prior->vertex.w > CLIP_PLANE_W_2;
//    for (; current != &input.vertices[input.numVertices]; ) {
    for (; current != &input[inputCount]; ) {
//        currentDot = (current->w < CLIP_PLANE_W) ? -1 : 1;
        currentDot = current->vertex.w > CLIP_PLANE_W_2;
        
//        if (previousDot * currentDot < 0) {
        if (previousDot != currentDot) {
//            factor = (CLIP_PLANE_W_2 - prior->w) / (prior->w - current->w);
            factor = (CLIP_PLANE_W_2 - prior->vertex.w) / (current->vertex.w - prior->vertex.w );
//            Coordinates4D diff = vectorSubtract(*current, *prior);
            Coordinates4D diff = vectorSubtract(current->vertex, prior->vertex);
            diff.x *= factor;
            diff.y *= factor;
            diff.z *= factor;
            diff.w *= factor;
//            diff = vectorAdd(*prior, diff);
            diff = vectorAdd(prior->vertex, diff);
            
//            output->vertices[output->numVertices] = diff;
//            output->numVertices++;
            
//            double alpha = factor;
//            double beta = 1.0 - alpha;
//            double invW1 = current->vertex.w/current->vertex.z;
//            double invW2 = prior->vertex.w/prior->vertex.z;
//            double currentDepth = 1.0/(alpha*invW1 + beta*invW2);
////            double pBaryDivisor = 1.0/(alpha*invW1 + beta*invW2);
//            alpha *= invW1 * currentDepth;
//            beta  *= invW2 * currentDepth;
            
            
            auto f1 = regist(*current);
            auto f2 = regist(*prior);
            auto o1 = regist(output[outputCountResult]);
            
//            if(alpha > 0 && beta > 0 && alpha < 1 && beta < 1)
                baryInterpolate2d(o1, f1, f2, factor, 1.0-factor);
//            else
//                output[outputCountResult] = *current;   // HACK just to get attribute to be non-0
            output[outputCountResult].vertex = diff;
            outputCountResult++;
        }
        
//        if (currentDot > 0) {
        if (currentDot ) {
//            output->vertices[output->numVertices] = *current;
//            output->numVertices++;
            output[outputCountResult] = *current;
            outputCountResult++;
        } //else {
          //  this->mvaddstring(0, 2, "W clipped");
        //}
        
        previousDot = currentDot;
        prior = current;
        current++;
    }
}

template <class T> void RenderPipeline::clipPlane(T* input, int inputCount, T* output, int& outputCountResult, int axis, int plane) {
//    output->numVertices = 0;
    outputCountResult = 0;
    
    double factor;
    
//    char previousDot;
//    char currentDot;
//    int previousDot;
//    int currentDot;
    bool previousDot;
    bool currentDot;
    
//    Coordinates4D *current = &input.vertices[0];
//    Coordinates3D *currentNormal = &input.normals[0];
//    ColorRGBA *currentColor = &input.colors[0];
    T *current = &input[0];
    
//    Coordinates4D *prior = &input.vertices[input.numVertices-1];
    T *prior = &input[inputCount-1];
//    previousDot = ((double*)prior)[axis] <= prior->w ;
    if(plane == 1) {
        previousDot = ((double*)&prior->vertex)[axis] <= prior->vertex.w ;
    } else {
        previousDot = -((double*)&prior->vertex)[axis] <= prior->vertex.w;
    }
//    for (; current != &input.vertices[input.numVertices]; ) {
    for (; current != &input[inputCount]; ) {
//        currentDot = ((double*)current)[axis] <= current->w;
        if(plane == 1) {
            currentDot = ((double*)&current->vertex)[axis] <= current->vertex.w;
        } else {
            currentDot = -((double*)&current->vertex)[axis] <= current->vertex.w;
        }
//            factor = current->w/current->x;
        // duplicate and shift
        
        if (previousDot != currentDot ) {
//            mvprintw(line++, 0, " - - Clip against plane %d", axis);
//            factor = (prior->w - ((double*)prior)[axis])/((prior->w - ((double*)prior)[axis]) - (current->w - ((double*)current)[axis])) ;
            if(plane == 1 ) {
                factor = (prior->vertex.w - ((double*)&prior->vertex)[axis])/((prior->vertex.w - ((double*)&prior->vertex)[axis]) - (current->vertex.w - ((double*)&current->vertex)[axis])) ;
            } else {
                factor = (prior->vertex.w + ((double*)&prior->vertex)[axis])/((prior->vertex.w + ((double*)&prior->vertex)[axis]) - (current->vertex.w + ((double*)&current->vertex)[axis])) ;
            }
//            Coordinates4D diff = vectorSubtract(*current, *prior);
            Coordinates4D diff = vectorSubtract(current->vertex, prior->vertex);
            diff.x *= factor;
            diff.y *= factor;
            diff.z *= factor;
            diff.w *= factor;
//            if(axis == 2 && plane == -1)
//                diff.w *= -1;
            diff = vectorAdd(prior->vertex, diff);
            
//            output->vertices[output->numVertices] = diff;
//            output->normals[output->numVertices] = *currentNormal; // this should be interpolated as well
//            output->colors[output->numVertices] = *currentColor; // this should be interpolated as well
//            output->numVertices++;
//            double alpha = factor;
//            double beta = 1.0 - alpha;
//            double invW1 = current->vertex.w/current->vertex.z;
//            double invW2 = prior->vertex.w/prior->vertex.z;
//            double currentDepth = 1.0/(alpha*invW1 + beta*invW2);
//            double pBaryDivisor = 1.0/(alpha*invW1 + beta*invW2);
//            alpha *= invW1 * currentDepth;
//            beta  *= invW2 * currentDepth;
            
            
            auto f1 = regist(*current);
            auto f2 = regist(*prior);
            auto o1 = regist(output[outputCountResult]);
            
//            baryInterpolate2d(o1, f1, f2, alpha, beta);
            baryInterpolate2d(o1, f1, f2, factor, 1.0-factor);
//            output[outputCountResult] = *current;   // HACK just to get attribute to be non-0
            output[outputCountResult].vertex = diff;
            // TODO: interpolate rest of output[outputCountResult] using factor
            outputCountResult++;
        }
        
        if (currentDot ) {
//                mvprintw(line++, 0, " - - Vertex Insertion");
//            output->vertices[output->numVertices] = *current;
//            output->normals[output->numVertices] = *currentNormal;
//            output->colors[output->numVertices] = *currentColor;
//            output->numVertices++;
            output[outputCountResult] = *current;
            outputCountResult++;
        } //else {
          //  this->mvaddstring(axis, 3+axis, "Plane clipped");
        //}
        
        previousDot = currentDot;
        prior = current;
        current++;
//        currentNormal++;
//        currentColor++;
    }
    return;
}

template <class T> void RenderPipeline::clipPolygon(T* input, int inputCount, T* output, int& outputCountResult) {
    T copy[20];
//    mvprintw(line++, 0, "Polygon %d");
//    for (int i = 0; i < input.numVertices; i++) {
//        if (input.vertices[i].w <= 0) {
//            mvprintw(line++, 0, "Polygon W clip needed: %f", input.vertices[i].w);
//            return 0;
//        }
//    }
    
    // The following functions are from https://fabiensanglard.net/polygon_codec/
    
//    clipW(input, output, line);    // I don't think this is necessary (well it is, but it's an ultra rare corner case?)
    
    
// For debugging:
    outputCountResult = inputCount;
    for(int i = 0; i < inputCount; i++)
        copy[i] = input[i];
    
    
    
//    return;
    
    clipW(copy, outputCountResult, output, outputCountResult);
    
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 0, 1);
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 0, -1);
    
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 1, 1);
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 1, -1);
    
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 2, 1);
    
    for(int i = 0; i < outputCountResult; i++)
        copy[i] = output[i];
    clipPlane(copy, outputCountResult, output, outputCountResult, 2, -1);
    
   // return 1;

}

#endif
