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


//template <class T> void clipPolygon(T* input, int inputCount, T* output, int& outputCountResult);

class RenderPipeline {
public:
	FrameBuffer* fbo;
	FrameBuffer* depthBuffer;
	
	double depthClearColor;
	
	std::queue<pthread_t*> renderThreads;
	
//	DepthBuffer tempD;
    void* userData;
    
    Mat4D viewport;
    
//    void (*vertexShader)(...);
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
    template <class T, class U> void rasterizeShader(T* vertexInfo, U* uniformInfo, int triangleLayout[][3], int numTriangles, void* userData, void (*vertexShader)(U*, T&, T&));
	void rasterizePolygonsShader(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	void rasterizeThreaded(Polygon4D* polygons, int count, Mat4D& modelView, Mat4D& projection, Mat4D& viewport, void* userData, int &line);
	
	void setRenderBuffer(const int& x, const int& y, ColorRGBA& color);
	void setRenderBuffer(const int& index, ColorRGBA& color);
	
	void renderBufferToTerminal();
	void depthBufferToTerminal();
    
    template <class T> void clipPolygon(T* input, int inputCount, T* output, int& outputCountResult);
    template <class T> void clipW(T* input, int inputCount, T* output, int& outputCountResult);
    template <class T> void clipPlane(T* input, int inputCount, T* output, int& outputCountResult, int axis, int plane);
    
    
    void mvaddstring(int x, int y, const char* string);
};



template <class T, class U> void RenderPipeline::rasterizeShader(T* vertexInfo, U* uniformInfo, int triangleLayout[][3], int numTriangles, void* userData, void (*vertexShader)(U*, T&, T&)) {
    T scratch[3];
    int scratchLayout[][3] = {
        {0, 1, 2}
    };
    T scratchClipped[10];
    int clippedVertexCount;
    
    this->userData = userData;
    
    for (int t = 0; t < numTriangles; t++) {
        scratch[0] = vertexInfo[triangleLayout[t][0]];
        scratch[1] = vertexInfo[triangleLayout[t][1]];
        scratch[2] = vertexInfo[triangleLayout[t][2]];
        
        for(int i = 0; i < 3; i++) {
            //            scratch[i].vertex = matrixVectorMultiply(modelView, scratch[i].vertex);
            //            scratch[i].location = scratch[i].vertex;
            //            scratch[i].normal = matrixVectorMultiply(modelView, scratch[i].normal);
            //            scratch[i].vertex = matrixVectorMultiply(projection, scratch[i].vertex);
            //            vertexShader = (void*)myVertexShader;
            //            myVertexShader(uniformInfo, scratch[i], vertexInfo[triangleLayout[t][i]]);
            vertexShader(uniformInfo, scratch[i], vertexInfo[triangleLayout[t][i]]);
            //            scratch[i].vertex.y = -scratch[i].vertex.y;
            
            // For every vertex, perform a vertex shading:
            // Should be able to take inputs from a set of custom VertexInformation like model, projection, view matrics
            // Give to a user-defined shader which produces user-defined struct of vertex attributes
        }
        
        // Then.. clip geometry?
        clipPolygon(scratch, 3, scratchClipped, clippedVertexCount);
            
        for(int i = 0; i < clippedVertexCount; i++) {
            // Then.. apply a viewport and provide to triangle rasterizer
//            scratch[i].vertex = matrixVectorMultiply(viewport, scratch[i].vertex);
            scratchClipped[i].vertex = matrixVectorMultiply(viewport, scratchClipped[i].vertex);
        }
        
        
//            triangleFill(&scratch[0], &scratch[1], &scratch[2]);
        for(int i = 2; i < clippedVertexCount; i++) {
            triangleFill(&scratchClipped[0], &scratchClipped[i-1], &scratchClipped[i]);
        }
//        fbo->data[20*4 + 3] = '0'+clippedVertexCount;

    }
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
    
    double invW1 = 1.0/fragment1->vertex.w;
    double invW2 = 1.0/fragment2->vertex.w;
    double invW3 = 1.0/fragment3->vertex.w;
    
    const int Y1 = round(16.0f * fragment1->vertex.y*invW1);
        const int Y2 = round(16.0f * fragment2->vertex.y*invW2);
        const int Y3 = round(16.0f * fragment3->vertex.y*invW3);

    const int X1 = round(16.0f * fragment1->vertex.x*invW1);
    const int X2 = round(16.0f * fragment2->vertex.x*invW2);
    const int X3 = round(16.0f * fragment3->vertex.x*invW3);
    
    
//    double invDepth1 = 1.0/fragment1->vertex.z;
//    double invDepth2 = 1.0/fragment2->vertex.z;
//    double invDepth3 = 1.0/fragment3->vertex.z;
    double invDepth1 = fragment1->vertex.w/fragment1->vertex.z;
    double invDepth2 = fragment2->vertex.w/fragment2->vertex.z;
    double invDepth3 = fragment3->vertex.w/fragment3->vertex.z;

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
                       
                       
                       setWithShader2( pt, correctInvDepth, userData, (void*) &fragment);
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
            if(plane == 1) {
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
//    for(int i = 0; i < outputCountResult; i++)
//        copy[i] = output[i];
//    clipPlane(copy, outputCountResult, output, outputCountResult, 2, -1);
    
    
   // return 1;

}

#endif
