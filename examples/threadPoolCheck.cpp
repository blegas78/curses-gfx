#include <unistd.h>
#include <cstdio>
#include <chrono>

#include <ncurses.h>

#include "curses-gfx.h"
#include "curses-gfx-handler.h"


struct RenderData {
    int x, y;
    int q;
    unsigned int cycles;
};

std::mutex myMutex;

void function(void* data) {
    RenderData d = *(RenderData*) data;
//    printf("Thread started %f\n", d.time);
    for(int x = 0; x < d.q; x++) {
        for(int y = 0; y < d.q; y++) {
//            usleep(d.time*1000000);
            
            volatile unsigned int i = 0, j = 1;
            while(i++ != j++) {
                if(j == d.cycles) {
                    break;
                }
            }
//            char c = 'A' + 26*(0.5+0.5*sin(cos(y+d.y*log(pow(x+d.x,3.78642)))));
            
            std::unique_lock<std::mutex> lock(myMutex);
            mvaddch(y+d.y, x+d.x, '*');
        }
    }
    
}

int main(void) {
    double estimatedSingleThreadTime = 0;
    CursesGfxTerminal mCursesGfxTerminal;
    mCursesGfxTerminal.setupTerminal();
    RasterizerThreadPool mRasterizerThreadPool;
    const int numThreads = 8;
    const int numCols = 12;
    const int numRows = 5;
    const int q = 16;
    RenderData data[numCols*numRows];
    for(int r = 0; r < numRows; r++) {
        for(int i = 0; i < numCols; i++) {
            data[i + r*numCols].x = i*(q+1) + 5;
            data[i + r*numCols].y = r*(q+1) + 5;
            
            data[i + r*numCols].cycles = 3000000 * (1 + 0.75*sin(r*i));
            data[i + r*numCols].q = q;
        }
    }
    
    mRasterizerThreadPool.Start(numThreads);

    auto before = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < numCols*numRows; i++) {
        mRasterizerThreadPool.QueueJob(function, &data[i]);
    }
    while( mRasterizerThreadPool.busy()) {
        usleep(1000);
        std::unique_lock<std::mutex> lock(myMutex);
        refresh();
    }
    mRasterizerThreadPool.busyWait();
    auto now = std::chrono::high_resolution_clock::now();
//    printf("busyWait Finished\n");
    mRasterizerThreadPool.Stop();
//    printf("Complete.\n");
    
    
    std::chrono::duration<double, std::milli> float_ms = now - before;
    usleep(500000);
    
    // check single "thread" performance:
    before = std::chrono::high_resolution_clock::now();
    function(&data[0]);
    now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> single_float_ms = now - before;
    
    mCursesGfxTerminal.cleanupTerminal();
    printf("Time: %f Single thread estimate: %f Speed increase ratio: %f\n", float_ms.count()/1000.0, single_float_ms.count()/1000.0*numRows*numCols, numRows*numCols*single_float_ms.count()/float_ms.count());
    return 0;
}
