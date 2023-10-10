#include <unistd.h>
#include <cstdio>

#include "curses-gfx-handler.h"


void function(void* data) {
    double time = *(double*) data;
    printf("Thread started %f\n", time);
    usleep(time*1000000);
    printf("Thread done with time %f\n", time);
}

int main(void) {
    RasterizerThreadPool mRasterizerThreadPool;
    
    mRasterizerThreadPool.Start(3);
    usleep(100000);
    double time = 1.34;
    mRasterizerThreadPool.QueueJob(function, &time);
    double time2 = 2.34;
    mRasterizerThreadPool.QueueJob(function, &time2);
    printf("Threads queueueued\n");
//    while(mRasterizerThreadPool.busy()) {
//        printf(" - busy(): %d\n",mRasterizerThreadPool.busy());
//    usleep(100000);
//}

    mRasterizerThreadPool.busyWait();
    printf("busyWait Finished\n");
    mRasterizerThreadPool.Stop();
    printf("Complete.\n");
    return 0;
}
