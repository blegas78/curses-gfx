#include "curses-gfx-handler.h"

#include <unistd.h>

//bool RasterizerThreadPool::should_terminate = false;
//std::mutex RasterizerThreadPool::queue_mutex;                  // Prevents data races to the job queue
//std::condition_variable RasterizerThreadPool::mutex_condition; // Allows threads to wait on new jobs or termination
//std::vector<std::thread> RasterizerThreadPool::threads;
//std::queue<std::pair<std::function<void(void*)>,void*>> RasterizerThreadPool::jobs;

void RasterizerThreadPool::ThreadLoop() {
    while (true) {
//        std::function<void()> job;
//        printf("in a threadLoop\n");
        std::pair<std::function<void(void*)>,void*> job;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            mutex_condition.wait(lock, [this] {
                return !jobs.empty() || should_terminate;
            });
            if (should_terminate) {
                return;
            }
            
//            printf("threadloop got a popper\n");
            job = jobs.front();
            jobs.pop();
//            {
//                std::unique_lock<std::mutex> lock(busy_mutex);
//                activeThreadCount--;
//            }
        }
        job.first(job.second);
        {
            std::unique_lock<std::mutex> lock(busy_mutex);
            activeThreadCount--;
            //if(activeThreadCount == 0) {
            //}
        }
        mutex_cv_busy.notify_all();
//        if(activeThreadCount == 0 && jobs.empty()) {
//            mutex_cv_busy.notify_all();
//        }
    }
}

//void RasterizerThreadPool::QueueJob(const std::function<void()>& job) {
void RasterizerThreadPool::QueueJob(const std::function<void(void*)>& job, void* data) {
//    while(busy()) {};
    {
        std::unique_lock<std::mutex> lock(busy_mutex);
        activeThreadCount++;
    }
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
//        jobs.push(job);
        jobs.push(std::pair(job, data));
    }
    mutex_condition.notify_one();
}


bool RasterizerThreadPool::busy() {
    bool poolbusy;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        poolbusy = !jobs.empty();
    }
    return poolbusy;
}

void RasterizerThreadPool::busyWait() {
//    while(activeThreadCount != 0) {
//        usleep(10);
//    }
//    if(activeThreadCount != 0)
    {
        std::unique_lock<std::mutex> lock(busy_mutex);
//        if(activeThreadCount != 0) {
            mutex_cv_busy.wait(lock, [this] {
                return activeThreadCount == 0;
            });
//        }
    }
    
//    busy();
}

void RasterizerThreadPool::Start(uint32_t numThreads) {
    activeThreadCount = 0;
    should_terminate = false;
    const uint32_t num_threads = numThreads;//std::thread::hardware_concurrency(); // Max # of threads the system supports
    for (uint32_t ii = 0; ii < num_threads; ++ii) {
        threads.emplace_back(std::thread(&RasterizerThreadPool::ThreadLoop,this));
    }
}

void RasterizerThreadPool::Stop() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        should_terminate = true;
    }
    mutex_condition.notify_all();
    for (std::thread& active_thread : threads) {
        active_thread.join();
    }
    threads.clear();
}
