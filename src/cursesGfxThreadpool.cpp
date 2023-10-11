#include "curses-gfx-handler.h"

#include <unistd.h>


void RasterizerThreadPool::ThreadLoop() {
    while (true) {
        std::pair<std::function<void(void*)>,void*> job;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            mutex_condition.wait(lock, [this] {
                return !jobs.empty() || should_terminate;
            });
            if (should_terminate) {
                return;
            }
            
            job = jobs.front();
            jobs.pop();
        }
        job.first(job.second);
        {
            std::unique_lock<std::mutex> lock(busy_mutex);
            activeThreadCount--;
        }
        mutex_cv_busy.notify_all();
    }
}

void RasterizerThreadPool::QueueJob(const std::function<void(void*)>& job, void* data) {

    {
        std::unique_lock<std::mutex> lock(busy_mutex);
        activeThreadCount++;
    }
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        jobs.push(std::pair(job, data));
    }
    mutex_condition.notify_one();
}


bool RasterizerThreadPool::busy() {
    bool poolbusy;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        poolbusy = activeThreadCount != 0;//!jobs.empty();
    }
    return poolbusy;
}

void RasterizerThreadPool::busyWait() {
    {
        std::unique_lock<std::mutex> lock(busy_mutex);
            mutex_cv_busy.wait(lock, [this] {
                return activeThreadCount == 0;
            });
    }
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
