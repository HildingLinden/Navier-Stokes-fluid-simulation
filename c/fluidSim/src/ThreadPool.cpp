#include <immintrin.h>
#include "ThreadPool.h"

std::condition_variable doneCV;
std::mutex doneMux;
int numWorkersDone;

void Worker::threadLoop()
{
	while (true) {
		std::unique_lock<std::mutex> lock(sleepMux);
		sleepCV.wait(lock, [&] { return hasWork || shouldExit; });
		lock.unlock();
		if (shouldExit) return;
		
		work();

		{
			std::lock_guard<std::mutex> doneLock(doneMux);
			hasWork = false;
			numWorkersDone++;
		}

		doneCV.notify_all();
	}
}

void Worker::start(std::function<void()> work)
{
	{
		std::lock_guard<std::mutex> lock(sleepMux);
		this->work = work;
		hasWork = true;
	}

	sleepCV.notify_all();
}

void Worker::terminate()
{
	{
		std::lock_guard<std::mutex> lock(sleepMux);
		shouldExit = true;
	}

	sleepCV.notify_all();
}

void ThreadPool::init(int numWorkers)
{
	this->numWorkers = numWorkers;
	
	workers = new Worker[numWorkers];
	for (int i = 0; i < numWorkers; i++) {
		workers[i].thread = std::thread(&Worker::threadLoop, &workers[i]);
	}

	std::cout << numWorkers << " threads started\n";
}

void ThreadPool::computeOnThreads(int size, std::function<void(int, int)> work)
{
	if (numWorkers == 1) {
		workers[0].start([=] { work(1, size-1); });
	}
	else {
		int sizePerThread = (size - 2) / (numWorkers - 1);
		int remainingWork = (size - 2) % (numWorkers - 1);
		for (int i = 0; i < numWorkers - 1; i++) {
			workers[i].start([=] { work(i * sizePerThread + 1, (i + 1) * sizePerThread + 1); });
		}
		workers[numWorkers - 1].start([=] { work((numWorkers - 1) * sizePerThread + 1, (numWorkers - 1) * sizePerThread + 1 + remainingWork); });
	}

	std::unique_lock<std::mutex> lock(doneMux);
	doneCV.wait(lock, [&] { return numWorkersDone == numWorkers; });
	numWorkersDone = 0;
}

ThreadPool::~ThreadPool()
{
	for (int i = 0; i < numWorkers; i++) {
		workers[i].terminate();
	}
	for (int i = 0; i < numWorkers; i++) {
		workers[i].thread.join();
	}
	std::cout << numWorkers << " threads terminated\n";
}
