#include "ThreadPool.h"

Worker::Worker(ThreadPool &threadPool) : threadPool(threadPool) {
	thread = std::thread(&Worker::threadLoop, this);
}

Worker::~Worker() {
	{
		std::lock_guard<std::mutex> lock(sleepMux);
		shouldExit = true;
	}

	sleepCV.notify_all();

	thread.join();
}

void Worker::threadLoop() {
	while (true) {
		std::unique_lock<std::mutex> lock(sleepMux);
		sleepCV.wait(lock, [&] { return hasWork || shouldExit; });
		lock.unlock();
		if (shouldExit) { return; }

		work();

		hasWork = false;
		threadPool.threadDone();
	}
}

void Worker::start(std::function<void()> work) {
	{
		std::lock_guard<std::mutex> lock(sleepMux);
		this->work = std::move(work);
		hasWork = true;
	}

	sleepCV.notify_all();
}

void ThreadPool::init(int numWorkers) {
	this->numWorkers = numWorkers;

	for (int i = 0; i < numWorkers; i++) {
		workers.emplace_back(std::make_unique<Worker>(*this));
	}

	std::cout << numWorkers << " threads started\n";
}

void ThreadPool::computeOnThreads(int size, const std::function<void(int, int)> &work) {
	if (numWorkers == 1) {
		workers[0]->start([=] { work(1, size-1); });
	} else {
		int sizePerThread = (size - 2) / (numWorkers - 1);
		int remainingWork = (size - 2) % (numWorkers - 1);
		for (int i = 0; i < numWorkers - 1; i++) {
			workers[i]->start([=] { work(i * sizePerThread + 1, (i + 1) * sizePerThread + 1); });
		}
		workers[numWorkers - 1]->start([=] { work((numWorkers - 1) * sizePerThread + 1, (numWorkers - 1) * sizePerThread + 1 + remainingWork); });
	}

	std::unique_lock<std::mutex> lock(doneMux);
	doneCV.wait(lock, [&] { return numWorkersDone == numWorkers; });
	numWorkersDone = 0;
}

void ThreadPool::threadDone() {
	{
		std::lock_guard<std::mutex> doneLock(doneMux);
		numWorkersDone++;
	}

	doneCV.notify_all();
}
