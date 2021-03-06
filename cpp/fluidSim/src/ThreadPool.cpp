#include "ThreadPool.h"

Worker::Worker(ThreadPool &threadPool) : m_threadPool(threadPool) {
	m_thread = std::thread(&Worker::threadLoop, this);
}

Worker::~Worker() {
	{
		std::lock_guard<std::mutex> lock(m_sleepMux);
		m_shouldExit = true;
	}

	m_sleepCV.notify_all();

	m_thread.join();
}

void Worker::threadLoop() {
	while (true) {
		std::unique_lock<std::mutex> lock(m_sleepMux);
		m_sleepCV.wait(lock, [&] { return m_hasWork || m_shouldExit; });
		lock.unlock();
		if (m_shouldExit) { return; }

		m_work();

		m_hasWork = false;
		m_threadPool.threadDone();
	}
}

void Worker::start(std::function<void()> work) {
	{
		std::lock_guard<std::mutex> lock(m_sleepMux);
		m_work = std::move(work);
		m_hasWork = true;
	}

	m_sleepCV.notify_all();
}

void ThreadPool::init(int numWorkers) {
	m_numWorkers = numWorkers;

	for (int i = 0; i < numWorkers; i++) {
		m_workers.emplace_back(std::make_unique<Worker>(*this));
	}

	std::cout << numWorkers << " threads started\n";
}

void ThreadPool::computeOnThreads(int size, const std::function<void(int, int)> &work) {
	if (m_numWorkers == 1) {
		m_workers[0]->start([=] { work(1, size-1); });
	} else {
		const int sizePerThread = (size - 2) / (m_numWorkers - 1);
		const int remainingWork = (size - 2) % (m_numWorkers - 1);
		for (int i = 0; i < m_numWorkers - 1; i++) {
			m_workers[i]->start([=] { work(i * sizePerThread + 1, (i + 1) * sizePerThread + 1); });
		}
		m_workers[m_numWorkers - 1]->start([=] { work((m_numWorkers - 1) * sizePerThread + 1, (m_numWorkers - 1) * sizePerThread + 1 + remainingWork); });
	}

	std::unique_lock<std::mutex> lock(m_doneMux);
	m_doneCV.wait(lock, [&] { return m_numWorkersDone == m_numWorkers; });
	m_numWorkersDone = 0;
}

void ThreadPool::threadDone() {
	{
		std::lock_guard<std::mutex> doneLock(m_doneMux);
		m_numWorkersDone++;
	}

	m_doneCV.notify_all();
}
