#pragma once
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

class ThreadPool;

class Worker {
	bool m_shouldExit = false;
	bool m_hasWork = false;
	std::condition_variable m_sleepCV;
	std::mutex m_sleepMux;
	std::function<void()> m_work;
	ThreadPool &m_threadPool;
	std::thread m_thread;

public:
	explicit Worker(ThreadPool &threadPool);
	~Worker();										// Rule of 5
	Worker(const Worker &) = delete;				// Copy constructor
	Worker(Worker &&) = delete;						// Move constructor
	Worker &operator=(const Worker &) = delete;		// Copy assignment
	Worker &operator=(Worker &&) = delete;			// Move assignment
	void threadLoop();
	void start(std::function<void()> work);
};

class ThreadPool {
	std::vector<std::unique_ptr<Worker>> m_workers;
	int m_numWorkers = 0;
	std::condition_variable m_doneCV;
	std::mutex m_doneMux;
	int m_numWorkersDone = 0;

public:
	void init(int numWorkers);
	void computeOnThreads(int size, const std::function<void(int, int)> &work);
	void threadDone();
};
