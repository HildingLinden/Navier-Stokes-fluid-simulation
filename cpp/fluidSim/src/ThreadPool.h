#pragma once
#include <mutex>
#include <condition_variable>
#include <functional>
#include <thread>
#include <iostream>
#include <vector>

class ThreadPool;

class Worker {
	bool shouldExit = false;
	bool hasWork = false;
	std::condition_variable sleepCV;
	std::mutex sleepMux;
	std::function<void()> work;
	ThreadPool &threadPool;
	std::thread thread;

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
	std::vector<std::unique_ptr<Worker>> workers;
	int numWorkers = 0;
	std::condition_variable doneCV;
	std::mutex doneMux;
	int numWorkersDone = 0;

public:
	void init(int numWorkers);
	void computeOnThreads(int size, const std::function<void(int, int)> &work);
	void threadDone();
};
